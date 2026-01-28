#ifndef OPTIMIZATION_ALGORITHM_HPP
#define OPTIMIZATION_ALGORITHM_HPP

#include "placement_algorithm.hpp"
#include "constrained_dp_algorithm.hpp"

namespace optimization_algorithm {

template<class T> concept COLORABLE = stepwise_colorable::STEPWISE_COLORABLE<T>;

using std::size_t;
using std::vector;
using std::shared_ptr;
using std::unique_ptr;
using std::any;

template<placement_algorithm::PLACEMENT_ALGORITHM PlacementAlgorithm> class RecursivePlacement : public common::LogInfo{
public:
	RecursivePlacement(int verbose = common::LogInfo::DEFAULT_VERBOSE) : LogInfo(verbose){}

	common::AnnotatedBinaryTree& operator()(typename PlacementAlgorithm::Color& color, typename PlacementAlgorithm::ThreadPool& threadPool, common::AnnotatedBinaryTree& tree, std::ranges::range auto const& taxa) {
		size_t nLeaves = tree.leaves().size();
		LogInfo const& vlog = log.vlog<1>();
		for (size_t iTaxon : taxa) {
			if (nLeaves < 3) tree.emplaceRoot(common::AnnotatedBinaryTree::LEAF_ID, iTaxon);
			else {
				log << "Placing " << nLeaves + 1 << "-th taxon (" << common::taxonName2ID[iTaxon] << ")..." << std::endl;
				tree.displaySimpleNewick<double, double>(vlog << "Current tree: ", "", "") << std::endl;
				PlacementAlgorithm placement(color, threadPool, tree, verbose + 1);
				placement.place(iTaxon);
			}
			nLeaves++;
		}
		return tree;
	}
};

class ProcedureBase {
public:
	class Subprocedure {
	protected:
		bool declarationDone = false;
		bool runDone = false;
		any returned;

		virtual void declare() = 0;
		virtual void run() = 0;

		vector<shared_ptr<Subprocedure> > declarationDependencies, declarationDependenciesOnRun, runDependencies;
		vector<shared_ptr<Subprocedure> >& dd = declarationDependencies;
		vector<shared_ptr<Subprocedure> >& dr = declarationDependenciesOnRun;
		vector<shared_ptr<Subprocedure> >& rr = runDependencies;

	public:
		void recursivelyDeclare() {
			if (declarationDone) return;
			declarationDone = true;
			for (shared_ptr<Subprocedure>& dep : declarationDependencies) {
				dep->recursivelyDeclare();
			}
			for (shared_ptr<Subprocedure>& dep : declarationDependenciesOnRun) {
				dep->recursivelyDeclare();
				dep->recursivelyRun();
			}
			declare();
		}

		void recursivelyRun() {
			recursivelyDeclare();
			if (runDone) return;
			runDone = true;
			for (shared_ptr<Subprocedure>& dep : runDependencies) {
				dep->recursivelyRun();
			}
			run();
		}

		template<class T> void prepareReturn(T&& t) {
			returned = std::make_shared<T>(std::forward<T>(t));
		}

		template<class T> shared_ptr<T> returnValue() try{
			recursivelyRun();
			return std::any_cast<shared_ptr<T> >(returned);
		}
		catch (std::bad_any_cast const& e) {
			throw std::runtime_error("Error: Subprocedure return value type mismatch.");
		}
	};

	using S = shared_ptr<Subprocedure>;
	
	template<class T> static T& get(S& s) { return *(s->returnValue<T>()); } 

	template<typename T> class Input : public virtual Subprocedure {
	public:
		Input(shared_ptr<T> t) {
			returned = t;
		}

		Input(T* t) {
			shared_ptr<T> p(t);
			returned = p;
		}

		Input(T& t) : Input(&t) {}

		Input(T t) : Input(new T(std::move(t))) {}

	protected:
		virtual void declare() override {}

		virtual void run() override {}
	};

	template<typename T, typename... Args> static shared_ptr<Input<T> > pInput(Args... args) {
		shared_ptr<Input<T> > p(new Input<T>(new T(std::forward<Args>(args)...)));
		return p;
	}
};

template<COLORABLE C> struct DefaultProcedureAttributes {
	using Color = C;
	using score_t = Color::score_t;
	template<typename T> using Scheduler = thread_pool::SimpleScheduler<T>;
	using Placement = placement_algorithm::StepwiseColorPlacement<placement_algorithm::StepwiseColorPlacementDefaultAttributes<Color> >;
	using Random = common::Random<std::mt19937_64>;
};

template<typename Attributes> class Procedure : public ProcedureBase {
	using Log = common::LogInfo;
	using score_t = Attributes::score_t;
	using Tree = common::AnnotatedBinaryTree;
	using Random = Attributes::Random;
	using Color = Attributes::Color;
	template<typename T> using Scheduler = thread_pool::SimpleScheduler<T>;
	using ThreadPool = thread_pool::ThreadPool<score_t, Scheduler>;
	using Placement = Attributes::Placement;
	using Recursive = RecursivePlacement<Placement>;
	using ConstrainedDP = constrained_dp_algorithm::ConstrainedDP<score_t, Random>;
	static score_t constexpr ZERO = Placement::ZERO;
	
public:
	class RecursivePlacementSubprocedure : public virtual Subprocedure {
	public:
		RecursivePlacementSubprocedure(S color, S tp, S tree, S taxa, S verbose) {
			dr = { verbose };
			rr = { color, tp, tree, taxa };
		}

	protected:
		unique_ptr<Recursive> rp;
		
		virtual void declare() override {
			auto& verbose = get<int>(dr[0]);
			rp.reset(new Recursive(verbose));
		}

		virtual void run() override {
			auto& color = get<Color>(rr[0]);
			auto& tp = get<ThreadPool>(rr[1]);
			auto& tree = get<Tree>(rr[2]);
			auto& taxa = get<vector<size_t> >(rr[3]);
			(*rp)(color, tp, tree, taxa);
			returned = rr[2]->returnValue<Tree>();
		}
	};

	class ScoreTreeSubprocedure : public virtual Subprocedure {
	public:
		ScoreTreeSubprocedure(S color, S tp, S tree, S verbose) {
			dr = { color, tp, tree, verbose };
		}

	protected:
		unique_ptr<Placement> placement;

		virtual void declare() override {
			auto& color = get<Color>(dr[0]);
			auto& tp = get<ThreadPool>(dr[1]);
			auto& tree = get<Tree>(dr[2]);
			auto& verbose = get<int>(dr[3]);
			placement.reset(new Placement(color, tp, tree, verbose));
		}

		virtual void run() override {
			*placement << "Score: " << placement->scoreTree() << std::endl;
			returned = dr[2]->returnValue<Tree>();
		}
	};

	class RandomTaxonOrderSubprocedure : public virtual Subprocedure {
	public:
		RandomTaxonOrderSubprocedure(S random, S nTaxa) {
			dr = { random, nTaxa };
		}

	protected:
		virtual void declare() override {}

		virtual void run() override {
			auto& random = get<Random>(dr[0]);
			auto& nTaxa = get<size_t>(dr[1]);
			auto taxa = random.randomTaxonOrder(nTaxa);
			prepareReturn(std::move(taxa));
		}
	};

	class RTO_RP_ST_Subprocedure : public virtual Subprocedure {
	public:
		RTO_RP_ST_Subprocedure(S random, S nTaxa, S color, S tp, S tree, S verbose) {
			dr = { random, nTaxa, color, tp, tree, verbose };
		}

	protected:
		shared_ptr<RandomTaxonOrderSubprocedure> rto;
		shared_ptr<RecursivePlacementSubprocedure> rp;
		unique_ptr<ScoreTreeSubprocedure> st;

		virtual void declare() override {
			auto& random = dr[0];
			auto& nTaxa = dr[1];
			auto& color = dr[2];
			auto& tp = dr[3];
			auto& tree = dr[4];
			auto& verbose = dr[5];
			rto.reset(new RandomTaxonOrderSubprocedure(random, nTaxa));
			rp.reset(new RecursivePlacementSubprocedure(color, tp, tree, rto, verbose));
			st.reset(new ScoreTreeSubprocedure(color, tp, rp, verbose));
		}

		virtual void run() override {
			returned = st->template returnValue<Tree>();
		}

	};
	
	class ConstrainedDPPlacementSubprocedure : public virtual Subprocedure {
	public:
		ConstrainedDPPlacementSubprocedure(S random, S nTaxa, S verbose, vector<S> trees) {
			dr = { random, nTaxa, verbose };
			rr = trees;
		}

	protected:
		unique_ptr<ConstrainedDP> dp;

		virtual void declare() override {
			auto& random = get<Random>(dr[0]);
			auto& nTaxa = get<size_t>(dr[1]);
			auto& verbose = get<int>(dr[2]);
			dp.reset(new ConstrainedDP(random, nTaxa, verbose));
		}

		virtual void run() override {
			for (auto& treeSubproc : rr) {
				auto& tree = get<Tree>(treeSubproc);
				dp->addTree(tree);
			}
			auto resultTree = dp->optimalUnrootedTree(ZERO);
			prepareReturn(std::move(resultTree));
		}
	};
	

	shared_ptr<Tree> oldRecursivePlacementProcedure(auto& data, size_t nTaxa, size_t nThreads, int verbose) {
		Random random;
		Log log(verbose);
		auto taxa = random.randomTaxonOrder(nTaxa);
		Color color(&data);
		ThreadPool tp(nThreads, 0, data.nElements);
		Tree tree;
		Recursive rp(verbose);
		rp(color, tp, tree, taxa);
		Placement placement(color, tp, tree, verbose);
		tree.displaySimpleNewick<double, double>(log << "Recursive placement tree: ", "", "") << std::endl;
		log << "Recursive placement score: " << placement.scoreTree() << std::endl;
		shared_ptr<Tree> pTree(new Tree(std::move(tree)));
		return pTree;
	}

	shared_ptr<Tree> recursivePlacementProcedure(auto& data, size_t _nTaxa, size_t nThreads, int _verbose) {
		S random = pInput<Random>();
		S nTaxa = pInput<size_t>(_nTaxa);
		S color = pInput<Color>(&data);
		S tp = pInput<ThreadPool>(nThreads, 0, data.nElements);
		S tree = pInput<Tree>();
		S verbose = pInput<int>(_verbose);
		S rto_rp_st( new RTO_RP_ST_Subprocedure(random, nTaxa, color, tp, tree, verbose) );
		return rto_rp_st->template returnValue<Tree>();
	}

	shared_ptr<Tree> constrainedDPProcedure(auto& data, size_t _nTaxa, size_t nTrees, size_t nThreads, int _verbose) {
		S random = pInput<Random>();
		S nTaxa = pInput<size_t>(_nTaxa);
		S color = pInput<Color>(&data);
		S tp = pInput<ThreadPool>(nThreads, 0, data.nElements);
		vector<S> trees;
		for (size_t iTree: std::views::iota((size_t) 0, nTrees)) trees.push_back(pInput<Tree>());
		S verbose = pInput<int>(_verbose);
		S taxa(new RandomTaxonOrderSubprocedure(random, nTaxa));
		vector<S> rto_rp_st;
		for (S& tree: trees) rto_rp_st.emplace_back( new RTO_RP_ST_Subprocedure(random, nTaxa, color, tp, tree, verbose) );
		S cdp(new ConstrainedDPPlacementSubprocedure(random, nTaxa, verbose, rto_rp_st));
		S sc(new ScoreTreeSubprocedure(color, tp, cdp, verbose));
		return sc->template returnValue<Tree>();
	}
};


};

#endif
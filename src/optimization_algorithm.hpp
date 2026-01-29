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
		LogInfo const v1log(verbose + 1);
		LogInfo const v2log(verbose + 2);
		log << "Recursive placement..." << std::endl;
		for (size_t iTaxon : taxa) {
			if (nLeaves < 3) tree.emplaceRoot(common::AnnotatedBinaryTree::LEAF_ID, iTaxon);
			else {
				v1log << "Placing " << nLeaves + 1 << "-th taxon (" << common::taxonName2ID[iTaxon] << ")" << std::endl;
				tree.displaySimpleNewick(v2log << "Current tree: ") << std::endl;
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
	
	template<class T> static T& get(S& s) try { return *(s->returnValue<T>()); }
	catch (...) { throw std::runtime_error("Error: Error occur in S::get."); }

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
public:
	using Log = common::LogInfo;
	using score_t = Attributes::score_t;
	using Tree = common::AnnotatedBinaryTree;
	using Random = Attributes::Random;
	using Color = Attributes::Color;
	using Data = Color::SharedConstData;
	template<typename T> using Scheduler = thread_pool::SimpleScheduler<T>;
	using ThreadPool = thread_pool::ThreadPool<score_t, Scheduler>;
	using Placement = Attributes::Placement;
	using Recursive = RecursivePlacement<Placement>;
	using ConstrainedDP = constrained_dp_algorithm::ConstrainedDP<score_t, Random>;
	static score_t constexpr ZERO = Placement::ZERO;
	static score_t constexpr EPSILON = Color::EPSILON;

	class RandomTaxonOrderSubprocedure : public virtual Subprocedure {
		size_t nTaxa;
	public:
		template<typename... Args>RandomTaxonOrderSubprocedure(S random, size_t nTaxa, Args... args) : nTaxa(nTaxa) {
			rr = { random, args... };
		}

	protected:
		virtual void declare() override {}

		virtual void run() override {
			auto& random = get<Random>(rr[0]);
			auto taxa = random.randomTaxonOrder(nTaxa);
			prepareReturn(std::move(taxa));
		}
	};

	class RecursivePlacementSubprocedure : public virtual Subprocedure {
		int verbose;

	public:
		template<typename... Args>RecursivePlacementSubprocedure(S color, S tp, S tree, S taxa, int verbose, Args... args) : verbose(verbose) {
			rr = { color, tp, tree, taxa, args... };
		}

	protected:
		unique_ptr<Recursive> rp;
		
		virtual void declare() override {;
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
		int verbose;

	public:
		template<typename... Args>ScoreTreeSubprocedure(S color, S tp, S tree, int verbose, Args... args): verbose(verbose) {
			dr = { color, tp, tree };
			rr = { args... };
		}

	protected:
		unique_ptr<Placement> placement;

		virtual void declare() override {
			auto& color = get<Color>(dr[0]);
			auto& tp = get<ThreadPool>(dr[1]);
			auto& tree = get<Tree>(dr[2]);
			placement.reset(new Placement(color, tp, tree, verbose));
		}

		virtual void run() override {
			*placement << "Score: " << placement->scoreTree() << std::endl;
			shared_ptr<Tree> tree = dr[2]->returnValue<Tree>();
			common::LogInfo const vlog(verbose + 1);
			tree->displaySimpleNewick(vlog) << std::endl;
			returned = tree;
		}
	};

	class RTO_RP_ST_Subprocedure : public virtual Subprocedure {
	public:
		template<typename... Args> RTO_RP_ST_Subprocedure(S random, S color, S tp, S tree, size_t nTaxa, int verbose, Args... args){
			S rto(new RandomTaxonOrderSubprocedure(random, nTaxa));
			S rp(new RecursivePlacementSubprocedure(color, tp, tree, rto, verbose));
			S st(new ScoreTreeSubprocedure(color, tp, rp, verbose));
			rr = { st, args... };
		}

	protected:
		virtual void declare() override {}

		virtual void run() override {
			returned = rr[0]->template returnValue<Tree>();
		}

	};

	class DeclareConstrainedDPSubprocedure : public virtual Subprocedure {
		size_t nTaxa;
		int verbose;
		
	public:
		template<typename... Args> DeclareConstrainedDPSubprocedure(S random, size_t nTaxa, int verbose, Args... args) : nTaxa(nTaxa), verbose(verbose) {
			rr = { random, args... };
		}

	protected:
		virtual void declare() override {}

		virtual void run() override {
			auto& random = get<Random>(rr[0]);
			returned = shared_ptr<ConstrainedDP>(new ConstrainedDP(random, nTaxa, verbose));
		}
	};
	
	class RunConstrainedDPSubprocedure : public virtual Subprocedure {
		vector<S> trees;

	public:
		template<typename... Args> RunConstrainedDPSubprocedure(S dcdp, vector<S> trees, Args... args) : trees(trees) {
			rr = { dcdp, args... };
			for (S const& tree: trees) rr.push_back(tree);
		}

	protected:
		virtual void declare() override {}

		virtual void run() override {
			auto& cdp = get<ConstrainedDP>(rr[0]);
			cdp << "Running constrained dynamic programming algorithm..." << std::endl;
			for (auto& treeSubproc : trees) {
				auto& tree = get<Tree>(treeSubproc);
				cdp.addTree(tree);
			}
			auto resultTree = cdp.optimalUnrootedTree(ZERO);
			prepareReturn(std::move(resultTree));
		}
	};

	class ConstrainedDPSubprocedure : public virtual Subprocedure {
		size_t nTaxa;
		int verbose;
		vector<S> trees;

	public:
		template<typename... Args> ConstrainedDPSubprocedure(S random, vector<S> trees, size_t nTaxa, int verbose, Args... args) : trees(trees), nTaxa(nTaxa), verbose(verbose) {
			dr = { random };
			rr = { args... };
			for (S const& tree : trees) rr.push_back(tree);
		}

	protected:
		unique_ptr<ConstrainedDP> dp;

		virtual void declare() override {
			auto& random = get<Random>(dr[0]);
			dp.reset(new ConstrainedDP(random, nTaxa, verbose));
		}

		virtual void run() override {
			for (auto& treeSubproc : trees) {
				auto& tree = get<Tree>(treeSubproc);
				dp->addTree(tree);
			}
			auto resultTree = dp->optimalUnrootedTree(ZERO);
			prepareReturn(std::move(resultTree));
		}
	};

	class RCDP_ST_Subprocedure : public virtual Subprocedure{
	public:
		template<typename... Args> RCDP_ST_Subprocedure(S dcdp, vector<S> trees, S color, S tp, int verbose, Args... args) {
			S rcdp(new RunConstrainedDPSubprocedure(dcdp, trees));
			S st(new ScoreTreeSubprocedure(color, tp, rcdp, verbose));
			rr = { st, args... };
		}

	protected:
		virtual void declare() override {}

		virtual void run() override {
			returned = rr[0]->template returnValue<Tree>();
		}
	};

	class RTO_RP_ST_RCDP_ST_Subprocedure : public virtual Subprocedure {
	public:
		template<typename... Args> RTO_RP_ST_RCDP_ST_Subprocedure(vector<S> trees, S random, S color, S tp, S dcdp, size_t nTaxa, int verbose, Args... args) {
			vector<S> rto_rp_st;
			for (S& tree: trees) rto_rp_st.emplace_back(new RTO_RP_ST_Subprocedure(random, color, tp, tree, nTaxa, verbose + 1));
			S rcdp_st(new RCDP_ST_Subprocedure(dcdp, rto_rp_st, color, tp, verbose));
			rr = { rcdp_st, args... };
		}

	protected:
		virtual void declare() override {}

		virtual void run() override {
			returned = rr[0]->template returnValue<Tree>();
		}

	};

	static shared_ptr<Tree> recursivePlacementProcedure(shared_ptr<Random> rnd, auto& data, size_t nTaxa, size_t nThreads, int verbose) {
		S random = shared_ptr<Input<Random> >(new Input<Random>(rnd));
		S color = pInput<Color>(&data);
		S tp = pInput<ThreadPool>(nThreads, 0, data.nElements);
		S tree = pInput<Tree>();
		S rto_rp_st( new RTO_RP_ST_Subprocedure(random, color, tp, tree, nTaxa, verbose) );
		return rto_rp_st->template returnValue<Tree>();
	}

	static shared_ptr<Tree> constrainedDPProcedure(S& dcdp, S& random, Data const& data, size_t nTaxa, size_t nTrees, size_t nThreads, int verbose) {
		S color = pInput<Color>(&data);
		S tp = pInput<ThreadPool>(nThreads, 0, data.nElements);
		vector<S> trees;
		for (size_t iTree : std::views::iota((size_t)0, nTrees)) trees.push_back(pInput<Tree>());
		S rto_rp_st_rcdp_sc(new RTO_RP_ST_RCDP_ST_Subprocedure(trees, random, color, tp, dcdp, nTaxa, verbose));
		return rto_rp_st_rcdp_sc->template returnValue<Tree>();
	}


	template<typename... Args> static shared_ptr<Tree> heuristSearch(Data const& data, size_t nTaxa, size_t r, size_t s, size_t nThreads, int verbose, shared_ptr<Tree> f(S&, S&, Data const&, size_t, size_t, size_t, int, Args...), Args... args) {
		S random = pInput<Random>();
		S dcdp(new DeclareConstrainedDPSubprocedure(random, nTaxa, verbose));
		shared_ptr<Tree> tree = f(dcdp, random, data, nTaxa, r, nThreads, verbose + 1, args...);
		score_t lastScore, newScore = tree->get<score_t>(Tree::SCORE);
		do {
			lastScore = newScore;
			tree = f(dcdp, random, data, nTaxa, s, nThreads, verbose + 1, args...);
			newScore = tree->get<score_t>(Tree::SCORE);
		} while (newScore > lastScore + EPSILON);
		return tree;
	}
};


};

#endif
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
using std::cerr;
using std::endl;
using common::LogInfo;

template<placement_algorithm::PLACEMENT_ALGORITHM PlacementAlgorithm> class RecursivePlacement : public LogInfo{
public:
	RecursivePlacement(int verbose = LogInfo::DEFAULT_VERBOSE) : LogInfo(verbose){}

	common::AnnotatedBinaryTree& operator()(typename PlacementAlgorithm::Color& color, typename PlacementAlgorithm::ThreadPool& threadPool, common::AnnotatedBinaryTree& tree, std::ranges::range auto const& taxa) {
		size_t nLeaves = tree.leaves().size();
		LogInfo const v1log(verbose + 1);
		LogInfo const v2log(verbose + 2);
		log() << "Recursive placement..." << endl;
		for (size_t iTaxon : taxa) {
			if (nLeaves < 3) tree.emplaceRoot(common::AnnotatedBinaryTree::LEAF_ID, iTaxon);
			else {
				v1log.log() << "Placing " << nLeaves + 1 << "-th taxon (" << common::taxonName2ID[iTaxon] << ")" << endl;
				tree.displaySimpleNewick(v2log.log() << "Current tree: ") << endl;
				PlacementAlgorithm placement(color, threadPool, tree, verbose + 1);
				placement.place(iTaxon);
			}
			nLeaves++;
		}
		return tree;
	}
};


using R = std::any;

class ProcedureBase {
public:
	class Subprocedure {
	protected:
		vector<shared_ptr<Subprocedure> > dep;
		bool finished = false;
	
	public:
		template<typename... Args> Subprocedure(Args... args): dep(args...) {}

		virtual ~Subprocedure() = default;

		virtual vector<shared_ptr<Subprocedure> > dependencies() { return dep; }
		
		virtual bool run() { return true; }

		virtual R result() { return R(); }

		bool recursivelyRun() {
			if (finished) return true;
			bool dependencyDone = true;
			for (auto& s : dependencies()) dependencyDone &= s->recursivelyRun();
			if (!dependencyDone) return false;
			return finished = run();
		}

		void addDependency() {}

		template<typename... Args> void addDependency(shared_ptr<Subprocedure> s, Args... args) { 
			dep.push_back(s);
			if constexpr (sizeof...(args) > 0) addDependency(args...);
		}

		template<typename T> static R prepResult(T&& result) {
			return shared_ptr<T>(new T(result));
		}

		template<typename T, typename... Args> static R prepResult(Args... args) {
			return shared_ptr<T>(new T(args...));
		}
		
		R getResult() {
			recursivelyRun();
			// if (!finished) finished = recursivelyRun(); // May have new thread jobs!
			// if (!finished) throw std::logic_error("Subprocedure: Dead lock!");
			R res = result();
			return res;
		}
	};
	using S = shared_ptr<Subprocedure>;

	template<typename T> class Input2 : public virtual Subprocedure {
		R res;

	public:
		virtual bool run() override { return true; }

		virtual R result() override { return res; }

		Input2(T&& t) : res(shared_ptr<T>(new T(t))) {}

		template<typename... Args> Input2(Args... args) : res(shared_ptr<T>(new T(std::forward<Args>(args)...))) {}
	};

	template<typename T, typename... Args> static shared_ptr<Input2<T> > sInput(Args... args) {
		Input2<T> *p = new Input2<T>(std::forward<Args>(args)...);
		shared_ptr<Input2<T> > s(p);
		return s;
	}

	class Block : public virtual Subprocedure {
	protected:
		S start;
		vector<S> units;

		Block(Subprocedure *start): start(start) { addDependency(this->start); }

		Block(S const& start) : start(start) { addDependency(this->start); }

		void initialize(){}

		template<typename... Args> void initialize(S unit, Args... args){
			initialize(std::forward<Args>(args)...);
			unit->addDependency(start);
			units.push_back(unit);
		}

		template<typename... Args> void initialize(vector<S> _units, Args... args) {
			initialize(std::forward<Args>(args)...);
			for (S& unit: _units){
				unit->addDependency(start);
				units.push_back(unit);
			}
		}

		template<typename... Args> Block(Subprocedure* start, Args... args) : Block(start) {
			initialize(std::forward<Args>(args)...);
		}

		template<typename... Args> Block(S const& start, Args... args) : Block(start) {
			initialize(std::forward<Args>(args)...);
		}
	public:
		Block() : start(new Subprocedure) {}

		template<typename... Args> Block(Args... args): Block() {
			initialize(std::forward<Args>(args)...);
		}

		virtual vector<shared_ptr<Subprocedure> > dependencies() override {
			vector<shared_ptr<Subprocedure> > res = dep;
			
			for (S& unit : units) res.push_back(unit);
			/*
			for (S& unit : units) {
				for (S& s : unit->dependencies()) {
					res.push_back(s);
				}
			}
			*/
			return res;
		}

		virtual R result() override { return units[0]->result(); }
	};

	class DeclareAndCall : public virtual Block {
	public:
		virtual S declaredType() { return start; }

		virtual R result() = 0;
	};

	class NewThread : public virtual Block {
	public:
		class Start : public virtual Subprocedure {
		public:
			bool lock = true;

			virtual bool run() { return lock; }
		};

	protected:
		Start* realStart = new Start(); // weak pointer
		unique_ptr<std::thread> thrd;
		LogInfo log;

	public:
		template<typename... Args> NewThread(int verbose, Args... args): log(verbose), Block(realStart, std::forward<Args>(args)...) {}

		virtual bool run() override {
			realStart->lock = false;
			thrd.reset(new std::thread([&] {
				log << "Spawn new thread (" << std::this_thread::get_id() << " for paralleled task...";
				recursivelyRun();
			}));
			return true;
		}

		virtual R result() override {
			thrd->join();
			return units[0]->result();
		}
	};

	template<class T, class Derived> requires std::derived_from<Derived, Subprocedure> static T& get(shared_ptr<Derived>& s) try {
		shared_ptr<T> p = std::any_cast<shared_ptr<T>>(s->getResult());
		if (!p) {
			std::cerr << "Fatal bug: Empty result! Forgetting dependencies?" << endl;
			exit(-1);
		}
		return *p;
	}
	catch (...) { throw std::runtime_error("Error: Error occur in S::get."); }
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
	using Log = LogInfo;
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

	class RandomTaxonOrderSubprocedure2 : public virtual Subprocedure {
		size_t nTaxa;
		S random;
		R res;

	public:
		template<typename... Args> RandomTaxonOrderSubprocedure2(S random, size_t nTaxa, Args... args) : nTaxa(nTaxa), random(random) {
			addDependency(random, args...);
		}

		virtual bool run() override {
			auto& random = get<Random>(this->random);
			res = prepResult(random.randomTaxonOrder(nTaxa));
			return true;
		}

		virtual R result() override { return res; }
	};

	class RecursivePlacementCall : public virtual DeclareAndCall {
		using Declared = Input2<Recursive>;

		S color, threadpool, tree, taxa;
		R res;

		shared_ptr<Declared> declared;

	public:
		virtual bool run() override {
			auto& color = get<Color>(this->color);
			auto& threadpool = get<ThreadPool>(this->threadpool);
			auto& tree = get<Tree>(this->tree);
			auto& taxa = get<vector<size_t> >(this->taxa);
			auto& rp = get<Recursive>(this->declared);
			rp(color, threadpool, tree, taxa);
			return true;
		}

		virtual R result() { return tree->getResult(); }

		template<typename... Args> RecursivePlacementCall(Subprocedure* start, S color, S threadpool, S tree, S taxa, Args... args) : Block(start), color(color), threadpool(threadpool), tree(tree), taxa(taxa) {
			declared = std::dynamic_pointer_cast<Declared>(this->start);
			addDependency(color, threadpool, tree, taxa, args...);
		}

		template<typename... Args> RecursivePlacementCall(int verbose, S color, S threadpool, S tree, S taxa, Args... args) : RecursivePlacementCall(new Declared(verbose, std::forward<Args>(args)...), color, threadpool, tree, taxa, args...) {}
	};

	class ScoreTreeCall : public virtual DeclareAndCall {
	public:
		class Declared : public virtual Subprocedure {
			int verbose;
			S color, threadpool, tree;
			R res;

		public:
			template<typename... Args> Declared(S color, S threadpool, S tree, int verbose, Args... args) : color(color), threadpool(threadpool), tree(tree), verbose(verbose) {
				addDependency(color, threadpool, tree, args...);
			}

			virtual bool run() override {
				auto& color = get<Color>(this->color);
				auto& threadpool = get<ThreadPool>(this->threadpool);
				auto& tree = get<Tree>(this->tree);
				//res = prepResult<Placement>(color, threadpool, tree, verbose);
				res = prepResult(Placement(color, threadpool, tree, verbose));
				return true;
			}

			virtual R result() override { return res; }
			friend ScoreTreeCall;
		};

	private:
		shared_ptr<Declared> declared;

	public:

		virtual bool run() override {
			Placement& placement = get<Placement>(declared);
			placement.log() << "Score: " << placement.scoreTree() << endl;
			return true;
		}

		virtual R result() { 
			std::any tree = declared->tree->getResult();
			Log const vlog(declared->verbose + 1);
			get<Tree>(declared->tree).displaySimpleNewick(vlog.log()) << endl;
			return tree;
		}

		template<typename... Args> ScoreTreeCall(Subprocedure* start, Args... args) : Block(start) {
			declared = std::dynamic_pointer_cast<Declared>(this->start);
			addDependency(args...);
		}

		template<typename... Args> ScoreTreeCall(S color, S threadpool, S tree, int verbose, Args... args) : ScoreTreeCall(new Declared(color, threadpool, tree, verbose, std::forward<Args>(args)...)) {}
	};

	class RTO_RP_ST_Block : public virtual Block {
	public:
		template<typename... Args> RTO_RP_ST_Block(S random, S color, S threadpool, S tree, size_t nTaxa, int verbose, Args... args) {
			addDependency(std::forward<Args>(args)...);
			S taxa(new RandomTaxonOrderSubprocedure2(random, nTaxa));
			S rpTree(new RecursivePlacementCall(verbose, color, threadpool, tree, taxa));
			S stTree(new ScoreTreeCall(color, threadpool, rpTree, verbose));
			initialize(taxa, rpTree, stTree);
		}

		template<typename... Args> RTO_RP_ST_Block(S random, S color, S threadpool, size_t nTaxa, int verbose, Args... args) : RTO_RP_ST_Block(random, color, threadpool, sInput<Tree>(), nTaxa, verbose, std::forward<Args>(args)...) {}

		template<typename... Args> RTO_RP_ST_Block(S random, S color, S tree, size_t nThread, size_t iElementBegin, size_t iElementEnd, size_t nTaxa, int verbose, Args... args) : RTO_RP_ST_Block(random, color, sInput<ThreadPool>(nThread, iElementBegin, iElementEnd), tree, nTaxa, verbose, std::forward<Args>(args)...) {}

		template<typename... Args> RTO_RP_ST_Block(S random, S color, size_t nThread, size_t iElementBegin, size_t iElementEnd, size_t nTaxa, int verbose, Args... args) : RTO_RP_ST_Block(random, color, sInput<ThreadPool>(nThread, iElementBegin, iElementEnd), nTaxa, verbose, std::forward<Args>(args)...) {}
	};

	class NT_RTO_RP_ST_Block : public virtual NewThread {
		template<typename... Args> NT_RTO_RP_ST_Block(int verbose, Args... args) : NewThread(verbose, S(new RTO_RP_ST_Block(std::forward<Args>(args)...))) {}

		template<typename... Args> NT_RTO_RP_ST_Block(S random, S color, size_t nThread, size_t iElementBegin, size_t iElementEnd, size_t nTaxa, int verbose, Args... args) : NewThread(verbose, S(new RTO_RP_ST_Block(random, color, nThread, iElementBegin, iElementEnd, nTaxa, verbose, std::forward<Args>(args)...))) {}
	};

	class ConstrainedDPCall : public virtual DeclareAndCall {
	public:	
		class Declared : public virtual Subprocedure {
			size_t nTaxa;
			int verbose;
			S random;
			R res;

		public:
			template<typename... Args> Declared(S random, size_t nTaxa, int verbose, Args... args) : random(random), nTaxa(nTaxa), verbose(verbose) {
				addDependency(random, args...);
			}

			virtual bool run() override {
				auto& random = get<Random>(this->random);
				res = prepResult<ConstrainedDP>(random, nTaxa, verbose);
				return true;
			}

			virtual R result() override { return res; }
		};

	private:
		shared_ptr<Declared> declared;
		vector<S> trees;
		R res;

	public:

		virtual bool run() override {
			auto& cdp = get<ConstrainedDP>(declared);
			for (S& tree: trees) {
				cdp.addTree(get<Tree>(tree));
			}
			res = prepResult(cdp.optimalUnrootedTree(ZERO));
			return true;
		}

		virtual R result() override { return res; }

		template<typename... Args> ConstrainedDPCall(S start, const vector<S> &trees, Args... args) : Block(start), trees(trees) {
			declared = std::dynamic_pointer_cast<Declared>(this->start);
			for (S const& tree: trees) addDependency(tree);
			addDependency(args...);
		}

		template<typename... Args> ConstrainedDPCall(Subprocedure* start, const vector<S>& trees, Args... args) : ConstrainedDPCall(S(start), trees, args...) {}

		template<typename... Args> ConstrainedDPCall(S random, size_t nTaxa, int verbose, const vector<S>& trees, Args... args) : ScoreTreeCall(new Declared(random, nTaxa, verbose, std::forward<Args>(args)...)) {}
	};
	
	class CCDP_ST_Block : public virtual Block {
	public:
		template<typename... Args> CCDP_ST_Block(S dCDP, vector<S> trees, S color, S threadpool, int verbose, Args... args) : Block(new Subprocedure(args...)) {
			S cdpTree(new ConstrainedDPCall(dCDP, trees));
			S stTree(new ScoreTreeCall(color, threadpool, cdpTree, verbose));
			initialize(cdpTree, stTree);
		}
	};

	class RTO_RP_ST_CCDP_ST_Block : public virtual Block {
	public:
		template<typename... Args> RTO_RP_ST_CCDP_ST_Block(S dCDP, S random, S color, S threadpool, size_t nTaxa, size_t nTrees, int verbose, Args... args) : Block(new Subprocedure(args...)) {
			vector<S> trees;
			for (size_t iTree : std::views::iota((size_t)0, nTrees)) trees.emplace_back(new RTO_RP_ST_Block(random, color, threadpool, nTaxa, verbose + 1));
			S tree(new CCDP_ST_Block(dCDP, trees, color, threadpool, verbose));
			initialize(trees, tree);
		}
	};
	
	static Tree recursivePlacementProcedure(S random, Data const& data, size_t nTaxa, size_t nTrees, size_t nThreads, int verbose, size_t iElementBegin, size_t iElementEnd) {
		S color = sInput<Color>(&data);
		S threadpool = sInput<ThreadPool>(nThreads, iElementBegin, iElementEnd);
		S tree = sInput<Tree>();
		S workflow(new RTO_RP_ST_Block(random, color, threadpool, tree, nTaxa, verbose));
		return get<Tree>(workflow);
	}
	
	static Tree recursivePlacementProcedure(S random, Data const& data, size_t nTaxa, size_t nTrees, size_t nThreads, int verbose) {
		return recursivePlacementProcedure(random, data, nTaxa, nTrees, nThreads, verbose, 0, data.nElements);
	}
	
	static Tree constrainedDPProcedure(S dCDP, S random, Data const& data, size_t nTaxa, size_t nTrees, size_t nThreads, int verbose, size_t iElementBegin, size_t iElementEnd) {
		Log log(verbose), vlog(verbose + 1);
		log.log() << "Running costrained dynamic programming (CDP) procedure..." << endl;
		vlog.log() << "#Rounds: " << nTrees << endl;
		vlog.log() << "#Threads: " << nThreads << endl;
		S color = sInput<Color>(&data);
		S threadpool = sInput<ThreadPool>(nThreads, iElementBegin, iElementEnd);
		S tree(new RTO_RP_ST_CCDP_ST_Block(dCDP, random, color, threadpool, nTaxa, nTrees, verbose));
		return get<Tree>(tree);
	}

	static Tree constrainedDPProcedure(S dCDP, S random, Data const& data, size_t nTaxa, size_t nTrees, size_t nThreads, int verbose) {
		return constrainedDPProcedure(dCDP, random, data, nTaxa, nTrees, nThreads, verbose, 0, data.nElements);
	}

	static Tree defaultProcedure(S dCDP, S random, Data const& data, size_t nTaxa, size_t nTrees, size_t nThreads, int verbose) {
		return constrainedDPProcedure(dCDP, random, data, nTaxa, nTrees, nThreads, verbose);
	}

	template<typename... Args> static Tree heuristSearch(Data const& data, size_t nTaxa, size_t r, size_t s, size_t nThreads, int verbose, Tree f(S, S, Data const&, size_t, size_t, size_t, int, Args...), Args... args) {
		S random = sInput<Random>();
		S dCDP(new ConstrainedDPCall::Declared(random, nTaxa, verbose));
		Tree tree = f(dCDP, random, data, nTaxa, r, nThreads, verbose + 1, args...);
		score_t lastScore, newScore = tree.get<score_t>(Tree::SCORE);
		do {
			lastScore = newScore;
			tree = f(dCDP, random, data, nTaxa, s, nThreads, verbose + 1, args...);
			newScore = tree.get<score_t>(Tree::SCORE);
		} while (newScore > lastScore + EPSILON);
		return tree;
	}
};


};

#endif
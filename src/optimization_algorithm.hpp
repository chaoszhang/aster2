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
				tree.displaySimpleNewick(v2log.log() << "Current tree: ") << endl;
				v1log.log() << "Placing " << nLeaves + 1 << "-th taxon (" << common::taxonName2ID[iTaxon] << ")" << endl;
				PlacementAlgorithm placement(color, threadPool, tree, verbose + 1);
				placement.place(iTaxon);
			}
			nLeaves++;
		}
		tree.displaySimpleNewick(v1log.log() << "Recursively placed tree: ") << endl;
		return tree;
	}
};

template<COLORABLE C> struct DefaultProcedureAttributes {
	using Color = C;
	using score_t = Color::score_t;
	template<typename T> using Scheduler = thread_pool::SimpleScheduler<T>;
	using Placement = placement_algorithm::StepwiseColorPlacement<placement_algorithm::StepwiseColorPlacementDefaultAttributes<Color> >;
	using Random = common::Random<std::mt19937_64>;
};

template<typename Attributes> class Procedure {
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

	struct RP {
		struct Prereq {
			int verbose;
			vector<size_t> taxa;

			Prereq(int verbose, const vector<size_t>& taxa) : verbose(verbose), taxa(taxa) {}
		} p;

		RP(Prereq const& req) : p(req) {}

		Tree operator()(Color& color, ThreadPool& threadpool, Tree tree) {
			Log log(p.verbose);
			log.log() << "Running recursive placement (RP) procedure..." << endl;
			Recursive rp(p.verbose);
			rp(color, threadpool, tree, p.taxa);
			return tree;
		}

		Tree operator()(Color& color, ThreadPool& threadpool) {
			Tree tree;
			return (*this)(color, threadpool, tree, p.taxa);
		}
	};

	struct TP {
		struct Prereq {
			int verbose;
			size_t nThreads, iElementBegin, iElementEnd;

			Prereq(int verbose, size_t nThreads, size_t iElementBegin, size_t iElementEnd) : verbose(verbose), nThreads(nThreads), iElementBegin(iElementBegin), iElementEnd(iElementEnd) {}
		} p;

		TP(Prereq const& req) : p(req) {}

		shared_ptr<ThreadPool> operator()() {
			Log log(p.verbose), vlog(p.verbose + 1);
			log.log() << "Generating thread pool (TP) with " << p.nThreads << " thread(s)..." << endl;
			vlog.log() << "Covering data range [" << p.iElementBegin << "," << p.iElementEnd << ")" << endl;
			shared_ptr<ThreadPool> threadpool(new ThreadPool(p.nThreads, p.iElementBegin, p.iElementEnd) );
			return threadpool;
		}
	};

	template<typename X> struct TP_X {
		struct Prereq {
			typename X::Prereq x;
			typename TP::Prereq tp;
			Prereq(X::Prereq const& x, typename TP::Prereq const& tp) : x(x), tp(tp) {}
		} p;

		TP_X(Prereq const& req) : p(req) {}

		Tree operator()(Color& color, Tree& tree) {
			TP tp(p.tp);
			X x(p.x);
			shared_ptr<ThreadPool> threadpool = tp();
			return x(color, *threadpool, tree);
		}

		Tree operator()(Color& color) {
			Tree tree;
			return (*this)(color, tree);
		}
	};

	template<typename X> struct Parallel_X {
		struct Prereq {
			int verbose;
			vector<typename X::Prereq> jobs;

			Prereq(int verbose, const vector<typename X::Prereq>& jobs) : verbose(verbose), jobs(jobs) {}
		} p;

		Parallel_X(Prereq const& req) : p(req) {}

		vector<Tree> operator()(Color& color) {
			size_t nJobs = p.jobs.size();
			vector<std::thread> thrds;
			vector<Tree> trees(nJobs);
			for (size_t iJob : std::views::iota((size_t)1, nJobs)) {
				Tree& tree = trees[iJob];
				thrds.emplace_back([this, iJob, &color, &tree]() {
					X job(p.jobs[iJob]);
					tree = job(color);
				});
			}
			X job0(p.jobs[0]);
			trees[0] = job0(color);
			for (auto& t : thrds) t.join();
			return trees;
		}
	};

	using TP_RP = TP_X<RP>;
	using Parallel_TP_RP = Parallel_X<TP_RP>;

	static typename Parallel_TP_RP::Prereq prereq_Parallel_TP_RP(Random &random, size_t nThreads, size_t iTotalElementBegin, size_t iTotalElementEnd, size_t nParallel, size_t nTaxa, int verbose) {
		vector<typename TP_RP::Prereq> jobs;
		size_t nElements = iTotalElementEnd - iTotalElementBegin;
		for (size_t iParallel : std::views::iota((size_t)0, nParallel)) {
			size_t iThreadBegin = iParallel * nThreads / nParallel;
			size_t iThreadEnd = (iParallel + 1) * nThreads / nParallel;
			size_t iElementBegin = iParallel * nElements / nParallel + iTotalElementBegin;
			size_t iElementEnd = (iParallel + 1) * nElements / nParallel + iTotalElementBegin;

			typename RP::Prereq rp(verbose + 1, random.randomTaxonOrder(nTaxa));
			typename TP::Prereq tp(verbose + 2, iThreadEnd - iThreadBegin, iElementBegin, iElementEnd);
			jobs.emplace_back(rp, tp);
		}
		typename Parallel_TP_RP::Prereq p(verbose, jobs);
		return p;
	}

	struct ST {
		struct Prereq {
			int verbose;

			Prereq(int verbose) : verbose(verbose) {}
		} p;

		ST(Prereq const& req) : p(req) {}

		Tree operator()(Color& color, ThreadPool& threadpool, Tree& tree) {
			Placement placement(color, threadpool, tree, p.verbose);
			placement.log() << "Score: " << placement.scoreTree() << endl;
			Log vlog(p.verbose + 1);
			tree.displaySimpleNewick(vlog.log()) << endl;
			return tree;
		}
	};

	struct TP_Sequencial_ST_CDP_ST {
		struct Prereq {
			int verbose;
			typename TP::Prereq tp;

			Prereq(int verbose, typename TP::Prereq const& tp) : verbose(verbose), tp(tp) {}
		} p;

		TP_Sequencial_ST_CDP_ST(Prereq const& req) : p(req) {}

		Tree operator()(Color& color, vector<Tree>& trees, ConstrainedDP& cdp) {
			TP tp(p.tp);
			shared_ptr<ThreadPool> threadpool = tp();
			for (Tree& tree : trees) {
				typename ST::Prereq stp(p.verbose + 1);
				ST st(stp);
				cdp.addTree(st(color, *threadpool, tree));
			}
			typename ST::Prereq stp(p.verbose);
			ST st(stp);
			Tree tree = cdp.optimalUnrootedTree(ZERO);
			return st(color, *threadpool, tree);
		}
	};

	static typename TP_Sequencial_ST_CDP_ST::Prereq prereq_TP_Sequencial_ST_CDP_ST(size_t nThreads, size_t iElementBegin, size_t iElementEnd, int verbose) {
		typename TP::Prereq tp(verbose + 1, nThreads, iElementBegin, iElementEnd);
		typename TP_Sequencial_ST_CDP_ST::Prereq p(verbose, tp);
		return p;
	}

	static Tree subsampleProcedure2(ConstrainedDP &cdp, Random &random, Color& color, Data const& data, size_t nTaxa, size_t nTrees, size_t nThreads, int verbose) {
		Log::setShowThread();
		Log log(verbose), vlog(verbose + 1);
		log.log() << "Running subsample procedure..." << endl;
		vlog.log() << "#Rounds: " << nTrees << endl;
		vlog.log() << "#Threads: " << nThreads << endl;
		size_t minElements = ARG.get<size_t>("subsample-min");
		size_t nElements = data.nElements;

		size_t _nPartitions = 1;
		while (_nPartitions <= nTrees && nElements / (_nPartitions + 1) >= minElements) _nPartitions++;
		size_t nRuns = (nTrees + _nPartitions - 1) / _nPartitions;
		// cerr << "nRuns = " << nRuns << endl;
		vector<Tree> trees;
		for (size_t iRun : std::views::iota((size_t)0, nRuns)) {
			size_t iTreeBegin = iRun * nTrees / nRuns;
			size_t iTreeEnd = (iRun + 1) * nTrees / nRuns;
			size_t nBatches = iTreeEnd - iTreeBegin;
			size_t nPartitions = (nBatches + nThreads - 1) / nThreads;
			// cerr << "nPartitions = " << nPartitions << endl;
			int newVerbose = (nBatches > 1) ? verbose + 1 : verbose;
			for (size_t iPartition : std::views::iota((size_t)0, nPartitions)) {
				size_t iElementBegin = iPartition * nElements / nPartitions;
				size_t iElementEnd = (iPartition + 1) * nElements / nPartitions;
				size_t iTreeBegin2 = iPartition * nBatches / nPartitions + iTreeBegin;
				size_t iTreeEnd2 = (iPartition + 1) * nBatches / nPartitions + iTreeBegin;
				//cerr << "iTreeBegin = " << iTreeBegin2 << endl;
				//cerr << "iTreeEnd = " << iTreeEnd2 << endl;
				//cerr << "iElementBegin = " << iElementBegin << endl;
				//cerr << "iElementEnd = " << iElementEnd << endl;
				vlog.log() << "Running " << iTreeEnd2 - iTreeBegin2 << " subsampled tree(s) partitioning element(s) [" << iElementBegin << "," << iElementEnd << ")..." << endl;
				typename Parallel_TP_RP::Prereq p = prereq_Parallel_TP_RP(random, nThreads, iElementBegin, iElementEnd, iTreeEnd2 - iTreeBegin2, nTaxa, verbose);
				Parallel_TP_RP job(p);

				for (Tree& tree : job(color)) trees.push_back(tree);
			}
		}
		typename TP_Sequencial_ST_CDP_ST::Prereq p = prereq_TP_Sequencial_ST_CDP_ST(nThreads, 0, data.nElements, verbose);
		TP_Sequencial_ST_CDP_ST job(p);
		return job(color, trees, cdp);
	}

	static Tree defaultProcedure2(ConstrainedDP& cdp, Random& random, Color& color, Data const& data, size_t nTaxa, size_t nTrees, size_t nThreads, int verbose) {
		return subsampleProcedure2(cdp, random, color, data, nTaxa, nTrees, nThreads, verbose);
	}

	static Tree heuristSearch2(Data const& data, size_t nTaxa, size_t r, size_t s, size_t nThreads, int verbose, Tree f(ConstrainedDP&, Random&, Color&, Data const&, size_t, size_t, size_t, int)) {
		Random random;
		ConstrainedDP cdp(random, nTaxa, verbose);
		Color color(&data);
		Tree tree = f(cdp, random, color, data, nTaxa, r, nThreads, verbose + 1);
		score_t lastScore, newScore = tree.get<score_t>(Tree::SCORE);
		do {
			lastScore = newScore;
			tree = f(cdp, random, color, data, nTaxa, s, nThreads, verbose + 1);
			newScore = tree.get<score_t>(Tree::SCORE);
		} while (newScore > lastScore + EPSILON);
		return tree;
	}
};


};

#endif
#ifndef OPTIMIZATION_ALGORITHM_HPP
#define OPTIMIZATION_ALGORITHM_HPP

#include "placement_algorithm.hpp"
#include "nni_algorithm.hpp"
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

ChangeLog logTwoStepPlacement("TwoStepPlacement",
	"2026-02-05", "Chao Zhang", "Initial version", "minor");

template<placement_algorithm::PLACEMENT_ALGORITHM PlacementAlgorithm, nni_algorithm::NNI_ALGORITHM NNIAlgorithm, bool UNROOTED = true> class TwoStepPlacement : public LogInfo {
	using Tree = common::AnnotatedBinaryTree;

	size_t nThreads, iElementBegin, iElementEnd;

public:
	TwoStepPlacement(size_t nThreads, size_t iElementBegin, size_t iElementEnd, int verbose = LogInfo::DEFAULT_VERBOSE) : LogInfo(verbose), nThreads(nThreads), iElementBegin(iElementBegin), iElementEnd(iElementEnd){}

	Tree operator()(typename PlacementAlgorithm::Color& color, std::ranges::range auto const& taxa) {
		LogInfo const v1log(verbose + 1);
		LogInfo const v2log(verbose + 2);
		LogInfo const v3log(verbose + 3);
		LogInfo const v4log(verbose + 4);
		Tree backbone;
		vector<size_t> taxonOrder;
		for (size_t iTaxon : taxa) taxonOrder.push_back(iTaxon);
		size_t nLeaves = taxonOrder.size();
		size_t nBackboneLeaves = std::sqrt(nLeaves * std::log2(nLeaves));
		log() << "Two-step placement..." << endl;
		{
			v1log.log() << "Building backbone..." << endl;
			RecursivePlacement<PlacementAlgorithm> backboneRP(verbose + 2);
			typename PlacementAlgorithm::ThreadPool threadpool(nThreads, iElementBegin, iElementEnd);
			backboneRP(color, threadpool, backbone, taxonOrder | std::views::take(nBackboneLeaves));
		}
		{
			v2log.log() << "Applying NNI moves to backbone..." << endl;
			typename NNIAlgorithm::ThreadPool threadpool3(nThreads, iElementBegin, iElementEnd);
			NNIAlgorithm nniAlg(color, threadpool3, backbone, verbose + 1);
			nniAlg.performNNI();
		}
		{
			v1log.log() << "Determining the branches of placement on backbone..." << endl;
			typename PlacementAlgorithm::ThreadPool threadpool(nThreads, iElementBegin, iElementEnd);
			std::unordered_map<Tree::Node*, vector<size_t> > subtreeTaxa;
			std::unordered_map<size_t, Tree::Node*> branchOfPlacement;
			{
				for (size_t iTaxon : taxonOrder | std::views::drop(nBackboneLeaves)) {
					PlacementAlgorithm placement(color, threadpool, backbone, verbose + 1);
					v2log.log() << "Determining the branch of placement for " << common::taxonName2ID[iTaxon] << "..." << endl;
					Tree::Node* branch = placement.optimalPlacement(iTaxon);
					if (UNROOTED && (branch == backbone.root()->leftChild() || branch == backbone.root()->rightChild())) branch = backbone.root();
					branchOfPlacement[iTaxon] = branch;
					subtreeTaxa[branch].push_back(iTaxon);
				}
			}
			v1log.log() << "Building per branch subtrees..." << endl;
			std::unordered_map<Tree::Node*, std::tuple<Tree, Tree::Node*, Tree::Node*> > subtrees;
			vector<size_t> orphans;
			for (size_t iTaxon : taxonOrder | std::views::drop(nBackboneLeaves)) {
				Tree::Node* backboneBranch = branchOfPlacement[iTaxon];
				if (subtrees.contains(backboneBranch)) continue;
				
				v2log.log() << "Building the per branch subtree where " << common::taxonName2ID[iTaxon] << " was placed..." << endl;
				Tree subtree = backbone.deepCopy();
				vector<size_t> const& taxa2place = subtreeTaxa[backboneBranch];
				Tree::Node *branch, *sister;
				std::unordered_set<Tree::Node*> badBranches;
				if (backboneBranch->isRoot()) {
					for (Tree::Node* node : subtree.nodes()) {
						if (node == subtree.root()) continue;
						if (UNROOTED && node == subtree.root()->leftChild()) continue;
						if (UNROOTED && node == subtree.root()->rightChild()) continue;
						badBranches.insert(node);
					}
					if (UNROOTED) {
						branch = subtree.root()->leftChild();
						sister = subtree.root()->rightChild();
					}
					else {
						branch = subtree.root();
						sister = nullptr;
					}
				}
				else {
					vector<bool> location;
					for (Tree::Node* node = backboneBranch; !node->isRoot(); node = node->parent()) {
						location.push_back(node == node->parent()->leftChild());
					}
					branch = subtree.root();
					for (bool isLeft : location | std::views::reverse) {
						branch = (isLeft) ? branch->leftChild() : branch->rightChild();
					}
					for (Tree::Node* node : subtree.nodes()) {
						if (node == branch) continue;
						badBranches.insert(node);
					}
					sister = (location[0]) ? branch->parent()->rightChild() : branch->parent()->leftChild();
				}
				for (size_t jTaxon : taxa2place) {
					v3log.log() << "Try placing " << common::taxonName2ID[jTaxon] << "..." << endl;
					PlacementAlgorithm placement(color, threadpool, subtree, verbose + 2);
					Tree::Node* placedBranch = placement.optimalPlacement(jTaxon);
					if (badBranches.contains(placedBranch)) {
						v3log.log() << "Failed! Postponed for later..." << endl;
						orphans.push_back(jTaxon);
					}
					else placedBranch->emplaceAbove(Tree::LEAF_ID, jTaxon);
					subtree.displaySimpleNewick(v4log.log()) << endl;
				}
				subtrees[backboneBranch] = {subtree, branch, sister};
			}
			v1log.log() << "Merging per branch subtrees..." << endl;
			for (auto const& element : subtrees) {
				Tree::Node* backboneBranch = element.first;
				auto& [subtree, branch, sister] = element.second;
				if (backboneBranch->isRoot()) {
					if (UNROOTED) {
						// branch == oldSubtree.root()->leftChild();
						// sister == oldSubtree.root()->rightChild();
						branch->swapLeftChildren(backbone.root()->leftChild());
						branch->swapRightChildren(backbone.root()->leftChild());
						sister->swapLeftChildren(backbone.root()->rightChild());
						sister->swapRightChildren(backbone.root()->rightChild());
						subtree.root()->swap(backbone.root());
					}
					else {
						// branch == oldSubtree.root();
						branch->swapLeftChildren(backbone.root());
						branch->swapRightChildren(backbone.root());
						branch->parent()->swap(backbone.root());
						subtree.root()->swap(backbone.root());
					}
				}
				else {
					if (backboneBranch == backboneBranch->parent()->rightChild()) backboneBranch->parent()->swapChildren();
					if (sister == sister->parent()->leftChild()) sister->parent()->swapChildren();
					branch->swapLeftChildren(backboneBranch);
					branch->swapRightChildren(backboneBranch);
					sister->parent()->leftChild()->swap(backboneBranch);
				}
			}

			v1log.log() << "Placing postponed taxa if any..." << endl;
			{
				RecursivePlacement<PlacementAlgorithm> orphanRP(verbose + 1);
				orphanRP(color, threadpool, backbone, orphans);
			}
		}
		backbone.displaySimpleNewick(v1log.log() << "Two-step placed tree: ") << endl;
		return backbone;
	}
};

template<COLORABLE C> struct DefaultProcedureAttributes {
	using Color = C;
	using score_t = Color::score_t;
	template<typename T> using Scheduler = thread_pool::SimpleScheduler<T>;
	using Placement = placement_algorithm::StepwiseColorPlacement<placement_algorithm::StepwiseColorPlacementDefaultAttributes<Color> >;
	using NNIAlg = nni_algorithm::StepwiseColorNNI<nni_algorithm::StepwiseColorNNIDefaultAttributes<Color> >;
	using Random = common::Random<std::mt19937_64>;
};

ChangeLog logProcedure("optimization_algorithm::Procedure",
	"2026-02-13", "Chao Zhang", "Switch to sequential subsample", "patch");

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
	using ThreadPool3 = thread_pool::ThreadPool<std::array<score_t, 3>, Scheduler>;
	using Recursive = RecursivePlacement<Placement>;
	using NNIAlg = Attributes::NNIAlg;
	using TwoStep = TwoStepPlacement<Placement, NNIAlg>;
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

	struct TSP {
		struct Prereq {
			int verbose;
			size_t nThreads, iElementBegin, iElementEnd;
			vector<size_t> taxa;

			Prereq(int verbose, size_t nThreads, size_t iElementBegin, size_t iElementEnd, const vector<size_t>& taxa) : verbose(verbose), nThreads(nThreads), iElementBegin(iElementBegin), iElementEnd(iElementEnd), taxa(taxa) {}
		} p;

		TSP(Prereq const& req) : p(req) {}

		Tree operator()(Color& color) {
			TwoStep tsp(p.nThreads, p.iElementBegin, p.iElementEnd, p.verbose);
			return tsp(color, p.taxa);
		}
	};

	template<typename X> struct TPX {
		struct Prereq {
			int verbose;
			size_t nThreads, iElementBegin, iElementEnd;

			Prereq(int verbose, size_t nThreads, size_t iElementBegin, size_t iElementEnd) : verbose(verbose), nThreads(nThreads), iElementBegin(iElementBegin), iElementEnd(iElementEnd) {}
		} p;

		TPX(Prereq const& req) : p(req) {}

		shared_ptr<X> operator()() {
			Log log(p.verbose), vlog(p.verbose + 1);
			log.log() << "Generating thread pool (TP) with " << p.nThreads << " thread(s)..." << endl;
			vlog.log() << "Covering data range [" << p.iElementBegin << "," << p.iElementEnd << ")" << endl;
			shared_ptr<X> threadpool(new X(p.nThreads, p.iElementBegin, p.iElementEnd) );
			return threadpool;
		}
	};

	using TP = TPX<ThreadPool>;
	using TP3 = TPX<ThreadPool3>;

	template<typename TPx, typename X> struct TPx_X {
		struct Prereq {
			typename X::Prereq x;
			typename TPx::Prereq tp;
			Prereq(X::Prereq const& x, typename TPx::Prereq const& tp) : x(x), tp(tp) {}
		} p;

		TPx_X(Prereq const& req) : p(req) {}

		Tree operator()(Color& color, Tree& tree) {
			TPx tp(p.tp);
			X x(p.x);
			auto threadpool = tp();
			return x(color, *threadpool, tree);
		}

		Tree operator()(Color& color) {
			Tree tree;
			return (*this)(color, tree);
		}
	};

	using TP_RP = TPx_X<TP, RP>;

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

	using Parallel_TP_RP = Parallel_X<TP_RP>;
	using Parallel_TSP = Parallel_X<TSP>;

	static typename TP_RP::Prereq prereq_TP_RP(Random& random, size_t nThreads, size_t iElementBegin, size_t iElementEnd, size_t nTaxa, int verbose) {
		typename RP::Prereq rp(verbose, random.randomTaxonOrder(nTaxa));
		typename TP::Prereq tp(verbose + 1, nThreads, iElementBegin, iElementEnd);
		typename TP_RP::Prereq p(rp, tp);
		return p;
	}

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

	static typename TSP::Prereq prereq_TSP(Random& random, size_t nThreads, size_t iElementBegin, size_t iElementEnd, size_t nTaxa, int verbose) {
		typename TSP::Prereq p(verbose, nThreads, iElementBegin, iElementEnd, random.randomTaxonOrder(nTaxa));
		return p;
	}

	static typename Parallel_TSP::Prereq prereq_Parallel_TSP(Random& random, size_t nThreads, size_t iTotalElementBegin, size_t iTotalElementEnd, size_t nParallel, size_t nTaxa, int verbose) {
		vector<typename TSP::Prereq> jobs;
		size_t nElements = iTotalElementEnd - iTotalElementBegin;
		for (size_t iParallel : std::views::iota((size_t)0, nParallel)) {
			size_t iThreadBegin = iParallel * nThreads / nParallel;
			size_t iThreadEnd = (iParallel + 1) * nThreads / nParallel;
			size_t iElementBegin = iParallel * nElements / nParallel + iTotalElementBegin;
			size_t iElementEnd = (iParallel + 1) * nElements / nParallel + iTotalElementBegin;

			typename TSP::Prereq tsp(verbose + 1, iThreadEnd - iThreadBegin, iElementBegin, iElementEnd, random.randomTaxonOrder(nTaxa));
			jobs.emplace_back(tsp);
		}
		typename Parallel_TSP::Prereq p(verbose, jobs);
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

	template<typename TPx, typename X> struct TPx_Sequential_X {
		struct Prereq {
			int verbose;
			typename TPx::Prereq tp;

			Prereq(int verbose, typename TPx::Prereq const& tp) : verbose(verbose), tp(tp) {}
		} p;

		TPx_Sequential_X(Prereq const& req) : p(req) {}

		vector<Tree> operator()(Color& color, vector<Tree>& input) {
			TPx tp(p.tp);
			auto threadpool = tp();
			vector<Tree> output;
			for (Tree& tree : input) {
				typename X::Prereq xp(p.verbose + 1);
				X x(xp);
				output.push_back(x(color, *threadpool, tree));
			}
			return output;
		}
	};

	using TP_ST = TPx_X<TP, ST>;

	using TP_Sequential_ST = TPx_Sequential_X<TP, ST>;

	static typename TP_ST::Prereq prereq_TP_ST(size_t nThreads, size_t iElementBegin, size_t iElementEnd, int verbose) {
		typename TP::Prereq tp(verbose + 1, nThreads, iElementBegin, iElementEnd);
		typename TP_ST::Prereq p(verbose, tp);
		return p;
	}

	static typename TP_Sequential_ST::Prereq prereq_TP_Sequential_ST(size_t nThreads, size_t iElementBegin, size_t iElementEnd, int verbose) {
		typename TP::Prereq tp(verbose + 1, nThreads, iElementBegin, iElementEnd);
		typename TP_Sequential_ST::Prereq p(verbose, tp);
		return p;
	}

	struct TP_Sequential_ST_CDP_ST {
		struct Prereq {
			int verbose;
			typename TP::Prereq tp;

			Prereq(int verbose, typename TP::Prereq const& tp) : verbose(verbose), tp(tp) {}
		} p;

		TP_Sequential_ST_CDP_ST(Prereq const& req) : p(req) {}

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

	static typename TP_Sequential_ST_CDP_ST::Prereq prereq_TP_Sequential_ST_CDP_ST(size_t nThreads, size_t iElementBegin, size_t iElementEnd, int verbose) {
		typename TP::Prereq tp(verbose + 1, nThreads, iElementBegin, iElementEnd);
		typename TP_Sequential_ST_CDP_ST::Prereq p(verbose, tp);
		return p;
	}

	struct NNI {
		struct Prereq {
			int verbose;

			Prereq(int verbose) : verbose(verbose) {}
		} p;

		NNI(Prereq const& req) : p(req) {}

		Tree operator()(Color& color, ThreadPool3& threadpool, Tree& tree) {
			NNIAlg nni(color, threadpool, tree, p.verbose);
			nni.log() << "Performing NNI search..." << endl;
			nni.performNNI();
			Log vlog(p.verbose + 1);
			tree.displaySimpleNewick(vlog.log()) << endl;
			return tree;
		}
	};

	using TP3_NNI = TPx_X<TP3, NNI>;

	using TP3_Sequential_NNI = TPx_Sequential_X<TP3, NNI>;

	static typename TP3_NNI::Prereq prereq_TP3_NNI(size_t nThreads, size_t iElementBegin, size_t iElementEnd, int verbose) {
		typename TP3::Prereq tp(verbose + 1, nThreads, iElementBegin, iElementEnd);
		typename TP3_NNI::Prereq p(verbose, tp);
		return p;
	}

	static typename TP3_Sequential_NNI::Prereq prereq_TP3_Sequential_NNI(size_t nThreads, size_t iElementBegin, size_t iElementEnd, int verbose) {
		typename TP3::Prereq tp(verbose + 1, nThreads, iElementBegin, iElementEnd);
		typename TP3_Sequential_NNI::Prereq p(verbose, tp);
		return p;
	}

	struct ST_NNI_ST {
		struct Prereq {
			int verbose;

			Prereq(int verbose) : verbose(verbose) {}
		} p;

		ST_NNI_ST(Prereq const& req) : p(req) {}

		std::array<Tree, 2> operator()(Color& color, ThreadPool& threadpool, Tree& tree) {
			typename ST::Prereq stp(p.verbose);
			typename NNI::Prereq nnip(p.verbose);
			ST st(stp);
			Tree tree1 = st(color, threadpool, tree);
			Tree tree2 = tree1;
			NNI nni(nnip);
			tree2 = nni(color, threadpool, tree2);
			ST st2(stp);
			tree2 = st2(color, threadpool, tree2);
			return { tree1, tree2 };
		}
	};

	static Tree subsampleProcedure(ConstrainedDP &cdp, Random &random, Color& color, Data const& data, size_t nTaxa, size_t nTrees, size_t nThreads, int verbose) {
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
				if (ARG.has("no-two-step") || nTaxa < 200) {
					typename Parallel_TP_RP::Prereq p = prereq_Parallel_TP_RP(random, nThreads, iElementBegin, iElementEnd, iTreeEnd2 - iTreeBegin2, nTaxa, verbose);
					Parallel_TP_RP job(p);
					for (Tree& tree : job(color)) trees.push_back(tree);
				}
				else {
					typename Parallel_TSP::Prereq p = prereq_Parallel_TSP(random, nThreads, iElementBegin, iElementEnd, iTreeEnd2 - iTreeBegin2, nTaxa, verbose);
					Parallel_TSP job(p);
					for (Tree& tree : job(color)) trees.push_back(tree);
				}
			}
		}
		typename TP_Sequential_ST::Prereq tp_sq_stp = prereq_TP_Sequential_ST(nThreads, 0, data.nElements, verbose);
		typename TP3_Sequential_NNI::Prereq tp3_sq_nnip = prereq_TP3_Sequential_NNI(nThreads, 0, data.nElements, verbose);
		TP_Sequential_ST tp_sq_st(tp_sq_stp);
		vector<Tree> trees2 = tp_sq_st(color, trees);
		for (Tree& tree : trees2) cdp.addTree(tree);
		TP3_Sequential_NNI tp3_sq_nni(tp3_sq_nnip);
		vector<Tree> trees3 = tp3_sq_nni(color, trees2);
		TP_Sequential_ST tp_sq_st2(tp_sq_stp);
		vector<Tree> trees4 = tp_sq_st2(color, trees3);
		for (Tree& tree : trees4) cdp.addTree(tree);
		Tree tree = cdp.optimalUnrootedTree(ZERO);
		typename TP_ST::Prereq tp_stp = prereq_TP_ST(nThreads, 0, data.nElements, verbose);
		TP_ST tp_st(tp_stp);
		Tree tree2 = tp_st(color, tree);
		typename TP3_NNI::Prereq tp3_nnip = prereq_TP3_NNI(nThreads, 0, data.nElements, verbose);
		TP3_NNI tp3_nni(tp3_nnip);
		Tree tree3 = tp3_nni(color, tree2);
		TP_ST tp_st2(tp_stp);
		Tree tree4 = tp_st2(color, tree2);
		return tree4;
	}

	static Tree subsampleProcedure2(ConstrainedDP& cdp, Random& random, Color& color, Data const& data, size_t nTaxa, size_t nTrees, size_t nThreads, int verbose) {
		// Log::setShowThread();
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
		for (size_t iRun : std::views::iota((size_t)0, nRuns)) {
			size_t iTreeBegin = iRun * nTrees / nRuns;
			size_t iTreeEnd = (iRun + 1) * nTrees / nRuns;
			size_t nPartitions = iTreeEnd - iTreeBegin;
			for (size_t iPartition : std::views::iota((size_t)0, nPartitions)) {
				size_t iElementBegin = iPartition * nElements / nPartitions;
				size_t iElementEnd = (iPartition + 1) * nElements / nPartitions;
				//cerr << "iElementBegin = " << iElementBegin << endl;
				//cerr << "iElementEnd = " << iElementEnd << endl;
				Tree tree;
				if (ARG.has("no-two-step") || nTaxa < 200) {
					typename TP_RP::Prereq p = prereq_TP_RP(random, nThreads, iElementBegin, iElementEnd, nTaxa, verbose + 1);
					TP_RP job(p);
					tree = job(color);
				}
				else {
					typename TSP::Prereq p = prereq_TSP(random, nThreads, iElementBegin, iElementEnd, nTaxa, verbose + 1);
					TSP job(p);
					tree = job(color);
				}
				typename TP_ST::Prereq tp_stp = prereq_TP_ST(nThreads, 0, data.nElements, verbose + 1);
				typename TP3_NNI::Prereq tp3_nnip = prereq_TP3_NNI(nThreads, 0, data.nElements, verbose + 1);
				TP_ST tp_st(tp_stp);
				Tree tree2 = tp_st(color, tree);
				cdp.addTree(tree2);
				TP3_NNI tp3_nni(tp3_nnip);
				Tree tree3 = tp3_nni(color, tree2);
				TP_ST tp_st2(tp_stp);
				Tree tree4 = tp_st2(color, tree3);
				cdp.addTree(tree4);
			}
		}
		{
			Tree tree = cdp.optimalUnrootedTree(ZERO);
			typename TP_ST::Prereq tp_stp = prereq_TP_ST(nThreads, 0, data.nElements, verbose);
			TP_ST tp_st(tp_stp);
			Tree tree2 = tp_st(color, tree);
			typename TP3_NNI::Prereq tp3_nnip = prereq_TP3_NNI(nThreads, 0, data.nElements, verbose);
			TP3_NNI tp3_nni(tp3_nnip);
			Tree tree3 = tp3_nni(color, tree2);
			TP_ST tp_st2(tp_stp);
			Tree tree4 = tp_st2(color, tree2);
			return tree4;
		}
	}

	static Tree defaultProcedure(ConstrainedDP& cdp, Random& random, Color& color, Data const& data, size_t nTaxa, size_t nTrees, size_t nThreads, int verbose) {
		return subsampleProcedure2(cdp, random, color, data, nTaxa, nTrees, nThreads, verbose);
	}

	static Tree heuristSearch(Data const& data, size_t nTaxa, size_t r, size_t s, size_t nThreads, int verbose, Tree f(ConstrainedDP&, Random&, Color&, Data const&, size_t, size_t, size_t, int)) {
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
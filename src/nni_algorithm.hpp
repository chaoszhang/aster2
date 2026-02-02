#ifndef NNI_ALGORITHM_HPP
#define NNI_ALGORITHM_HPP

#include "common.hpp"
#include "stepwise_colorable.hpp"
#include "threadpool.hpp"

namespace nni_algorithm {

using std::size_t;
using std::string;
using std::vector;
using std::array;
using std::cerr;
using std::endl;

template<class T> concept QUADRIPARTITION_SCORE_ALGORITHM = requires(T t, typename T::Color & color, typename T::ThreadPool & threadPool, common::AnnotatedBinaryTree & tree, std::size_t index, typename T::score_t score, void annotator(common::AnnotatedBinaryTree::Node*, typename T::score_t, typename T::score_t, typename T::score_t), int verbose) {
	requires std::derived_from<T, common::LogInfo>;
	requires std::same_as<array<typename T::score_t, 3>, typename T::ThreadPool::score_t>;
	requires thread_pool::THREAD_POOL<typename T::ThreadPool>;
	{ T{ color, threadPool, tree, verbose } };
	{ t.labelTree(annotator) };
	{ t.labelTree() };
};

template<class T> concept NNI_ALGORITHM = requires {
	requires QUADRIPARTITION_SCORE_ALGORITHM<T>;
	requires std::same_as<typename T::score_t, typename T::Color::score_t>;
	requires std::integral<typename T::score_t> || std::floating_point<typename T::score_t>;
};

template<class Attributes> concept STEPWISE_COLOR_QUADRIPARTITION_SCORE_ATTRIBUTES = requires(typename Attributes::score_t score){
	requires stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE<typename Attributes::Color>;
	requires thread_pool::SCHEDULER<Attributes::template Scheduler, typename Attributes::Color::score_t>;
	{ score + score } noexcept -> std::convertible_to<typename Attributes::score_t const>;
	{ score - score } noexcept -> std::convertible_to<typename Attributes::score_t const>;
	{ Attributes::ZERO } noexcept -> std::convertible_to<typename Attributes::score_t const>;
};

template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE SC, typename Stats> struct StepwiseColorQuadripartitionScoreDefaultAttributes {
	using Color = SC;
	using score_t = Stats;
	template<typename T> using Scheduler = thread_pool::SimpleScheduler<T>;
	static inline score_t const ZERO = score_t::ZERO;
};

template<class Attributes> concept STEPWISE_COLOR_NNI_ATTRIBUTES = requires{
	requires STEPWISE_COLOR_QUADRIPARTITION_SCORE_ATTRIBUTES<Attributes>;
	requires std::same_as<typename Attributes::score_t, typename Attributes::Color::score_t>;
};

template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE SC> struct StepwiseColorNNIDefaultAttributes {
	using Color = SC;
	using score_t = Color::score_t;
	template<typename T> using Scheduler = thread_pool::SimpleScheduler<T>;
	static inline score_t constexpr ZERO = Color::ZERO;
};

ChangeLog logStepwiseColorQuadripartitionScore("StepwiseColorQuadripartitionScore",
	"2026-02-02", "Chao Zhang", "Initial code", "minor");

template<STEPWISE_COLOR_QUADRIPARTITION_SCORE_ATTRIBUTES Attributes> class StepwiseColorQuadripartitionScore : public common::LogInfo {
public:
	using Color = Attributes::Color;
	using score_t = Attributes::score_t;
	using ThreadPool = thread_pool::ThreadPool<array<score_t, 3>, Attributes::template Scheduler>;
	using Tree = common::AnnotatedBinaryTree;
	static inline string const LEFT_RIGHT = Tree::QUADRIPARTITION_SCORE;
	static inline string const LEFT_OUTGROUP = Tree::QUADRIPARTITION_ALTERNATIVE_1_SCORE;
	static inline string const RIGHT_OUTGROUP = Tree::QUADRIPARTITION_ALTERNATIVE_2_SCORE;
	static inline string const LEAF_ID = Tree::LEAF_ID;
	static inline score_t const ZERO = Attributes::ZERO;

protected:
	Color& stepwiseColor;
	ThreadPool& threadPool;
	Tree& tree;
	thread_pool::Instruction<array<score_t, 3> > instr;
	std::queue<array<score_t, 3> > results;
	std::unordered_map<Tree::Node*, size_t> leafId;

	void setColor(vector<size_t> const& leaves, size_t iColor) noexcept {
		Color& localStepwiseColor = stepwiseColor;
		instr.mapFuncs.emplace_back([leaves, iColor, &localStepwiseColor](size_t iElement) noexcept {
			for (size_t iTaxon : leaves) localStepwiseColor.elementSetTaxonColor(iElement, iTaxon, iColor);
		});
	}

	void clearColor(vector<size_t> const& leaves, size_t iColor) noexcept {
		Color& localStepwiseColor = stepwiseColor;
		instr.mapFuncs.emplace_back([leaves, iColor, &localStepwiseColor](size_t iElement) noexcept {
			for (size_t iTaxon : leaves) localStepwiseColor.elementClearTaxonColor(iElement, iTaxon, iColor);
		});
	}

	void clearSetColor(vector<size_t> const& leaves, size_t iColor, size_t jColor) noexcept {
		Color& localStepwiseColor = stepwiseColor;
		instr.mapFuncs.emplace_back([leaves, iColor, jColor, &localStepwiseColor](size_t iElement) noexcept {
			for (size_t iTaxon : leaves) {
				localStepwiseColor.elementClearTaxonColor(iElement, iTaxon, iColor);
				localStepwiseColor.elementSetTaxonColor(iElement, iTaxon, jColor);
			}
		});
	}

	void computeScore() noexcept {
		Color& localStepwiseColor = stepwiseColor;
		std::function<array<score_t, 3>(size_t)> func = [&localStepwiseColor](size_t iElement) noexcept -> array<score_t, 3> { return localStepwiseColor.elementQuadripartitionScores(iElement); };
		instr.mapFuncs.emplace_back(func);
		instr.reduceFuncs.emplace_back([](array<score_t, 3> a, array<score_t, 3> b) noexcept -> array<score_t, 3> { return { a[0] + b[0], a[1] + b[1], a[2] + b[2] }; });
		instr.zeros.push_back({ZERO , ZERO , ZERO});
	}

	array<score_t, 3> getComputedScore() noexcept {
		array<score_t, 3> res = results.front();
		results.pop();
		return res;
	}

	void leaves(Tree::Node* node, vector<size_t>& leafVec) noexcept {
		if (node->isLeaf()) leafVec.push_back(leafId.at(node));
		else {
			leaves(node->leftChild(), leafVec);
			leaves(node->rightChild(), leafVec);
		}
	}

	vector<size_t> leaves(Tree::Node* node) noexcept {
		vector<size_t> leafVec;
		leaves(node, leafVec);
		return leafVec;
	}

	void runCachedInstructions() {
		for (array<score_t, 3> result : threadPool(instr)) results.push(result);
		instr.mapFuncs.clear();
		instr.reduceFuncs.clear();
		instr.zeros.clear();
	}

	template<bool runSetter, bool runGetter> void colorSubtree(Tree::Node* node, void annotator(Tree::Node*, score_t, score_t, score_t), score_t EPSILON) noexcept {
		if (node->isLeaf()) {
			if constexpr (runSetter) clearSetColor(leaves(node), 0, 2);
			return;
		}

		vector<size_t> rcLeaves, rcrcLeaves, lcrcLeaves;

		rcLeaves = leaves(node->rightChild());
		
		// a:0 b:0 c:0 d:0
		colorSubtree<runSetter, runGetter>(node->rightChild(), annotator, EPSILON); // a:0 b:0 c:2 d:2

		if constexpr (runSetter) clearSetColor(rcLeaves, 2, 0); // a:0 b:0 c:0 d:0

		colorSubtree<runSetter, runGetter>(node->leftChild(), annotator, EPSILON); // a:2 b:2 c:0 d:0

		// Never move the following lines above colorSubtrees! lc & rc often changes!
		Tree::Node* lc = node->leftChild();
		Tree::Node* rc = node->rightChild();
		if (!rc->isLeaf()) rcrcLeaves = leaves(rc->rightChild());
		if (!lc->isLeaf()) lcrcLeaves = leaves(lc->rightChild());

		if constexpr (runSetter) {
			clearSetColor(rcLeaves, 0, 1); // a:2 b:2 c:1 d:1
			clearSetColor(rcrcLeaves, 1, 3); // a:2 b:2 c:1 d:3
			if (!rc->isLeaf()) computeScore(); // (o-c|ab-d, o-ab|c-d, o-d|c-ab)
			clearSetColor(rcrcLeaves, 3, 1); // a:2 b:2 c:1 d:1
			clearSetColor(lcrcLeaves, 2, 3); // a:2 b:3 c:1 d:1
			if (!lc->isLeaf()) computeScore(); // (o-cd|a-b, o-a|cd-b, o-b|cd-a)
			clearSetColor(lcrcLeaves, 3, 2); // a:2 b:2 c:1 d:1
			clearSetColor(rcLeaves, 1, 2); // a:2 b:2 c:2 d:2
		}

		if constexpr (runSetter && runGetter) runCachedInstructions(); // For NNI only

		array<score_t, 3> lcScores = { ZERO, ZERO, ZERO }, rcScores = { ZERO, ZERO, ZERO };
		if constexpr (runGetter) {
			if (!rc->isLeaf()) {
				rcScores = getComputedScore(); // (o-c|ab-d, o-ab|c-d, o-d|c-ab)
				auto [oc, cd, od] = rcScores;
				annotator(rc, cd, oc, od);
			}
			if (!lc->isLeaf()) {
				lcScores = getComputedScore(); // (o-cd|a-b, o-a|cd-b, o-b|cd-a)
				auto [ab, oa, ob] = lcScores;
				annotator(lc, ab, oa, ob);
			}
		}

		// For NNI only
		if constexpr (runSetter && runGetter) nni(node, lcScores, rcScores, EPSILON);
	}

	int nniCounter = 0;
	void nni(Tree::Node* node, array<score_t, 3> lcScores, array<score_t, 3> rcScores, score_t EPSILON) noexcept {
		auto [oc, cd, od] = rcScores; // (o-c|ab-d, o-ab|c-d, o-d|c-ab)
		auto [ab, oa, ob] = lcScores; // (o-cd|a-b, o-a|cd-b, o-b|cd-a)

		Tree::Node* lc = node->leftChild();
		Tree::Node* rc = node->rightChild();

		score_t maxDelta = ZERO;
		Tree::Node* newOutgroup = nullptr;
		if (oa > ab + maxDelta + EPSILON) {
			maxDelta = oa - ab;
			newOutgroup = lc->leftChild();
		}
		if (ob > ab + maxDelta + EPSILON) {
			maxDelta = ob - ab;
			newOutgroup = lc->rightChild();
		}
		if (oc > cd + maxDelta + EPSILON) {
			maxDelta = oc - cd;
			newOutgroup = rc->leftChild();
		}
		if (od > cd + maxDelta + EPSILON) {
			maxDelta = od - cd;
			newOutgroup = rc->rightChild();
		}

		if (newOutgroup) {
			newOutgroup->nni();
			Tree::Node* newSubtreeRoot = newOutgroup->parent();
			node->makeLeftHeavyByLeafCount();
			clearSetColor(leaves(newOutgroup), 2, 0); // a:2 b:2 c:2 d:2 -> a:0 b:2 c:2 d:2 or alike
			nniColor(node, EPSILON);
			clearSetColor(leaves(newOutgroup), 0, 2); // a:0 b:2 c:2 d:2 or alike -> a:0 b:2 c:2 d:2
			newSubtreeRoot->makeLeftHeavyByLeafCount();
			nniColor(newSubtreeRoot, EPSILON);
		}
	}

	void nniColor(Tree::Node* node, score_t EPSILON) {
		if (node->isLeaf()) return;

		Tree::Node* lc = node->leftChild();
		Tree::Node* rc = node->rightChild();

		vector<size_t> rcLeaves, rcrcLeaves, lcrcLeaves;

		rcLeaves = leaves(rc);
		if (!rc->isLeaf()) rcrcLeaves = leaves(rc->rightChild());
		if (!lc->isLeaf()) lcrcLeaves = leaves(lc->rightChild());

		//a:2 b:2 c:2 d:2
		clearSetColor(rcLeaves, 2, 1); // a:2 b:2 c:1 d:1
		clearSetColor(rcrcLeaves, 1, 3); // a:2 b:2 c:1 d:3
		if (!rc->isLeaf()) computeScore(); // (o-c|ab-d, o-ab|c-d, o-d|c-ab)
		clearSetColor(rcrcLeaves, 3, 1); // a:2 b:2 c:1 d:1
		clearSetColor(lcrcLeaves, 2, 3); // a:2 b:3 c:1 d:1
		if (!lc->isLeaf()) computeScore(); // (o-cd|a-b, o-a|cd-b, o-b|cd-a)
		clearSetColor(lcrcLeaves, 3, 2); // a:2 b:2 c:1 d:1
		clearSetColor(rcLeaves, 1, 2); // a:2 b:2 c:2 d:2

		runCachedInstructions();
		
		array<score_t, 3> lcScores = { ZERO, ZERO, ZERO }, rcScores = { ZERO, ZERO, ZERO };
		if (!rc->isLeaf()) rcScores = getComputedScore(); // (o-c|ab-d, o-ab|c-d, o-d|c-ab)
		if (!lc->isLeaf()) lcScores = getComputedScore(); // (o-cd|a-b, o-a|cd-b, o-b|cd-a)

		nni(node, lcScores, rcScores, EPSILON);
	}

	template<bool NNI> void colorTree(void annotator(Tree::Node*, score_t, score_t, score_t), score_t EPSILON = ZERO) noexcept {
		vector<size_t> leaves;
		for (auto const& e : leafId) leaves.push_back(e.second);
		setColor(leaves, 0);
		if constexpr (NNI) colorSubtree<true, true>(tree.root(), annotator, EPSILON);
		else colorSubtree<true, false>(tree.root(), annotator, EPSILON);
		clearColor(leaves, 2);
		runCachedInstructions();
		if constexpr (!NNI) colorSubtree<false, true>(tree.root(), annotator, EPSILON);
	}

	static void dummyAnnotator(Tree::Node*, score_t, score_t, score_t) noexcept {}

	static void defaultAnnotator(Tree::Node* node, score_t left_right, score_t left_outgroup, score_t right_outgroup) {
		node->set(LEFT_RIGHT, left_right);
		node->set(LEFT_OUTGROUP, left_outgroup);
		node->set(RIGHT_OUTGROUP, right_outgroup);
	}

public:
	StepwiseColorQuadripartitionScore(Color& stepwiseColor, ThreadPool& threadPool, Tree& tree, int verbose = common::LogInfo::DEFAULT_VERBOSE) : stepwiseColor(stepwiseColor), threadPool(threadPool), tree(tree), LogInfo(verbose) {
		tree.makeLeftHeavyByLeafCount();
		for (Tree::Node* node : tree.leaves()) leafId[node] = std::any_cast<size_t>(node->get(LEAF_ID));
	}

	void labelTree(void annotator(Tree::Node*, score_t, score_t, score_t) = defaultAnnotator) noexcept { colorTree<false>(annotator); }
};

ChangeLog logStepwiseColorNNI("StepwiseColorNNI",
	"2026-02-02", "Chao Zhang", "Initial code", "minor");

template<STEPWISE_COLOR_NNI_ATTRIBUTES Attributes> class StepwiseColorNNI : public StepwiseColorQuadripartitionScore<Attributes> {
public:
	using Color = Attributes::Color;
	using score_t = Color::score_t;
	using ThreadPool = thread_pool::ThreadPool<array<score_t, 3>, Attributes::template Scheduler>;
	using Tree = common::AnnotatedBinaryTree;
	static inline score_t const ZERO = Color::ZERO;
	static inline score_t constexpr EPSILON = Color::EPSILON;

public:
	StepwiseColorNNI(Color& stepwiseColor, ThreadPool& threadPool, Tree& tree, int verbose = common::LogInfo::DEFAULT_VERBOSE) : StepwiseColorQuadripartitionScore<Attributes>(stepwiseColor, threadPool, tree, verbose) {}

	void performNNI() noexcept { this->template colorTree<true>(this->dummyAnnotator, EPSILON); }
};

};
#endif
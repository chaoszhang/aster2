#ifndef PLACEMENT_ALGORITHM_HPP
#define PLACEMENT_ALGORITHM_HPP

#include "common.hpp"
#include "stepwise_colorable.hpp"
#include "threadpool.hpp"

namespace placement_algorithm{

using std::size_t;
using std::string;
using std::vector;

template<class T> concept PLACEMENT_ALGORITHM = requires(T t, typename T::Color & color, typename T::ThreadPool & threadPool, common::AnnotatedBinaryTree & tree, std::size_t index, typename T::score_t score, int verbose) {
	requires std::derived_from<T, common::LogInfo>;
	requires std::integral<typename T::score_t> || std::floating_point<typename T::score_t>;
	requires std::same_as<typename T::score_t, typename T::Color::score_t>;
	requires std::same_as<typename T::score_t, typename T::ThreadPool::score_t>;
	requires thread_pool::THREAD_POOL<typename T::ThreadPool>;
	{ T{ color, threadPool, tree, verbose } };
	{ t.scoreTree() } noexcept -> std::same_as<typename T::score_t>;
	{ t.optimalPlacement(index) } noexcept -> std::same_as<typename T::Tree::Node*>;
	{ t.place(index) } noexcept;
};

template<class Attributes> concept STEPWISE_COLOR_PLACEMENT_ATTRIBUTES = requires{
	requires stepwise_colorable::STEPWISE_COLORABLE<typename Attributes::Color>;
	requires thread_pool::SCHEDULER<Attributes::template Scheduler, typename Attributes::Color::score_t>;
};

template<stepwise_colorable::STEPWISE_COLORABLE SC> struct StepwiseColorPlacementDefaultAttributes{
	using Color = SC;
	template<typename T> using Scheduler = thread_pool::SimpleScheduler<T>;
};

template<STEPWISE_COLOR_PLACEMENT_ATTRIBUTES Attributes> class StepwiseColorPlacement : public common::LogInfo {
public:
	using Color = Attributes::Color;
	using score_t = Color::score_t;
	using ThreadPool = thread_pool::ThreadPool<score_t, Attributes::template Scheduler>;
	using Tree = common::AnnotatedBinaryTree;
	static inline string const NO_PLACE = Tree::TRIPARTITION_SCORE;
	static inline string const PLACE_ABOVE = "PLACE_ABOVE";
	static inline string const PLACE_LEFT = "PLACE_LEFT";
	static inline string const PLACE_RIGHT = "PLACE_RIGHT";
	static inline string const NEW_NODE = "NEW_NODE";
	static inline string const SCORE_OVER = "SCORE_OVER";
	static inline string const SCORE_WITHIN = "SCORE_WITHIN";
	static inline string const LEAF_ID = Tree::LEAF_ID;
	static inline score_t constexpr ZERO = Color::ZERO;

private:
	template<bool isSetter> using SetterOnlyVector = typename std::conditional<isSetter, vector<size_t>, bool>::type;

	Color &stepwiseColor;
	ThreadPool &threadPool;
	Tree &tree;
	thread_pool::Instruction<score_t> instr;
	std::queue<score_t> results;
	std::unordered_map<Tree::Node*, size_t> leafId;
	
	template<bool isSetter> void setColor(SetterOnlyVector<isSetter> const &leaves, size_t iColor) noexcept{
		if constexpr (isSetter){
			Color &localStepwiseColor = stepwiseColor;
			instr.mapFuncs.emplace_back([leaves, iColor, &localStepwiseColor](size_t iElement) noexcept{
				for (size_t iTaxon: leaves) localStepwiseColor.elementSetTaxonColor(iElement, iTaxon, iColor);
			});
		}
	}
	
	template<bool isSetter> void clearColor(SetterOnlyVector<isSetter> const &leaves, size_t iColor) noexcept{
		if constexpr (isSetter){
			Color &localStepwiseColor = stepwiseColor;
			instr.mapFuncs.emplace_back([leaves, iColor, &localStepwiseColor](size_t iElement) noexcept{
				for (size_t iTaxon: leaves) localStepwiseColor.elementClearTaxonColor(iElement, iTaxon, iColor);
			});
		}
	}
	
	template<bool isSetter> void clearSetColor(SetterOnlyVector<isSetter> const &leaves, size_t iColor, size_t jColor) noexcept{
		if constexpr (isSetter){
			Color &localStepwiseColor = stepwiseColor;
			instr.mapFuncs.emplace_back([leaves, iColor, jColor, &localStepwiseColor](size_t iElement) noexcept {
				for (size_t iTaxon: leaves) {
					localStepwiseColor.elementClearTaxonColor(iElement, iTaxon, iColor);
					localStepwiseColor.elementSetTaxonColor(iElement, iTaxon, jColor);
				}
			});
		}
	}
	
	template<bool isSetter> score_t computeScore() noexcept{
		Color &localStepwiseColor = stepwiseColor;
		if constexpr (isSetter) {
			std::function<score_t(size_t)> func = [&localStepwiseColor](size_t iElement) noexcept -> score_t{ return localStepwiseColor.elementScore(iElement); };
			instr.mapFuncs.emplace_back(func);
			instr.reduceFuncs.emplace_back([](score_t a, score_t b) noexcept -> score_t{ return a + b; });
			instr.zeros.push_back(ZERO);
			return ZERO;
		}
		else {
			score_t res = results.front();
			results.pop();
			return res;
		}
	}
	
	template<bool isSetter, bool isPlacement> SetterOnlyVector<isSetter> colorSubtree(Tree::Node* node, SetterOnlyVector<isSetter> const& newTaxon) noexcept{
		SetterOnlyVector<isSetter> leaves;
		score_t score;
		if (node->isLeaf()){
			if constexpr(isSetter) leaves.push_back(leafId.at(node));
			clearSetColor<isSetter>(leaves, 0, 2);
		}
		else {
			SetterOnlyVector<isSetter> rcLeaves = colorSubtree<isSetter, isPlacement>(node->rightChild(), newTaxon);
			clearSetColor<isSetter>(rcLeaves, 2, 0);
			leaves = colorSubtree<isSetter, isPlacement>(node->leftChild(), newTaxon);
			if constexpr(isSetter) for (size_t iTaxon : rcLeaves) leaves.push_back(iTaxon); 
			clearSetColor<isSetter>(rcLeaves, 0, 1);
			if constexpr(isPlacement) {
				clearSetColor<isSetter>(newTaxon, 0, 1);
				score = computeScore<isSetter>();
				if constexpr(!isSetter) node->set(PLACE_RIGHT, score);
				clearSetColor<isSetter>(newTaxon, 1, 2);
				score = computeScore<isSetter>();
				if constexpr (!isSetter) node->set(PLACE_LEFT, score);
				clearSetColor<isSetter>(newTaxon, 2, 0);
				score = computeScore<isSetter>();
				if constexpr (!isSetter) node->set(PLACE_ABOVE, score);
			}
			else {
				score = computeScore<isSetter>();
				if constexpr(!isSetter) node->set(NO_PLACE, score);
			}
			clearSetColor<isSetter>(rcLeaves, 1, 2);
		}
		if constexpr(isPlacement) {
			clearSetColor<isSetter>(newTaxon, 0, 1);
			score = computeScore<isSetter>();
			if constexpr(!isSetter) node->set(NEW_NODE, score);
			clearSetColor<isSetter>(newTaxon, 1, 0);
		}

		return leaves;
	}
	
	template<bool isPlacement> void labelTree(size_t iNewTaxon = -1) noexcept{
		vector<size_t> leaves, newTaxon;
		if constexpr (isPlacement) newTaxon.push_back(iNewTaxon);
		for (auto const &e : leafId) leaves.push_back(e.second);
		setColor<true>(leaves, 0);
		setColor<true>(newTaxon, 0);
		colorSubtree<true, isPlacement>(tree.root(), newTaxon);
		clearColor<true>(leaves, 2);
		clearColor<true>(newTaxon, 0);
		for (score_t result: threadPool(instr)) results.push(result);
		colorSubtree<false, isPlacement>(tree.root(), true);
	}
	
	score_t scoreSubtree(Tree::Node* node) noexcept{
		if (node->isLeaf()) return ZERO;
		return scoreSubtree(node->leftChild()) + scoreSubtree(node->rightChild()) + node->get<score_t>(NO_PLACE);
	}

	std::tuple<score_t, score_t, Tree::Node*> optimalPlacementSubtree(Tree::Node* node) noexcept{
		score_t newNodeScore = node->get<score_t>(NEW_NODE);
		if (node->isLeaf()) return std::make_tuple(ZERO, newNodeScore, node);

		std::tuple<score_t, score_t, Tree::Node*> left = optimalPlacementSubtree(node->leftChild());
		std::tuple<score_t, score_t, Tree::Node*> right = optimalPlacementSubtree(node->rightChild());
		score_t placeLeftScore = node->get<score_t>(PLACE_LEFT);
		score_t placeRightScore = node->get<score_t>(PLACE_RIGHT);
		score_t placeAboveScore = node->get<score_t>(PLACE_ABOVE);

		score_t dpOver = std::get<0>(left) + std::get<0>(right) + placeAboveScore;
		score_t dpAbove = dpOver + newNodeScore;
		score_t dpLeft = std::get<1>(left) + std::get<0>(right) + placeLeftScore;
		score_t dpRight = std::get<0>(left) + std::get<1>(right) + placeRightScore;
		score_t dpMax = std::max({ dpAbove, dpLeft, dpRight });
		return std::make_tuple(dpOver, dpMax, (dpAbove == dpMax) ? node : (dpLeft == dpMax) ? std::get<2>(left) : std::get<2>(right));
	}

public:
	StepwiseColorPlacement(Color& stepwiseColor, ThreadPool& threadPool, Tree& tree, int verbose = common::LogInfo::DEFAULT_VERBOSE) : stepwiseColor(stepwiseColor), threadPool(threadPool), tree(tree), LogInfo(verbose) {
		tree.makeLeafHeavyByLeafCount();
		for (Tree::Node* node : tree.leaves()) leafId[node] = std::any_cast<size_t>(node->get(LEAF_ID));
	}

	score_t scoreTree() noexcept {
		labelTree<false>();
		return scoreSubtree(tree.root());
	}

	Tree::Node* optimalPlacement(size_t iNewTaxon) noexcept {
		labelTree<true>(iNewTaxon);
		return std::get<2>(optimalPlacementSubtree(tree.root()));
	}

	void place(size_t iNewTaxon) noexcept {
		optimalPlacement(iNewTaxon)->emplaceAbove(LEAF_ID, iNewTaxon);
	}
};


};
#endif
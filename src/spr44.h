#ifndef SPR44_H
#define SPR44_H

#include "common.hpp"
#include "stepwise_colorable.hpp"
#include "threadpool.hpp"
#include "astral.hpp"
#include "caster.hpp"
#include <set>
#include <map>
#include <algorithm>

namespace spr44_algorithm {

	using std::size_t;
	using std::string;
	using std::vector;
	using std::array;
	using std::cerr;
	using std::endl;

	using std::tuple;
	using std::make_tuple;
	using std::get;
	using std::unordered_map;
	using std::map;

	template<class T> concept QUADRIPARTITION_SCORE_ALGORITHM = requires(T t, typename T::Color & color, typename T::ThreadPool & threadPool, common::AnnotatedBinaryTree & tree, std::size_t index, typename T::score_t score, void annotator(common::AnnotatedBinaryTree::Node*, typename T::score_t, typename T::score_t, typename T::score_t), int verbose) {
		requires std::derived_from<T, common::LogInfo>;
		requires std::same_as<array<typename T::score_t, 3>, typename T::ThreadPool::score_t>;
		requires thread_pool::THREAD_POOL<typename T::ThreadPool>;
	{ T{ color, threadPool, tree, verbose } };
	{ t.labelTree() };
	};

	template<class T> concept SPR44_ALGORITHM = requires {
		requires QUADRIPARTITION_SCORE_ALGORITHM<T>;
		requires std::same_as<typename T::score_t, typename T::Color::score_t>;
		requires std::integral<typename T::score_t> || std::floating_point<typename T::score_t>;
	};

	template<class Attributes> concept STEPWISE_COLOR_QUADRIPARTITION_SCORE_ATTRIBUTES = requires(typename Attributes::score_t score) {
		requires stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE<typename Attributes::Color>;
		requires thread_pool::SCHEDULER<Attributes::template Scheduler, typename Attributes::Color::score_t>;
	{ Attributes::ZERO } noexcept -> std::convertible_to<typename Attributes::score_t const>;
	};

	template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE SC, typename Support> struct StepwiseColorQuadripartitionScoreDefaultAttributes_SPR44 {
		using Color = SC;
		using score_t = Support;
		template<typename T> using Scheduler = thread_pool::SimpleScheduler<T>;

		static inline score_t const ZERO = score_t::ZERO;

		static array<score_t, 3> map(Color& color, size_t iElement) { return score_t::map(color, iElement); }

		static score_t reduce(score_t const& a, score_t const& b) { return score_t::reduce(a, b); }

		static void annotate(common::AnnotatedBinaryTree::Node* node, score_t left, score_t right, score_t outgroup, string key) { score_t::annotate(node, left, right, outgroup, key); }
	};

	template<class Attributes> concept STEPWISE_COLOR_SPR44_ATTRIBUTES = requires{
		requires STEPWISE_COLOR_QUADRIPARTITION_SCORE_ATTRIBUTES<Attributes>;
		requires std::same_as<typename Attributes::score_t, typename Attributes::Color::score_t>;
	};

	template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE SC> struct StepwiseColorSPR44DefaultAttributes {
		using Color = SC;
		using score_t = Color::score_t;
		template<typename T> using Scheduler = thread_pool::SimpleScheduler<T>;
		static inline score_t constexpr ZERO = Color::ZERO;

		static array<score_t, 3> map(Color& color, size_t iElement) { return color.elementQuadripartitionScores(iElement); }

		static score_t reduce(score_t const& a, score_t const& b) { return a + b; }

		static void annotate(common::AnnotatedBinaryTree::Node*, score_t, score_t, score_t, string) noexcept {}
	};


	template<STEPWISE_COLOR_QUADRIPARTITION_SCORE_ATTRIBUTES Attributes> class StepwiseColorQuadripartitionScore_SPR44 : public common::LogInfo {
	public:
		using Color = Attributes::Color;
		using score_t = Attributes::score_t;
		using ThreadPool = thread_pool::ThreadPool<array<score_t, 3>, Attributes::template Scheduler>;
		using Tree = common::AnnotatedBinaryTree;
		static inline string const PLACE44 = "PLACE44"; //对于某种切割方式，每一个节点下的四色四叶树 
		static inline string const PlaceAns = "PlaceAns"; //对于某种切割方式，放置在不同节点的得分
		static inline string const LEAF_ID = Tree::LEAF_ID;
		static inline score_t const ZERO = Attributes::ZERO;
		static inline bool const IS_ROOTED = Color::IS_ROOTED;

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
					if constexpr (stepwise_colorable::ELEMENT_CLEAR_AND_SET_TAXON_COLOR<Color>) {
						localStepwiseColor.elementClearAndSetTaxonColor(iElement, iTaxon, iColor, jColor);
					}
					else {
						localStepwiseColor.elementClearTaxonColor(iElement, iTaxon, iColor);
						localStepwiseColor.elementSetTaxonColor(iElement, iTaxon, jColor);
					}
				}
				});
		}

		void computeScore() noexcept {
			Color& localStepwiseColor = stepwiseColor;
			std::function<array<score_t, 3>(size_t)> func = [&localStepwiseColor](size_t iElement) noexcept -> array<score_t, 3> { return Attributes::map(localStepwiseColor, iElement); };
			instr.mapFuncs.emplace_back(func);
			instr.reduceFuncs.emplace_back([](array<score_t, 3> a, array<score_t, 3> b) noexcept -> array<score_t, 3> { return { Attributes::reduce(a[0], b[0]), Attributes::reduce(a[1], b[1]), Attributes::reduce(a[2], b[2]) }; });
			instr.zeros.push_back({ ZERO , ZERO , ZERO });
		}

		array<score_t, 3> getComputedScore() noexcept {
			array<score_t, 3> res = results.front();
			results.pop();
			return res;
		}

		void leaves(Tree::Node* node, vector<size_t>& leafVec) noexcept {
			if (node->get<bool>("prune")) return;
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

		template<bool SPR44> void colorSubtree(Tree::Node* node, void annotator(Tree::Node*, score_t, score_t, score_t, string), string key1) noexcept {
			if (node->get<bool>("prune")) return;
			if (node->isLeaf()) {
				clearSetColor(leaves(node), 0, 2);
				return;
			}

			vector<size_t> rcLeaves = {};
			rcLeaves = leaves(node->rightChild());

			colorSubtree<SPR44>(node->rightChild(), annotator, key1);
			clearSetColor(rcLeaves, 2, 0);
			colorSubtree<SPR44>(node->leftChild(), annotator, key1);

			clearSetColor(rcLeaves, 0, 1);
			computeScore();
			clearSetColor(rcLeaves, 1, 2);

			runCachedInstructions();

			array<score_t, 3> Scores = { ZERO, ZERO, ZERO };
			Scores = getComputedScore(); // (left, right, outgroup)
			auto [left, right, outgroup] = Scores;
			if constexpr (SPR44) node->set(key1, Scores); // For SPR44 only
			annotator(node, left, right, outgroup, key1);
			return;
		}

		template<bool SPR44> void colorTree(void annotator(Tree::Node*, score_t, score_t, score_t, string), vector<Tree::Node*> leafId_Slice, string key1) noexcept {
			if constexpr (stepwise_colorable::SET_QUADRIPARTITION_MODE<Color>) stepwiseColor.setQuadripartitionMode(true);

			vector<size_t> leaves = {};
			for (auto const& e : leafId) {
				leaves.push_back(e.second);
			}
			setColor(leaves, 0);
			vector<size_t> leaves_slice = {};
			for (auto const& e : leafId_Slice) {
				leaves_slice.push_back(e->get<size_t>(LEAF_ID));
			}
			clearSetColor(leaves_slice, 0, 3);

			colorSubtree<SPR44>(tree.root(), annotator, key1);

			clearSetColor(leaves_slice, 3, 2);
			clearColor(leaves, 2);
			runCachedInstructions();

			if constexpr (stepwise_colorable::SET_QUADRIPARTITION_MODE<Color>) stepwiseColor.setQuadripartitionMode(false);
			return;
		}

		static void dummyAnnotator(Tree::Node*, score_t, score_t, score_t, string) noexcept {}

	protected:
		//计算这个切割方式下每一种放置位置的得分 
		vector<tuple<typename Color::score_t, Tree::Node*>> ans = {};
		vector<tuple<typename Color::score_t, Tree::Node*>> optimalPlacementSubtree(Tree::Node* node, int flag, string key1, string key2) noexcept {
			if (node->isRoot()) {
				node->set(key2, make_tuple(Color::ZERO, node));
				ans.push_back(node->get<tuple<typename Color::score_t, Tree::Node*>>(key2));
				if (!node->isLeaf() && !node->get<bool>("prune")) {
					optimalPlacementSubtree(node->leftChild(), 1, key1, key2);
					optimalPlacementSubtree(node->rightChild(), 0, key1, key2);
				}
				return ans;
			}

			if (flag == 1) {
				typename Color::score_t scores = get<0>((node->parent())->get<tuple<typename Color::score_t, Tree::Node*>>(key2)) - node->parent()->get<array<typename Color::score_t, 3>>(key1).at(2) + node->parent()->get<array<typename Color::score_t, 3>>(key1).at(0);
				node->set(key2, make_tuple(scores, node));
				ans.push_back(node->get<tuple<typename Color::score_t, Tree::Node*>>(key2));
			}
			else if (flag == 0) {
				typename Color::score_t scores = get<0>((node->parent())->get<tuple<typename Color::score_t, Tree::Node*>>(key2)) - node->parent()->get<array<typename Color::score_t, 3>>(key1).at(2) + node->parent()->get<array<typename Color::score_t, 3>>(key1).at(1);
				node->set(key2, make_tuple(scores, node));
				ans.push_back(node->get<tuple<typename Color::score_t, Tree::Node*>>(key2));
			}

			if (!node->isLeaf() && !node->get<bool>("prune")) {
				optimalPlacementSubtree(node->leftChild(), 1, key1, key2);
				optimalPlacementSubtree(node->rightChild(), 0, key1, key2);
			}
			return ans;
		}

		//遍历所有放置，找这个切割下的最好 
		template<bool SPR44> tuple<score_t, Tree::Node*> optimalPlacement(Tree& tree, Tree::Node* treeSlice, void annotator(Tree::Node*, score_t, score_t, score_t, string)) noexcept {
			vector<Tree::Node*> leafId_Slice = {};
			treeSlice->subtreeLeaves(leafId_Slice);

			int flag = node_to_num[treeSlice];
			string key1 = PLACE44 + std::to_string(flag);
			string key2 = PlaceAns + std::to_string(flag);

			colorTree<SPR44>(annotator, leafId_Slice, key1);

			if constexpr (SPR44) {
				ans.clear();
				optimalPlacementSubtree(tree.root(), -1, key1, key2);

				score_t standard_ans = 0;
				for (int i = 0;i < ans.size();i++) {
					Tree::Node* node_i = get<1>(ans[i]);
					if (node_i->get<bool>("prune")) {
						standard_ans = get<0>(ans[i]);
						break;
					}
				}

				auto it = std::max_element(
					ans.begin(),
					ans.end(),
					[](auto const& a, auto const& b) {
						return std::get<0>(a) < std::get<0>(b);
					}
				);
				get<0>(*it) = get<0>(*it) - standard_ans;
				return *it;
			}
			return {};
		}

	protected:
		vector<Tree::Node*> node_euler;
		vector<int> depth;
		map<Tree::Node*, int> node_to_num;
		map<int, Tree::Node*> num_to_node;
		int lca[410][410];
		int lg2[410] = { 0 };
		void dfs(Tree::Node* node, int dep, bool cleared) {
			if (!cleared) {
				node_euler.clear();
				node_to_num.clear();
				num_to_node.clear();
				depth.clear();
			}
			node_euler.push_back(node);
			depth.push_back(dep);
			if (!node_to_num.contains(node)) {
				node_to_num[node] = node_euler.size() - 1;
			}
			num_to_node[node_euler.size() - 1] = node;
			if (!node->isLeaf()) {
				dfs(node->leftChild(), dep + 1, 1);
				node_euler.push_back(node);
				num_to_node[node_euler.size() - 1] = node;
				depth.push_back(dep);
				dfs(node->rightChild(), dep + 1, 1);
				node_euler.push_back(node);
				num_to_node[node_euler.size() - 1] = node;
				depth.push_back(dep);
			}
			return;
		}

		void RMQ() {
			int m = node_euler.size();
			lg2[0] = -1;
			for (int i = 1; i <= m; ++i) {
				lg2[i] = lg2[i >> 1] + 1;
			}
			for (int i = 0; i <= m - 1; ++i) {
				lca[i][0] = i;
			}
			for (int j = 1; j <= lg2[m]; ++j) {
				for (int i = 0; i + (1 << j) - 1 <= (m - 1); ++i) {
					lca[i][j] = depth[lca[i][j - 1]] < depth[lca[i + (1 << (j - 1))][j - 1]] ? lca[i][j - 1] : lca[i + (1 << j - 1)][j - 1];
				}
			}
		}

		Tree::Node* getlca(Tree::Node* p, Tree::Node* r) {
			Tree::Node* lca1;
			Tree::Node* lca2;
			if (node_to_num[p] > node_to_num[r]) {
				int p_first = node_to_num[p];
				int r_first = node_to_num[r];
				int d = lg2[p_first - r_first + 1];
				lca1 = node_euler[lca[r_first][d]];
				lca2 = node_euler[lca[p_first - (1 << d) + 1][d]];
			}
			else {
				int p_first = node_to_num[p];
				int r_first = node_to_num[r];
				int d = lg2[r_first - p_first + 1];
				lca1 = node_euler[lca[p_first][d]];
				lca2 = node_euler[lca[r_first - (1 << d) + 1][d]];
			}
			return depth[node_to_num[lca1]] < depth[node_to_num[lca2]] ? lca1 : lca2;
		}

	protected:
		void setPrune(Tree::Node* node) {
			node->set("prune", false);
			if (node->isLeaf()) return;
			setPrune(node->leftChild());
			setPrune(node->rightChild());
			return;
		}

		//遍历所有切割方式，得到最好的切割和组合
		vector<tuple<score_t, Tree::Node*, Tree::Node*, bool>> final_ans;
		template<bool SPR44> vector<tuple<score_t, Tree::Node*, Tree::Node*, bool>> regraft(Tree& tree, Tree::Node* tree_change, void annotator(Tree::Node*, score_t, score_t, score_t, string)) {
			setPrune(tree.root());
			tree_change->set("prune", true);

			tuple<score_t, Tree::Node*> ans1 = optimalPlacement<SPR44>(tree, tree_change, annotator);

			if constexpr (SPR44) {
				Tree::Node* regraft_site = get<1>(ans1);
				bool flagprune = 0;
				if (!(regraft_site->get<bool>("prune"))) {
					if (regraft_site->isLeaf()) {
						flagprune = true;
					}
					else if (!(regraft_site->leftChild()->get<bool>("prune") || regraft_site->rightChild()->get<bool>("prune"))) {
						flagprune = true;
					}
				}
				final_ans.push_back(make_tuple(get<0>(ans1), get<1>(ans1), tree_change, flagprune));

				if (!tree_change->isLeaf()) {
					regraft<SPR44>(tree, tree_change->leftChild(), annotator);
					regraft<SPR44>(tree, tree_change->rightChild(), annotator);
				}
				return final_ans;
			}

			if (!tree_change->isLeaf()) {
				regraft<SPR44>(tree, tree_change->leftChild(), annotator);
				regraft<SPR44>(tree, tree_change->rightChild(), annotator);
			}
			return {};
		}

	protected:
		void set_end_num(Tree::Node* node) {
			node->set("end", false);
			node->set("num", false);
			if (node->isLeaf()) return;
			set_end_num(node->leftChild());
			set_end_num(node->rightChild());
			return;
		}

		void get_num(Tree::Node* node) {
			if (!node->get<bool>("end") && !node->isRoot()) {
				node->set("num", true);
				get_num(node->parent());
			}
			return;
		}

		void single_iBootstrap_support(Tree::Node* node, Tree::Node* treeSlice, size_t iBootstrap) noexcept {
			using score_t = typename Color::score_t;
			string key1 = "PLACE44" + std::to_string(node_to_num[treeSlice]) + std::to_string(iBootstrap);
			string key2 = "PlaceAns" + std::to_string(node_to_num[treeSlice]) + std::to_string(iBootstrap);

			setPrune(node);
			treeSlice->set("prune", true);

			ans.clear();
			ans = optimalPlacementSubtree(node, -1, key1, key2);

			typename Color::score_t standard_ans = 0.0;
			for (int i = 0;i < ans.size();i++) {
				Tree::Node* node_i = get<1>(ans[i]);
				if (node_i->get<bool>("prune")) {
					standard_ans = get<0>(ans[i]);
					break;
				}
			}

			//bootstrap
			for (int i = 0;i < ans.size();i++) {
				Tree::Node* node_i = get<1>(ans[i]);
				if (node_i->get<bool>("prune")) continue;
				if (get<0>(ans[i]) >= standard_ans + Color::EPSILON) {
					Tree::Node* node_regraft = get<1>(ans[i]);
					Tree::Node* lca = getlca(treeSlice, node_regraft);
					lca->set("end", true);
					if (!treeSlice->isRoot()) {
						get_num(treeSlice->parent());
					}
					//calculate the pruned
					if (!treeSlice->parent()->get<bool>("end")) {
						treeSlice->set("num", true);
					}
					//
					if (!node_regraft->isRoot()) {
						get_num(node_regraft);
					}
					lca->set("end", false);
				}
			}

			if (treeSlice->isLeaf()) return;
			single_iBootstrap_support(node, treeSlice->leftChild(), iBootstrap);
			single_iBootstrap_support(node, treeSlice->rightChild(), iBootstrap);
			return;
		}

	protected:
		void change_bool_to_int(Tree::Node* node) {
			if (!node->get<bool>("num")) {
				if (node->has("SPR44BlockBootStrap")) {
					int cnt = node->get<int>("SPR44BlockBootStrap");
					cnt++;
					node->set("SPR44BlockBootStrap", cnt);
				}
				else node->set("SPR44BlockBootStrap", 1);
			}
			else {
				if (!node->has("SPR44BlockBootStrap")) {
					node->set("SPR44BlockBootStrap", 0);
				}
			}
			if (node->isLeaf()) return;
			change_bool_to_int(node->leftChild());
			change_bool_to_int(node->rightChild());
			return;
		}

		void final_support_normalized(Tree::Node* node, size_t N_BOOTSTRAPS) {
			if (node->has<int>("SPR44BlockBootStrap")) {
				int cnt = node->get<int>("SPR44BlockBootStrap");
				node->set("SPR44BlockBootStrap", (double)cnt * 100 / N_BOOTSTRAPS);
			}
			else node->set("SPR44BlockBootStrap", 0);
			if (node->isLeaf()) return;
			final_support_normalized(node->leftChild(), N_BOOTSTRAPS);
			final_support_normalized(node->rightChild(), N_BOOTSTRAPS);
			return;
		}

		void final_support(common::AnnotatedBinaryTree::Node* node, size_t N_BOOTSTRAPS) {
			for (size_t iBootstrap : std::views::iota((size_t)0, N_BOOTSTRAPS)) {
				set_end_num(node);
				single_iBootstrap_support(node, node->leftChild(), iBootstrap);
				single_iBootstrap_support(node, node->rightChild(), iBootstrap);
				change_bool_to_int(node);
			}
			final_support_normalized(node, N_BOOTSTRAPS);
			return;
		}

	protected:
		template<bool SPR44> tuple<score_t, Tree::Node*, Tree::Node*, bool> optimalSpr44(void annotator(Tree::Node*, score_t, score_t, score_t, string), size_t N_BOOTSTRAPS) noexcept {
			dfs(tree.root(), 0, 0);
			final_ans.clear();
			final_ans = regraft<SPR44>(tree, tree.root()->leftChild(), annotator);
			final_ans = regraft<SPR44>(tree, tree.root()->rightChild(), annotator);

			if constexpr (SPR44) {
				auto it = std::max_element(
					final_ans.begin(),
					final_ans.end(),
					[](auto const& a, auto const& b) {
						return std::get<0>(a) < std::get<0>(b);
					}
				);
				return *it;
			}
			else {
				RMQ();
				final_support(tree.root(), N_BOOTSTRAPS);
			}
			return {};
		}

	public:
		StepwiseColorQuadripartitionScore_SPR44(Color& stepwiseColor, ThreadPool& threadPool, Tree& tree, int verbose = common::LogInfo::DEFAULT_VERBOSE) : stepwiseColor(stepwiseColor), threadPool(threadPool), tree(tree), LogInfo(verbose) {
			tree.makeLeftHeavyByLeafCount();
			for (Tree::Node* node : tree.leaves()) leafId[node] = std::any_cast<size_t>(node->get(LEAF_ID));
		}

		void labelTree(size_t N_BOOTSTRAPS) noexcept { optimalSpr44<false>(Attributes::annotate, N_BOOTSTRAPS); }
	};



	template<STEPWISE_COLOR_SPR44_ATTRIBUTES Attributes> class StepwiseColorSPR44 : public StepwiseColorQuadripartitionScore_SPR44<Attributes> {
	public:
		using Color = Attributes::Color;
		using score_t = Color::score_t;
		using ThreadPool = thread_pool::ThreadPool<array<score_t, 3>, Attributes::template Scheduler>;
		using Tree = common::AnnotatedBinaryTree;
		static inline score_t const ZERO = Color::ZERO;
		static inline score_t constexpr EPSILON = Color::EPSILON;

	public:
		StepwiseColorSPR44(Color& stepwiseColor, ThreadPool& threadPool, Tree& tree, int verbose = common::LogInfo::DEFAULT_VERBOSE) : StepwiseColorQuadripartitionScore_SPR44<Attributes>(stepwiseColor, threadPool, tree, verbose) {}

		//实际进行切割和组合 
		void spr_place() noexcept {
			auto result = this->template optimalSpr44<true>(this->dummyAnnotator, 0);
			this->log() << "changeNumber" << get<0>(result) << std::endl;

			Tree::Node* nodeprune = get<2>(result);
			Tree::Node* noderegraft = get<1>(result);

			if (get<3>(result)) {
				Tree sub = nodeprune->pruneAbove();
				noderegraft->regraftAbove(sub);
			}
			return;
		}
	};

};
#endif


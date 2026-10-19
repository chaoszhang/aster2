#ifndef SPR44_HPP
#define SPR44_HPP

#include "common.hpp"
#include "stepwise_colorable.hpp"
#include "threadpool.hpp"

namespace spr44_algorithm {

	using std::size_t;
	using std::string;
	using std::vector;
	using std::array;
	using std::cerr;
	using std::endl;

	using std::tuple;
	using std::pair;
	using std::make_tuple;
	using std::make_pair;
	using std::get;
	using std::unordered_map;
	using std::map;
	using std::shared_ptr;

	ChangeLog logSPR("logSPR",
		"2026-10-16", "Yujie Wang", "add spr", "minor",
		"2026-10-19", "Chao Zhang", "memory efficiency improvements by DeepSeek", "patch");

	template<class T> concept QUADRIPARTITION_SCORE_ALGORITHM = requires(T t, typename T::Color & color, typename T::ThreadPool & threadPool, common::AnnotatedBinaryTree & tree, std::size_t index, typename T::score_t score, void annotator(common::AnnotatedBinaryTree::Node*, typename T::score_t, typename T::score_t, typename T::score_t), int verbose) {
		requires std::derived_from<T, common::LogInfo>;
		requires std::same_as<array<typename T::score_t, 3>, typename T::ThreadPool::score_t>;
		requires thread_pool::THREAD_POOL<typename T::ThreadPool>;
	{ T{ color, threadPool, tree, verbose } };
	{ t.labelTree(size_t()) };
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


	// ---- EVAL FIX (M1) ---------------------------------------------------------
	// 旧实现把每个 (切割, 节点, bootstrap) 的三元组分数以 "PLACE44<cut><boot>" 为键
	// 塞进树节点的 std::map<std::string, array<score_t,3>>，条目数 = 切割数 × 节点数 ×
	// N_BOOTSTRAPS（n=201 时约 1.6e8 条 ≈ 21 GB），且永不释放。
	// 新实现：每个节点一份定长 vector<array<score_t,3>>，按 bootstrap 下标直接寻址；
	// 一个切割着色完成后立刻消费它的全部 bootstrap，下一个切割覆盖同一份内存。
	static inline std::string const BOOTSTRAP_VALUES_ATTR = "spr44_bootstrap_values";
	// 必须与 SPR44BlockBootstrap<Color>::MAIN_KEY 一致（见 quadripartition_support_spr44.hpp）
	static inline std::string const BOOTSTRAP_SUPPORT_KEY = "SPR44BlockBootstrap";

	template<typename score_t> using BootstrapValues = vector<std::array<score_t, 3> >;

	template<typename score_t>
	BootstrapValues<score_t>& bootstrapValues(common::AnnotatedBinaryTree::Node* node) {
		using BV = BootstrapValues<score_t>;
		if (!node->has<std::shared_ptr<BV> >(BOOTSTRAP_VALUES_ATTR)) {
			node->set(BOOTSTRAP_VALUES_ATTR, std::make_shared<BV>());
		}
		return *node->get<std::shared_ptr<BV> >(BOOTSTRAP_VALUES_ATTR);
	}

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

		// ---- EVAL FIX (M1): 注释路径的稠密状态（按节点编号索引，不写节点属性）----
		size_t nBootstraps = 0;
		size_t wordsPerNode = 0;
		vector<Tree::Node*> allNodes;
		std::unordered_map<Tree::Node*, size_t> nodeDenseId;
		vector<char> sprPrune, sprEnd, sprNum;
		vector<unsigned long long> markedBits; // allNodes.size() * wordsPerNode

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
			if (node->get<bool>("spr_prune")) return;
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
			if (node->get<bool>("spr_prune")) return;
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
			if constexpr (SPR44) node->set(PLACE44, Scores); // For SPR44 only
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
		vector<pair<typename Color::score_t, Tree::Node*> > ans = {};

		template<bool SPR44> vector<pair<typename Color::score_t, Tree::Node*> > optimalPlacementSubtree(Tree::Node* node, int flag, string key1, string key2) noexcept {
			using placeans_map = map<string, pair<typename Color::score_t, Tree::Node*> >;
			using place44_map = map<string, std::array<typename Color::score_t, 3> >;
			if constexpr (SPR44) {
				key1 = PLACE44;
				key2 = PlaceAns;
			}
			if (!node->has("placeans_of_node")) {
				shared_ptr<placeans_map> placeans_map_ptr(new placeans_map{});
				node->set("placeans_of_node", placeans_map_ptr);
			}

			if (node->isRoot()) {
				if constexpr (SPR44) {
					node->set(key2, make_pair(Color::ZERO, node));
					ans.push_back(node->get<pair<typename Color::score_t, Tree::Node*>>(key2));
				}
				else {
					placeans_map& map_p = *(node->get<shared_ptr<placeans_map> >("placeans_of_node"));
					map_p[key2] = { Color::ZERO, node };
					ans.push_back({ Color::ZERO, node });
				}
				if (!node->isLeaf() && !node->get<bool>("spr_prune")) {
					optimalPlacementSubtree<SPR44>(node->leftChild(), 1, key1, key2);
					optimalPlacementSubtree<SPR44>(node->rightChild(), 0, key1, key2);
				}
				return ans;
			}

			if (flag == 1) {
				if constexpr (SPR44) {
					typename Color::score_t scores = (node->parent())->get<pair<typename Color::score_t, Tree::Node*> >(key2).first - node->parent()->get<array<typename Color::score_t, 3> >(key1)[2] + node->parent()->get<array<typename Color::score_t, 3> >(key1)[0];
					pair<typename Color::score_t, Tree::Node*> ans_pair = { scores, node };
					node->set(key2, ans_pair);
					ans.push_back(ans_pair);
				}
				else {
					placeans_map& p_map_placeans = *(node->parent()->get<shared_ptr<placeans_map> >("placeans_of_node"));
					place44_map& p_map_place44 = *(node->parent()->get<shared_ptr<place44_map> >("place44_of_node"));
					typename Color::score_t scores = p_map_placeans[key2].first - p_map_place44[key1][2] + p_map_place44[key1][0];
					placeans_map& map_p = *(node->get<shared_ptr<placeans_map> >("placeans_of_node"));
					map_p[key2] = { scores, node };
					ans.push_back({ scores, node });
				}
			}
			else if (flag == 0) {
				if constexpr (SPR44) {
					typename Color::score_t scores = (node->parent())->get<pair<typename Color::score_t, Tree::Node*> >(key2).first - node->parent()->get<array<typename Color::score_t, 3> >(key1)[2] + node->parent()->get<array<typename Color::score_t, 3> >(key1)[1];
					pair<typename Color::score_t, Tree::Node*> ans_pair = { scores, node };
					node->set(key2, ans_pair);
					ans.push_back(ans_pair);
				}
				else {
					placeans_map& p_map_placeans = *(node->parent()->get<shared_ptr<placeans_map> >("placeans_of_node"));
					place44_map& p_map_place44 = *(node->parent()->get<shared_ptr<place44_map> >("place44_of_node"));
					typename Color::score_t scores = p_map_placeans[key2].first - p_map_place44[key1][2] + p_map_place44[key1][1];
					placeans_map& map_p = *(node->get<shared_ptr<placeans_map> >("placeans_of_node"));
					map_p[key2] = { scores, node };
					ans.push_back({ scores, node });
				}
			}

			if (!node->isLeaf() && !node->get<bool>("spr_prune")) {
				optimalPlacementSubtree<SPR44>(node->leftChild(), 1, key1, key2);
				optimalPlacementSubtree<SPR44>(node->rightChild(), 0, key1, key2);
			}
			return ans;
		}

		//遍历所有放置方式，找这个切割下的最好 
		template<bool SPR44> tuple<score_t, Tree::Node*> optimalPlacement(Tree& tree, Tree::Node* treeSlice, void annotator(Tree::Node*, score_t, score_t, score_t, string)) noexcept {
			vector<Tree::Node*> leafId_Slice = {};
			treeSlice->subtreeLeaves(leafId_Slice);

			int flag = node_to_num[treeSlice];
			string key1 = PLACE44 + std::to_string(flag);
			string key2 = PlaceAns + std::to_string(flag);

			colorTree<SPR44>(annotator, leafId_Slice, key1);

			if constexpr (SPR44) {
				ans.clear();
				optimalPlacementSubtree<SPR44>(tree.root(), -1, key1, key2);

				score_t standard_ans = 0;
				for (int i = 0; i < ans.size(); i++) {
					Tree::Node* node_i = ans[i].second;
					if (node_i->get<bool>("spr_prune")) {
						standard_ans = ans[i].first;
						break;
					}
				}

				auto it = std::max_element(
					ans.begin(),
					ans.end(),
					[](auto const& a, auto const& b) {
						return a.first < b.first;
					}
				);
				(*it).first -= standard_ans;
				return *it;
			}
			return {};
		}

	protected:
		vector<Tree::Node*> node_euler;
		vector<int> depth;
		map<Tree::Node*, int> node_to_num;
		map<int, Tree::Node*> num_to_node;
		vector<vector<int>> lca;
		vector<int> lg2;
		//int lca[410][410];
		//int lg2[410] = { 0 };
		void dfs(Tree::Node* node, int dep, bool cleared) {
			if (!cleared) {
				node_euler.clear();
				node_to_num.clear();
				num_to_node.clear();
				depth.clear();
				allNodes.clear();
				nodeDenseId.clear();
			}
			if (!nodeDenseId.contains(node)) {
				nodeDenseId[node] = allNodes.size();
				allNodes.push_back(node);
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
			lg2.push_back(-1);
			//lg2[0] = -1;
			for (int i = 1; i <= m; ++i) {
				lg2.push_back(lg2[i >> 1] + 1);
				//lg2[i] = lg2[i >> 1] + 1;
			}
			for (int i = 0; i <= m - 1; ++i) {
				lca.push_back({ i });
				//lca[i][0] = i;
			}
			for (int j = 1; j <= lg2[m]; ++j) {
				for (int i = 0; i + (1 << j) - 1 <= (m - 1); ++i) {
					lca[i].push_back(depth[lca[i][j - 1]] < depth[lca[i + (1 << (j - 1))][j - 1]] ? lca[i][j - 1] : lca[i + (1 << j - 1)][j - 1]);
					//lca[i][j] = depth[lca[i][j - 1]] < depth[lca[i + (1 << (j - 1))][j - 1]] ? lca[i][j - 1] : lca[i + (1 << j - 1)][j - 1];
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
			node->set("spr_prune", false);
			if (node->isLeaf()) return;
			setPrune(node->leftChild());
			setPrune(node->rightChild());
			return;
		}

		//遍历所有切割方式，得到最好的切割和组合
		vector<tuple<score_t, Tree::Node*, Tree::Node*, bool>> final_ans;
		template<bool SPR44> vector<tuple<score_t, Tree::Node*, Tree::Node*, bool>> regraft(Tree& tree, Tree::Node* tree_change, void annotator(Tree::Node*, score_t, score_t, score_t, string)) {
			setPrune(tree.root());
			tree_change->set("spr_prune", true);

			tuple<score_t, Tree::Node*> ans1 = optimalPlacement<SPR44>(tree, tree_change, annotator);

			if constexpr (SPR44) {
				Tree::Node* regraft_site = get<1>(ans1);
				bool flagprune = 0;
				if (!(regraft_site->get<bool>("spr_prune"))) {
					if (regraft_site->isLeaf()) {
						flagprune = true;
					}
					else if (!(regraft_site->leftChild()->get<bool>("spr_prune") || regraft_site->rightChild()->get<bool>("spr_prune"))) {
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
		// ---- EVAL FIX (M1): 注释路径的支撑度计算 -------------------------------
		// 旧的实现：single_iBootstrap_support 用 "PLACE44<cut><boot>" / "PlaceAns<cut><boot>"
		// 字符串键在节点属性 map 里存取；final_support 在每个 bootstrap 上把所有切割重走一遍。
		// 新实现：切割外层、bootstrap 内层（同一份着色结果被复用），支撑度用
		// "每节点 × 每 bootstrap 一位" 的 bitset 累积，最后一次性算百分比。
		void initBootstrapState(size_t nBootstrap) {
			nBootstraps = nBootstrap;
			wordsPerNode = (nBootstraps + 63) / 64;
			size_t const n = allNodes.size();
			sprPrune.assign(n, 0);
			sprEnd.assign(n, 0);
			sprNum.assign(n, 0);
			markedBits.assign(n * wordsPerNode, 0ULL);
			for (Tree::Node* node : allNodes) {
				bootstrapValues<typename Color::score_t>(node).assign(nBootstraps, std::array<typename Color::score_t, 3>{ Color::ZERO, Color::ZERO, Color::ZERO });
			}
			RMQ();
		}

		void syncPruneDense(Tree::Node* node) {
			sprPrune[nodeDenseId[node]] = node->get<bool>("spr_prune") ? 1 : 0;
			if (node->isLeaf()) return;
			syncPruneDense(node->leftChild());
			syncPruneDense(node->rightChild());
		}

		void get_num(Tree::Node* node) {
			while (!sprEnd[nodeDenseId[node]] && !node->isRoot()) {
				sprNum[nodeDenseId[node]] = 1;
				node = node->parent();
			}
		}

		// 放置得分只依赖父节点的得分与父节点的 (s0,s1,s2)，所以沿递归传给子节点即可，无需存储。
		void collectPlaceans(Tree::Node* node, typename Color::score_t parentScore, array<typename Color::score_t, 3> const& parentVals, int flag, size_t iBootstrap) {
			typename Color::score_t dp = Color::ZERO;
			if (!node->isRoot()) dp = parentScore - parentVals[2] + parentVals[(flag == 1) ? 0 : 1];
			ans.push_back({ dp, node });
			if (!node->isLeaf() && !sprPrune[nodeDenseId[node]]) {
				array<typename Color::score_t, 3> const& myVals = bootstrapValues<typename Color::score_t>(node)[iBootstrap];
				collectPlaceans(node->leftChild(), dp, myVals, 1, iBootstrap);
				collectPlaceans(node->rightChild(), dp, myVals, 0, iBootstrap);
			}
		}

		// 单个 (切割, bootstrap)：求每节点放置得分，按原逻辑标记，再把标记 OR 进 bitset。
		void markCutForBootstrap(Tree::Node* treeSlice, size_t iBootstrap) {
			std::fill(sprNum.begin(), sprNum.end(), (char)0);
			std::fill(sprEnd.begin(), sprEnd.end(), (char)0);
			ans.clear();
			array<typename Color::score_t, 3> const zero = { Color::ZERO, Color::ZERO, Color::ZERO };
			collectPlaceans(tree.root(), Color::ZERO, zero, -1, iBootstrap);
			typename Color::score_t standard = Color::ZERO;
			for (auto const& [s, nd] : ans) {
				if (sprPrune[nodeDenseId[nd]]) { standard = s; break; }
			}
			for (auto const& [s, nd] : ans) {
				if (sprPrune[nodeDenseId[nd]]) continue;
				if (s >= standard + Color::EPSILON) {
					Tree::Node* const lca = getlca(treeSlice, nd);
					sprEnd[nodeDenseId[lca]] = 1;
					if (!treeSlice->isRoot()) get_num(treeSlice->parent());
					if (!sprEnd[nodeDenseId[treeSlice->parent()]]) sprNum[nodeDenseId[treeSlice]] = 1;
					if (!nd->isRoot()) get_num(nd);
					sprEnd[nodeDenseId[lca]] = 0;
				}
			}
			unsigned long long const mask = 1ULL << (iBootstrap % 64);
			size_t const word = iBootstrap / 64;
			for (size_t i = 0; i < allNodes.size(); ++i) {
				if (sprNum[i]) markedBits[i * wordsPerNode + word] |= mask;
			}
		}

		// 一个切割：着色一次 -> 立刻把该切割的所有 bootstrap 消费掉。
		void annotateCuts(Tree::Node* treeSlice, void annotator(Tree::Node*, score_t, score_t, score_t, string)) {
			setPrune(tree.root());
			treeSlice->set("spr_prune", true);
			syncPruneDense(tree.root());
			{
				vector<Tree::Node*> sliceLeaves;
				treeSlice->subtreeLeaves(sliceLeaves);
				colorTree<false>(annotator, sliceLeaves, string());
			}
			for (size_t iBootstrap = 0; iBootstrap < nBootstraps; ++iBootstrap) {
				markCutForBootstrap(treeSlice, iBootstrap);
			}
			treeSlice->set("spr_prune", false);
			if (!treeSlice->isLeaf()) {
				annotateCuts(treeSlice->leftChild(), annotator);
				annotateCuts(treeSlice->rightChild(), annotator);
			}
		}

		void finishBootstrapSupport(size_t nBootstrap) {
			for (size_t i = 0; i < allNodes.size(); ++i) {
				size_t cnt = 0;
				for (size_t b = 0; b < nBootstraps; ++b) {
					if (((markedBits[i * wordsPerNode + b / 64] >> (b % 64)) & 1ULL) == 0) cnt++;
				}
				double const value = (nBootstrap == 0) ? 0.0 : (double)cnt * 100.0 / (double)nBootstrap;
				allNodes[i]->set(BOOTSTRAP_SUPPORT_KEY, value);
			}
		}

	protected:
		//得到每个切割方式最好的重接方式后进行不同切割方式的比较
		template<bool SPR44> tuple<score_t, Tree::Node*, Tree::Node*, bool> optimalSpr44(void annotator(Tree::Node*, score_t, score_t, score_t, string), size_t N_BOOTSTRAPS) noexcept {
			dfs(tree.root(), 0, 0);
			if constexpr (SPR44) {
				final_ans.clear();
				final_ans = regraft<SPR44>(tree, tree.root()->leftChild(), annotator);
				final_ans = regraft<SPR44>(tree, tree.root()->rightChild(), annotator);
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
				// EVAL FIX (M1): 切割外层 + bootstrap 内层，内存与数据集大小无关
				initBootstrapState(N_BOOTSTRAPS);
				annotateCuts(tree.root()->leftChild(), annotator);
				annotateCuts(tree.root()->rightChild(), annotator);
				finishBootstrapSupport(N_BOOTSTRAPS);
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
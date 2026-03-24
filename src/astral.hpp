#ifndef ASTRAL_HPP
#define ASTRAL_HPP

#include "stepwise_colorable.hpp"

namespace astral {

using std::size_t;
using std::views::iota;
using std::vector;
using std::array;
using std::string;
using std::unique_ptr;
using std::shared_ptr;

template<class Attributes> concept STEPWISE_COLOR_ATTRIBUTES = requires
{
	requires std::integral<typename Attributes::score_t> || std::floating_point<typename Attributes::score_t>;
	{ Attributes::ZERO } -> std::convertible_to<typename Attributes::score_t>;
	{ Attributes::EPSILON } -> std::convertible_to<typename Attributes::score_t>;
};

template<typename score_type, bool wSupport, bool wLength> struct StepwiseColorDefaultAttributes {
	using index_t = long long;
	using score_t = score_type;
	using support_t = std::conditional_t<wSupport, score_t, bool>;
	using length_t = std::conditional_t<wLength, score_t, bool>;

	static inline bool constexpr WEIGHTING_BY_SUPPORT = wSupport;
	static inline bool constexpr WEIGHTING_BY_LENGTH = wLength;
	static inline score_t constexpr ZERO = 0;
	static inline score_t constexpr EPSILON = [] { if constexpr (std::integral<score_t>) return 0; else return 1e-3; }();
};

template<STEPWISE_COLOR_ATTRIBUTES Attributes> class Color {
public:
	using score_t = Attributes::score_t;
	using support_t = Attributes::support_t;
	using length_t = Attributes::length_t;
	using index_t = Attributes::index_t;

	static inline bool constexpr WEIGHTING_BY_SUPPORT = Attributes::WEIGHTING_BY_SUPPORT;
	static inline bool constexpr WEIGHTING_BY_LENGTH = Attributes::WEIGHTING_BY_LENGTH;
	static inline bool constexpr IS_ROOTED = false;
	static inline score_t constexpr ZERO = Attributes::ZERO;
	static inline score_t constexpr EPSILON = Attributes::EPSILON;

	struct SharedConstData {
		using ParentClass = Color<Attributes>;

		struct Element {
			struct Node {
				unique_ptr<Node> lc, rc;
				index_t leafID = -1;
				support_t wSupport = 0;
				length_t wLength = 1;

				bool isLeaf() const noexcept { return !lc; }
			};

			shared_ptr<Node> root;
		};

		vector<Element> elements;
		size_t nTaxa = 0;
		
		size_t nElements() const noexcept { return elements.size(); }
	};

	struct Tree {
		struct Node {
			index_t parent = -1, lc = -1, rc = -1;
			support_t wSupport = 1;
			length_t wLength = 0;
			score_t score = 0;

			score_t x = 0, y = 0, z = 0, tx = 0, ty = 0, tz = 0;
			score_t x2a = 0, y2a = 0, z2a = 0, xya = 0, xza = 0, yza = 0;
			score_t x2b = 0, y2b = 0, z2b = 0, xyb = 0, xzb = 0, yzb = 0;
		};

		struct QNode {
			score_t w = 0, xy_zw = 0, xz_yw = 0, xw_yz = 0;

			score_t xwa = 0, ywa = 0, zwa = 0;
			score_t xwb = 0, ywb = 0, zwb = 0;

			score_t y_zw = 0, x_zw = 0, xy_w = 0, xy_z = 0;
			score_t z_yw = 0, x_yw = 0, xz_w = 0, xz_y = 0;
			score_t w_yz = 0, x_yz = 0, xw_z = 0, xw_y = 0;
		};

		vector<Node> nodes;
		vector<QNode> qnodes;
		vector<index_t> taxonID2nodeID; // taxonID2nodeID[iTaxon] -> iNode in nodes

		index_t build(typename SharedConstData::Element::Node const& sNode) noexcept {
			index_t lc = -1, rc = -1;
			if (!sNode.isLeaf()) {
				lc = build(*sNode.lc.get());
				rc = build(*sNode.rc.get());
			}
			index_t iNode = nodes.size();
			nodes.emplace_back();
			qnodes.emplace_back();
			Node& node = nodes.back();
			node.wSupport = sNode.wSupport;
			node.wLength = sNode.wLength;
			node.lc = lc;
			node.rc = rc;
			if (lc == -1) {
				taxonID2nodeID[sNode.leafID] = iNode;
			}
			else {
				nodes[lc].parent = iNode;
				nodes[rc].parent = iNode;
			}
			return iNode;
		}

		Tree(typename SharedConstData::Element const& element, size_t nTaxa) noexcept : taxonID2nodeID(nTaxa, -1) {
			build(*element.root);
		}

		void update(index_t iNode, bool quadMode) noexcept {
			if (iNode == -1) return;

			Node& w = nodes[iNode];
			Node& u = nodes[w.lc];
			Node& v = nodes[w.rc];

			length_t wLength = w.wLength;
			support_t wSupport = w.wSupport;

			w.x = (u.x + v.x) * wLength;
			w.y = (u.y + v.y) * wLength;
			w.z = (u.z + v.z) * wLength;

			w.x2a = u.x2a + v.x2a + u.x * v.x;
			w.y2a = u.y2a + v.y2a + u.y * v.y;
			w.z2a = u.z2a + v.z2a + u.z * v.z;
			w.xya = u.xya + v.xya + u.x * v.y + u.y * v.x;
			w.xza = u.xza + v.xza + u.x * v.z + u.z * v.x;
			w.yza = u.yza + v.yza + u.y * v.z + u.z * v.y;

			w.x2b = (u.x2b + v.x2b + u.x * v.x) * wSupport;
			w.y2b = (u.y2b + v.y2b + u.y * v.y) * wSupport;
			w.z2b = (u.z2b + v.z2b + u.z * v.z) * wSupport;
			w.xyb = (u.xyb + v.xyb + u.x * v.y + u.y * v.x) * wSupport;
			w.xzb = (u.xzb + v.xzb + u.x * v.z + u.z * v.x) * wSupport;
			w.yzb = (u.yzb + v.yzb + u.y * v.z + u.z * v.y) * wSupport;

			w.tx = (u.tx + v.tx + u.y * (v.z2a - v.z2b) + (u.z2a - u.z2b) * v.y + u.z * (v.y2a - v.y2b) + (u.y2a - u.y2b) * v.z) * wLength;
			w.ty = (u.ty + v.ty + u.x * (v.z2a - v.z2b) + (u.z2a - u.z2b) * v.x + u.z * (v.x2a - v.x2b) + (u.x2a - u.x2b) * v.z) * wLength;
			w.tz = (u.tz + v.tz + u.x * (v.y2a - v.y2b) + (u.y2a - u.y2b) * v.x + u.y * (v.x2a - v.x2b) + (u.x2a - u.x2b) * v.y) * wLength;

			w.score = u.x * v.tx + u.y * v.ty + u.z * v.tz + v.x * u.tx + v.y * u.ty + v.z * u.tz
				+ u.x2a * v.yza - u.x2b * v.yzb
				+ u.y2a * v.xza - u.y2b * v.xzb
				+ u.z2a * v.xya - u.z2b * v.xyb
				+ u.score + v.score;

			if (quadMode) {
				QNode& W = qnodes[iNode];
				QNode& U = qnodes[w.lc];
				QNode& V = qnodes[w.rc];

				W.w = (U.w + V.w) * wLength;

				W.xwa = U.xwa + V.xwa + u.x * V.w + U.w * v.x;
				W.ywa = U.ywa + V.ywa + u.y * V.w + U.w * v.y;
				W.zwa = U.zwa + V.zwa + u.z * V.w + U.w * v.z;

				W.xwb = (U.xwb + V.xwb + u.x * V.w + U.w * v.x) * wSupport;
				W.ywb = (U.ywb + V.ywb + u.y * V.w + U.w * v.y) * wSupport;
				W.zwb = (U.zwb + V.zwb + u.z * V.w + U.w * v.z) * wSupport;

				W.y_zw = (U.y_zw + V.y_zw + u.y * (V.zwa - V.zwb) + (U.zwa - U.zwb) * v.y) * wLength;
				W.x_zw = (U.x_zw + V.x_zw + u.x * (V.zwa - V.zwb) + (U.zwa - U.zwb) * v.x) * wLength;
				W.xy_w = (U.xy_w + V.xy_w + U.w * (v.xya - v.xyb) + (u.xya - u.xyb) * V.w) * wLength;
				W.xy_z = (U.xy_z + V.xy_z + u.z * (v.xya - v.xyb) + (u.xya - u.xyb) * v.z) * wLength;
				W.xy_zw = U.xy_zw + V.xy_zw + u.xya * V.zwa - u.xyb * V.zwb + U.zwa * v.xya - U.zwb * v.xyb
					+ u.x * V.y_zw + U.y_zw * v.x
					+ u.y * V.x_zw + U.x_zw * v.y
					+ u.z * V.xy_w + U.xy_w * v.z
					+ U.w * V.xy_z + U.xy_z * V.w;

				W.z_yw = (U.z_yw + V.z_yw + u.z * (V.ywa - V.ywb) + (U.ywa - U.ywb) * v.z) * wLength;
				W.x_yw = (U.x_yw + V.x_yw + u.x * (V.ywa - V.ywb) + (U.ywa - U.ywb) * v.x) * wLength;
				W.xz_w = (U.xz_w + V.xz_w + U.w * (v.xza - v.xzb) + (u.xza - u.xzb) * V.w) * wLength;
				W.xz_y = (U.xz_y + V.xz_y + u.y * (v.xza - v.xzb) + (u.xza - u.xzb) * v.y) * wLength;
				W.xz_yw = U.xz_yw + V.xz_yw + u.xza * V.ywa - u.xzb * V.ywb + U.ywa * v.xza - U.ywb * v.xzb
					+ u.x * V.z_yw + U.z_yw * v.x
					+ u.z * V.x_yw + U.x_yw * v.z
					+ u.y * V.xz_w + U.xz_w * v.y
					+ U.w * V.xz_y + U.xz_y * V.w;

				W.w_yz = (U.w_yz + V.w_yz + U.w * (v.yza - v.yzb) + (u.yza - u.yzb) * V.w) * wLength;
				W.x_yz = (U.x_yz + V.x_yz + u.x * (v.yza - v.yzb) + (u.yza - u.yzb) * v.x) * wLength;
				W.xw_z = (U.xw_z + V.xw_z + u.z * (V.xwa - V.xwb) + (U.xwa - U.xwb) * v.z) * wLength;
				W.xw_y = (U.xw_y + V.xw_y + u.y * (V.xwa - V.xwb) + (U.xwa - U.xwb) * v.y) * wLength;
				W.xw_yz = U.xw_yz + V.xw_yz + U.xwa * v.yza - U.xwb * v.yzb + u.yza * V.xwa - u.yzb * V.xwb
					+ u.x * V.w_yz + U.w_yz * v.x
					+ U.w * V.x_yz + U.x_yz * V.w
					+ u.y * V.xw_z + U.xw_z * v.y
					+ u.z * V.xw_y + U.xw_y * v.z;
			}

			update(w.parent, quadMode);
		}
	};

	void recursiveUpdate(size_t iElement, index_t iNode) noexcept {
		trees[iElement].update(iNode, quadMode);
		scores[iElement] = trees[iElement].nodes.back().score;
	}

	vector<Tree> trees;
	vector<score_t> scores;
	bool quadMode = false;

public:
	Color(SharedConstData const& data) noexcept : scores(data.nElements()) {
		for (auto const& element : data.elements) {
			trees.emplace_back(element, data.nTaxa);
		}
	}

	void elementSetTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept {
		Tree& tree = trees[iElement];
		index_t iNode = tree.taxonID2nodeID[iTaxon];
		typename Tree::Node& node = tree.nodes[iNode];
		// Assuming x2a, xya or alike are all zeros
		score_t wLength = node.wLength;
		if (iColor == 0) node.x += wLength;
		if (iColor == 1) node.y += wLength;
		if (iColor == 2) node.z += wLength;
		if (iColor == 3) tree.qnodes[iNode].w += wLength;
		recursiveUpdate(iElement, node.parent);
	}

	void elementClearTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept {
		Tree& tree = trees[iElement];
		index_t iNode = tree.taxonID2nodeID[iTaxon];
		typename Tree::Node& node = tree.nodes[iNode];
		// Assuming x2a, xya or alike are all zeros
		score_t wLength = node.wLength;
		if (iColor == 0) node.x -= wLength;
		if (iColor == 1) node.y -= wLength;
		if (iColor == 2) node.z -= wLength;
		if (iColor == 3) tree.qnodes[iNode].w -= wLength;
		recursiveUpdate(iElement, node.parent);
	}

	void elementClearAndSetTaxonColor(size_t iElement, size_t iTaxon, size_t iColor, size_t jColor) noexcept {
		Tree& tree = trees[iElement];
		index_t iNode = tree.taxonID2nodeID[iTaxon];
		typename Tree::Node& node = tree.nodes[iNode];
		// Assuming x2a, xya or alike are all zeros
		score_t wLength = node.wLength;
		if (iColor == 0) node.x -= wLength;
		if (iColor == 1) node.y -= wLength;
		if (iColor == 2) node.z -= wLength;
		if (iColor == 3) tree.qnodes[iNode].w -= wLength;
		if (jColor == 0) node.x += wLength;
		if (jColor == 1) node.y += wLength;
		if (jColor == 2) node.z += wLength;
		if (jColor == 3) tree.qnodes[iNode].w += wLength;
		recursiveUpdate(iElement, node.parent);
	}

	score_t elementScore(size_t iElement) const noexcept { return scores[iElement]; }

	array<score_t, 3> elementQuadripartitionScores(size_t iElement) const noexcept {
		typename Tree::QNode const& qnode = trees[iElement].qnodes.back();
		return { qnode.xy_zw, qnode.xz_yw, qnode.xw_yz };
	}

	void setQuadripartitionMode(bool mode) noexcept { quadMode = mode; }
};

namespace DriverHelper {

	using namespace std;

	namespace Newick{
		void skipWhitespace(string_view &sv) noexcept { while (!sv.empty() && isspace(sv.front())) sv.remove_prefix(1); }

		pair<string, string> splitLabelAndLength(string_view& sv) noexcept {
			skipWhitespace(sv);
			size_t i = 0;
			while (i < sv.size() && sv[i] != ':' && sv[i] != ',' && sv[i] != ')' && sv[i] != ';') i++;
			string label = string(sv.substr(0, i));
			sv.remove_prefix(i);
			string length;
			if (!sv.empty() && sv[0] == ':') {
				sv.remove_prefix(1);
				skipWhitespace(sv);
				i = 0;
				while (i < sv.size() && sv[i] != ',' && sv[i] != ')' && sv[i] != ';') i++;
				length = string(sv.substr(0, i));
				sv.remove_prefix(i);
			}
			return { label, length };
		}

		optional<double> parseSupport(string const& s) noexcept {
			try {
				return std::stod(s);
			}
			catch (...) {
				return std::nullopt;
			}
		}

		optional<double> parseLength(string const& s) noexcept {
			try {
				return std::stod(s);
			}
			catch (...) {
				return std::nullopt;
			}
		}
	}

	struct Node {
		vector<shared_ptr<Node> > children;
		string label, length;
		bool isReal = true;
	
		Node(string_view& sv) { // Assuming the input is a valid Newick format
			Newick::skipWhitespace(sv);
			if (sv[0] == '(') {
				sv.remove_prefix(1);
				while (true) {
					children.emplace_back(make_shared<Node>(sv));
					Newick::skipWhitespace(sv);
					if (sv[0] == ',') {
						sv.remove_prefix(1);
						continue;
					}
					else if (sv[0] == ')') {
						sv.remove_prefix(1);
						break;
					}
					else throw std::invalid_argument("Invalid Newick format");
				}
			}
			tie(label, length) = Newick::splitLabelAndLength(sv);
		}

		Node(shared_ptr<Node> lc, shared_ptr<Node> rc, bool isReal): isReal(isReal){
			children.push_back(lc);
			children.push_back(rc);
		}

		int huffman() noexcept {
			if (children.empty()) {
				common::taxonName2ID[label];
				return 1;
			}
			auto cmp = [](pair<shared_ptr<Node>,int> const &left, pair<shared_ptr<Node>,int> const &right) { return left.second < right.second; };
			priority_queue<pair<shared_ptr<Node>, int>, vector<pair<shared_ptr<Node>, int> >, decltype(cmp)> pq(cmp);
			for (auto& child : children) pq.push({child, child->huffman()});
			children.clear();
			while (pq.size() > 2) {
				auto [lc, hlc] = pq.top(); pq.pop();
				auto [rc, hrc] = pq.top(); pq.pop();
				pq.push({ make_shared<Node>(lc, rc, false), hlc + hrc });
			}
			{
				auto [lc, hlc] = pq.top(); pq.pop();
				auto [rc, hrc] = pq.top(); pq.pop();
				children.push_back(lc);
				children.push_back(rc);
				return hlc + hrc;
			}
		}

		void minMaxSupport(double &minSupport, double &maxSupport) const noexcept {
			if (children.empty()) return;
			for (auto& child : children) child->minMaxSupport(minSupport, maxSupport);
			if (isReal) {
				optional<double> support = Newick::parseSupport(label);
				if (support){
					minSupport = std::min(minSupport, *support);
					maxSupport = std::max(maxSupport, *support);
				}
			}
		}

		template<typename Color> Color::SharedConstData::Element::Node* convert(double minSupport, double maxSupport) noexcept {
			using DataNode = typename Color::SharedConstData::Element::Node;
			DataNode* node = new DataNode();

			if (!children.empty()) {
				node->lc.reset(children[0]->convert<Color>(minSupport, maxSupport));
				node->rc.reset(children[1]->convert<Color>(minSupport, maxSupport));
				if constexpr (Color::WEIGHTING_BY_SUPPORT) {
					optional<double> support = Newick::parseSupport(label);
					if (support) node->wSupport = 1 - (*support - minSupport) / (maxSupport - minSupport);
				}
			}
			else {
				node->leafID = common::taxonName2ID[label];
			}

			if constexpr (Color::WEIGHTING_BY_LENGTH) {
				optional<double> len = Newick::parseLength(length);
				if (len) node->wLength = exp(- *len);
			}

			return node;
		}
	};

	template<typename DataClass> DataClass read() {
		ifstream fin(ARG.get<string>("input"));
		string line;
		vector<shared_ptr<Node> > trees;
		double minSupport = 0.333, maxSupport = 1;
		while (getline(fin, line)) {
			string_view sv(line);
			shared_ptr<Node> root(new Node(sv));
			root->huffman();
			if constexpr (DataClass::ParentClass::WEIGHTING_BY_SUPPORT) root->minMaxSupport(minSupport, maxSupport);
			trees.push_back(root);
		}
		if (minSupport < 0.33) minSupport = 0;
		if (maxSupport > 2) {
			minSupport = 0;
			maxSupport = 100;
		}
		common::LogInfo log(1);
		if constexpr (DataClass::ParentClass::WEIGHTING_BY_SUPPORT){
			log.log() << "Min support value = " << minSupport << std::endl;
			log.log() << "Max support value = " << maxSupport << std::endl;
		}

		DataClass data;
		data.nTaxa = common::taxonName2ID.nTaxa();
		for (auto& tree : trees) {
			data.elements.emplace_back();
			data.elements.back().root.reset(tree->convert<typename DataClass::ParentClass>(minSupport, maxSupport));
		}
		return data;
	}

};

template<bool> class Driver : public common::LogInfo
{
	using string = std::string;

public:
	using DataClasses = std::variant<typename Color<StepwiseColorDefaultAttributes<double, true, true> >::SharedConstData, typename Color<StepwiseColorDefaultAttributes<double, true, false> >::SharedConstData, typename Color<StepwiseColorDefaultAttributes<double, false, true> >::SharedConstData, typename Color<StepwiseColorDefaultAttributes<unsigned long long, false, false> >::SharedConstData>;

	static std::pair<string, string> programNames() noexcept {
		return { "astral", "Accurate Species Tree ALgorithm (ASTRAL-IV)" };
	}

	static void addArguments() noexcept {
		ARG.addArgument('\0', "mode", "integer", "Weighting mode (1: hybrid weighting, 2: support only, 3: length only, 4: unweighted)", 0, true, true, "1");
	}

	static DataClasses getStepwiseColorSharedConstData() noexcept {
		size_t mode = ARG.get<size_t>("mode");
		if (mode == 1) return DriverHelper::read<typename Color<StepwiseColorDefaultAttributes<double, true, true> >::SharedConstData>();
		if (mode == 2) return DriverHelper::read<typename Color<StepwiseColorDefaultAttributes<double, true, false> >::SharedConstData>();
		if (mode == 3) return DriverHelper::read<typename Color<StepwiseColorDefaultAttributes<double, false, true> >::SharedConstData>();
		if (mode == 4) return DriverHelper::read<typename Color<StepwiseColorDefaultAttributes<unsigned long long, false, false> >::SharedConstData>();
		
		common::LogInfo vlog(-99);
		vlog.log() << "Error: Invalid mode: " << mode << std::endl;
		exit(-1);
	}
};

class Documentation : public common::DocumentationBase {
protected:
	string introduction() const noexcept override {
		return R"YOHANETYO(# Accurate Species Tree ALgorithm (ASTRAL-IV)
)YOHANETYO";
	}

	string input() const noexcept override {
		return R"YOHANETYO(# INPUT
)YOHANETYO";
	}

	string programName() const noexcept override { return "astral"; }

	string exampleInput() const noexcept override { return "genetrees.nw"; }
};

};
#endif // !ASTRAL_HPP

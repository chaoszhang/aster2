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
			support_t wSupport = 0;
			length_t wLength = 1;
			score_t score = 0;

			score_t x = 0, y = 0, z = 0, tx = 0, ty = 0, tz = 0, q = 0;
			score_t x2a = 0, y2a = 0, z2a = 0, xya = 0, xza = 0, yza = 0;
			score_t x2b = 0, y2b = 0, z2b = 0, xyb = 0, xzb = 0, yzb = 0;
		};

		vector<Node> nodes;
		vector<index_t> taxonID2nodeID; // taxonID2nodeID[iTaxon] -> iNode in nodes

		index_t build(typename SharedConstData::Element::Node const& sNode) noexcept {
			index_t lc = -1, rc = -1;
			if (!sNode.isLeaf()) {
				lc = build(*sNode.lc.get());
				rc = build(*sNode.rc.get());
			}
			index_t iNode = nodes.size();
			nodes.emplace_back();
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

		void update(index_t iNode) noexcept {
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

			w.tx = (u.tx + v.tx + u.y * (v.z2a - v.z2b) + (u.z2a - u.z2b) * v.y + u.z * (v.y2a - v.y2b) + (u.y2a - u.y2b) * v.z) * w.wLength;
			w.ty = (u.ty + v.ty + u.x * (v.z2a - v.z2b) + (u.z2a - u.z2b) * v.x + u.z * (v.x2a - v.x2b) + (u.x2a - u.x2b) * v.z) * w.wLength;
			w.tz = (u.tz + v.tz + u.x * (v.y2a - v.y2b) + (u.y2a - u.y2b) * v.x + u.y * (v.x2a - v.x2b) + (u.x2a - u.x2b) * v.y) * w.wLength;

			w.score = u.x * v.tx + u.y * v.ty + u.z * v.tz + v.x * u.tx + v.y * u.ty + v.z * u.tz
				+ u.x2a * v.yza - u.x2b * v.yzb
				+ u.y2a * v.xza - u.y2b * v.xzb
				+ u.z2a * v.xya - u.z2b * v.xyb
				+ u.score + v.score;

			update(w.parent);
		}
	};

	void recursiveUpdate(size_t iElement, index_t iNode) noexcept {
		trees[iElement].update(iNode);
		scores[iElement] = trees[iElement].nodes.back().score;
	}

	vector<Tree> trees;
	vector<score_t> scores;

public:
	Color(SharedConstData const& data) noexcept : scores(data.nElements()) {
		for (auto const& element : data.elements) {
			trees.emplace_back(element, data.nTaxa);
		}
	}

	void elementSetTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept {
		Tree& tree = trees[iElement];
		typename Tree::Node& node = tree.nodes[tree.taxonID2nodeID[iTaxon]];
		// Assuming x2a, xya or alike are all zeros
		if (iColor == 0) node.x++;
		if (iColor == 1) node.y++;
		if (iColor == 2) node.z++;
		recursiveUpdate(iElement, node.parent);
	}

	void elementClearTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept {
		Tree& tree = trees[iElement];
		typename Tree::Node& node = tree.nodes[tree.taxonID2nodeID[iTaxon]];
		// Assuming x2a, xya or alike are all zeros
		if (iColor == 0) node.x--;
		if (iColor == 1) node.y--;
		if (iColor == 2) node.z--;
		recursiveUpdate(iElement, node.parent);
	}

	void elementClearAndSetTaxonColor(size_t iElement, size_t iTaxon, size_t iColor, size_t jColor) noexcept {
		Tree& tree = trees[iElement];
		typename Tree::Node& node = tree.nodes[tree.taxonID2nodeID[iTaxon]];
		// Assuming x2a, xya or alike are all zeros
		if (iColor == 0) node.x--;
		if (iColor == 1) node.y--;
		if (iColor == 2) node.z--;
		if (jColor == 0) node.x++;
		if (jColor == 1) node.y++;
		if (jColor == 2) node.z++;
		recursiveUpdate(iElement, node.parent);
	}

	score_t elementScore(size_t iElement) const noexcept { return scores[iElement]; }
};

namespace DriverHelper {

using namespace std;

template<typename DataClass> DataClass read() {
	return DataClass();
}

};

template<bool> class Driver : public common::LogInfo
{
	using string = std::string;

public:
	//using DataClasses = std::variant<typename Color<StepwiseColorDefaultAttributes<double, true, true> >::SharedConstData, typename Color<StepwiseColorDefaultAttributes<double, false, true> >::SharedConstData, typename Color<StepwiseColorDefaultAttributes<double, true, false> >::SharedConstData, typename Color<StepwiseColorDefaultAttributes<unsigned long long, false, false> >::SharedConstData>;
	using DataClasses = std::variant<typename Color<StepwiseColorDefaultAttributes<double, true, true> >::SharedConstData>;

	static std::pair<string, string> programNames() noexcept {
		return { "astral", "Accurate Species Tree ALgorithm (ASTRAL-IV)" };
	}

	static void addArguments() noexcept {
		//ARG.addArgument('\0', "chunk", "integer", "The maximum number of sites in each local aligment block for parameter estimation", 0, true, true, "10000");
	}

	static DataClasses getStepwiseColorSharedConstData() noexcept {
		/*
		try { return DriverHelper::read<std::variant_alternative_t<0, DataClasses> >(); }
		catch (...) {}
		try { return DriverHelper::read<std::variant_alternative_t<1, DataClasses> >(); }
		catch (...) {}
		try { return DriverHelper::read<std::variant_alternative_t<2, DataClasses> >(); }
		catch (...) {}
		try { return DriverHelper::read<std::variant_alternative_t<3, DataClasses> >(); }
		catch (...) {}
		*/
		return DriverHelper::read<std::variant_alternative_t<0, DataClasses> >();
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

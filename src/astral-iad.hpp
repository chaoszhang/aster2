#ifndef ASTRAL_IAD_HPP
#define ASTRAL_IAD_HPP

#include "stepwise_colorable.hpp"
#include "quadripartition_support.hpp"
#include "quadripartition_support_spr44.hpp"

namespace astral_iad {

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

	template<typename score_type> struct StepwiseColorDefaultAttributes {
		using index_t = long long;
		using score_t = score_type;

		static inline score_t constexpr ZERO = 0;
		static inline score_t constexpr EPSILON = [] { if constexpr (std::integral<score_t>) return 0; else return 1e-6; }();
	};

	template<STEPWISE_COLOR_ATTRIBUTES Attributes> class Color {
	public:
		using score_t = Attributes::score_t;
		using index_t = Attributes::index_t;

		static inline bool constexpr IS_ROOTED = false;
		static inline score_t constexpr ZERO = Attributes::ZERO;
		static inline score_t constexpr EPSILON = Attributes::EPSILON;

		struct SharedConstData {
			using ParentClass = Color<Attributes>;

			struct Element {
				struct Node {
					unique_ptr<Node> lc, rc;
					index_t leafID = -1;
					score_t length = 0;

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
				index_t x = 0, y = 0, z = 0;
				score_t length = 0;
				score_t cox2 = 0, coy2 = 0, coz2 = 0, co1 = 0, co0 = 0;
			};

			struct QNode {
				// TODO
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
				node.length = sNode.length;
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
				
				w.x = u.x + v.x;
				w.y = u.y + v.y;
				w.z = u.z + v.z;
				
				score_t xyzl = w.x * w.y * w.z * w.length;

				w.cox2 = u.cox2 + v.cox2 + w.y * w.z * w.length;
				w.coy2 = u.coy2 + v.coy2 + w.x * w.z * w.length;
				w.coz2 = u.coz2 + v.coz2 + w.x * w.y * w.length;
				w.co1 = u.co1 + v.co1 - 2 * xyzl;
				w.co0 = u.co0 + v.co0 + (w.x + w.y + w.z + 3) * xyzl;

				if (quadMode) {
					// TODO
				}

				update(w.parent, quadMode);
			}
		};

		void recursiveUpdate(size_t iElement, index_t iNode) noexcept {
			trees[iElement].update(iNode, quadMode);
			typename Tree::Node& r = trees[iElement].nodes.back();
			scores[iElement] = (r.x * (r.x - 1) * r.cox2 + r.y * (r.y - 1) * r.coy2 + r.z * (r.z - 1) * r.coz2 + (r.x + r.y + r.z) * r.co1 + r.co0) / 2;
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
			if (iColor == 0) node.x++;
			if (iColor == 1) node.y++;
			if (iColor == 2) node.z++;
			//if (iColor == 3) tree.qnodes[iNode].w++;
			recursiveUpdate(iElement, node.parent);
		}

		void elementClearTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept {
			Tree& tree = trees[iElement];
			index_t iNode = tree.taxonID2nodeID[iTaxon];
			typename Tree::Node& node = tree.nodes[iNode];
			// Assuming x2a, xya or alike are all zeros
			if (iColor == 0) node.x--;
			if (iColor == 1) node.y--;
			if (iColor == 2) node.z--;
			//if (iColor == 3) tree.qnodes[iNode].w--;
			recursiveUpdate(iElement, node.parent);
		}

		void elementClearAndSetTaxonColor(size_t iElement, size_t iTaxon, size_t iColor, size_t jColor) noexcept {
			Tree& tree = trees[iElement];
			index_t iNode = tree.taxonID2nodeID[iTaxon];
			typename Tree::Node& node = tree.nodes[iNode];
			// Assuming x2a, xya or alike are all zeros
			if (iColor == 0) node.x--;
			if (iColor == 1) node.y--;
			if (iColor == 2) node.z--;
			//if (iColor == 3) tree.qnodes[iNode].w--;
			if (jColor == 0) node.x++;
			if (jColor == 1) node.y++;
			if (jColor == 2) node.z++;
			//if (jColor == 3) tree.qnodes[iNode].w++;
			recursiveUpdate(iElement, node.parent);
		}

		score_t elementScore(size_t iElement) const noexcept { return scores[iElement]; }

		/*
		array<score_t, 3> elementQuadripartitionScores(size_t iElement) const noexcept {
			typename Tree::QNode const& qnode = trees[iElement].qnodes.back();
			return { qnode.xy_zw, qnode.xz_yw, qnode.xw_yz };
		}

		void setQuadripartitionMode(bool mode) noexcept { quadMode = mode; }
		*/
	};

	namespace DriverHelper {

		using namespace std;

		namespace Newick {
			void skipWhitespace(string_view& sv) noexcept { while (!sv.empty() && isspace(sv.front())) sv.remove_prefix(1); }

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

			Node(shared_ptr<Node> lc, shared_ptr<Node> rc, bool isReal) : isReal(isReal) {
				children.push_back(lc);
				children.push_back(rc);
			}

			int huffman() noexcept {
				if (children.empty()) {
					common::taxonName2ID[label];
					return 1;
				}
				auto cmp = [](pair<shared_ptr<Node>, int> const& left, pair<shared_ptr<Node>, int> const& right) { return left.second < right.second; };
				priority_queue<pair<shared_ptr<Node>, int>, vector<pair<shared_ptr<Node>, int> >, decltype(cmp)> pq(cmp);
				for (auto& child : children) pq.push({ child, child->huffman() });
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

			template<typename Color> Color::SharedConstData::Element::Node* convert() noexcept {
				using DataNode = typename Color::SharedConstData::Element::Node;
				DataNode* node = new DataNode();

				if (!children.empty()) {
					node->lc.reset(children[0]->convert<Color>());
					node->rc.reset(children[1]->convert<Color>());
				}
				else {
					node->leafID = common::taxonName2ID[label];
				}

				optional<double> len = Newick::parseLength(length);
				if (len) node->length = *len;
				else node->length = 0;
				return node;
			}
		};

		template<typename DataClass> DataClass read() {
			ifstream fin(ARG.get<string>("input"));
			string line;
			vector<shared_ptr<Node> > trees;
			while (getline(fin, line)) {
				string_view sv(line);
				shared_ptr<Node> root(new Node(sv));
				root->huffman();
				trees.push_back(root);
			}
			DataClass data;
			data.nTaxa = common::taxonName2ID.nTaxa();
			for (auto& tree : trees) {
				data.elements.emplace_back();
				data.elements.back().root.reset(tree->convert<typename DataClass::ParentClass>());
			}
			return data;
		}

	};

	template<bool> class Driver : public common::LogInfo
	{
		using string = std::string;

	public:
		using DataClasses = std::variant<typename Color<StepwiseColorDefaultAttributes<double> >::SharedConstData>;
		//using BranchAnnotations = branch_annotation::BranchAnnotations<quadripartition_support::BA_QS, quadripartition_support::BA_LBB, quadripartition_support_spr44::BA_SBB, quadripartition_support::BA_CU_CASTER>;

		//static std::string constexpr DEFAULT_SUPPORT_ANNOTATION() noexcept { return "sbb"; };
		//static std::string constexpr DEFAULT_LENGTH_ANNOTATION() noexcept { return "cu"; };

		static std::pair<string, string> programNames() noexcept {
			return { "astral-iad", "Accurate Species Tree ALgorithm -- InterAnchor Distance (ASTRAL-IAD)" };
		}

		static void addArguments() noexcept {
			//ARG.addArgument('\0', "mode", "integer", "Weighting mode (1: hybrid weighting, 2: support only, 3: length only, 4: unweighted)", 0, true, true, "1");
		}

		static DataClasses getStepwiseColorSharedConstData() noexcept {
			//size_t mode = ARG.get<size_t>("mode");
			return DriverHelper::read<typename Color<StepwiseColorDefaultAttributes<double> >::SharedConstData>();

			//common::LogInfo vlog(-99);
			//vlog.log() << "Error: Invalid mode: " << mode << std::endl;
			//exit(-1);
		}
	};

	class Documentation : public common::DocumentationBase {
	protected:
		string introduction() const noexcept override {
			return R"YOHANETYO(# Accurate Species Tree ALgorithm -- InterAnchor Distance (ASTRAL-IAD)

)YOHANETYO";
		}

		string input() const noexcept override {
			return R"YOHANETYO(# INPUT

)YOHANETYO";
		}

		string programName() const noexcept override { return "astral"; }

		string exampleInput() const noexcept override { return "genetree.nw"; }
	};

};
#endif // !ASTRAL_IAD_HPP

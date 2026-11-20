#ifndef QDIST4P_HPP
#define QDIST4P_HPP

#include "stepwise_colorable.hpp"
#include "quadripartition_support.hpp"

namespace qdist4p {

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
		using ParentClass = Color;

		vector<vector<score_t> > dist;
		index_t nTaxa = 0;

		size_t nElements() const noexcept { return 1; }
	};

private:
	struct Node {
		index_t prev = -1, next = -1;
	};
	index_t X, Y, Z, W;
	index_t cntX = 0, cntY = 0, cntZ = 0, cntW = 0;

	vector<Node> nodes;
	vector<score_t> rowSum;

	SharedConstData const& sharedConstData;
	score_t sxy = 0, sxz = 0, syz = 0, swx = 0, swy = 0, swz = 0;

	void addNext(index_t from, index_t to) {
		index_t next = nodes[from].next;
		nodes[from].next = to;
		nodes[to].next = next;
		nodes[to].prev = from;
		nodes[next].prev = to;
	}

	void remove(index_t id) {
		index_t prev = nodes[id].prev;
		index_t next = nodes[id].next;
		nodes[prev].next = next;
		nodes[next].prev = prev;
		nodes[id].prev = id;
		nodes[id].next = id;
	}

	score_t listSum(index_t lst, index_t iTaxon) const noexcept {
		score_t total = 0;
		vector<score_t> const& row = sharedConstData.dist[iTaxon];
		for (index_t jTaxon = nodes[lst].next; jTaxon != lst; jTaxon = nodes[jTaxon].next) total += row[jTaxon];
		return total;
	}

public:
	void elementSetTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept {
		if (iElement != 0) return;
		vector<vector<score_t> > const& D = sharedConstData.dist;

		for (index_t jTaxon = nodes[X].next; jTaxon != X; jTaxon = nodes[jTaxon].next) rowSum[jTaxon] += D[iTaxon][jTaxon];
		for (index_t jTaxon = nodes[Y].next; jTaxon != Y; jTaxon = nodes[jTaxon].next) rowSum[jTaxon] += D[iTaxon][jTaxon];
		for (index_t jTaxon = nodes[Z].next; jTaxon != Z; jTaxon = nodes[jTaxon].next) rowSum[jTaxon] += D[iTaxon][jTaxon];
		for (index_t jTaxon = nodes[W].next; jTaxon != W; jTaxon = nodes[jTaxon].next) rowSum[jTaxon] += D[iTaxon][jTaxon];
		score_t sumX = listSum(X, iTaxon);
		score_t sumY = listSum(Y, iTaxon);
		score_t sumZ = listSum(Z, iTaxon);
		score_t sumW = listSum(W, iTaxon);
		rowSum[iTaxon] = sumX + sumY + sumZ + sumW;
		if (iColor == 0) {
			++cntX;
			addNext(X, iTaxon);
			sxy += sumY;
			sxz += sumZ;
			swx += sumW;
		} else if (iColor == 1) {
			++cntY;
			addNext(Y, iTaxon);
			sxy += sumX;
			syz += sumZ;
			swy += sumW;
		} else if (iColor == 2) {
			++cntZ;
			addNext(Z, iTaxon);
			sxz += sumX;
			syz += sumY;
			swz += sumW;
		} else if (iColor == 3) {
			++cntW;
			addNext(W, iTaxon);
			swx += sumX;
			swy += sumY;
			swz += sumZ;
		}
		sumX = listSum(X, iTaxon);
		sumY = listSum(Y, iTaxon);
		sumZ = listSum(Z, iTaxon);
		sumW = listSum(W, iTaxon);
	}

	void elementClearTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept {
		if (iElement != 0) return;
		vector<vector<score_t> > const& D = sharedConstData.dist;

		for (index_t jTaxon = nodes[X].next; jTaxon != X; jTaxon = nodes[jTaxon].next) rowSum[jTaxon] -= D[iTaxon][jTaxon];
		for (index_t jTaxon = nodes[Y].next; jTaxon != Y; jTaxon = nodes[jTaxon].next) rowSum[jTaxon] -= D[iTaxon][jTaxon];
		for (index_t jTaxon = nodes[Z].next; jTaxon != Z; jTaxon = nodes[jTaxon].next) rowSum[jTaxon] -= D[iTaxon][jTaxon];
		for (index_t jTaxon = nodes[W].next; jTaxon != W; jTaxon = nodes[jTaxon].next) rowSum[jTaxon] -= D[iTaxon][jTaxon];
		score_t sumX = listSum(X, iTaxon);
		score_t sumY = listSum(Y, iTaxon);
		score_t sumZ = listSum(Z, iTaxon);
		score_t sumW = listSum(W, iTaxon);
		rowSum[iTaxon] = 0;
		if (iColor == 0) {
			--cntX;
			remove(iTaxon);
			sxy -= sumY;
			sxz -= sumZ;
			swx -= sumW;
		}
		else if (iColor == 1) {
			--cntY;
			remove(iTaxon);
			sxy -= sumX;
			syz -= sumZ;
			swy -= sumW;
		}
		else if (iColor == 2) {
			--cntZ;
			remove(iTaxon);
			sxz -= sumX;
			syz -= sumY;
			swz -= sumW;
		}
		else if (iColor == 3) {
			--cntW;
			remove(iTaxon);
			swx -= sumX;
			swy -= sumY;
			swz -= sumZ;
		}
		sumX = listSum(X, iTaxon);
		sumY = listSum(Y, iTaxon);
		sumZ = listSum(Z, iTaxon);
		sumW = listSum(W, iTaxon);
	}
	
	void elementClearAndSetTaxonColor(size_t iElement, size_t iTaxon, size_t iColor, size_t jColor) noexcept {
		if (iElement != 0) return;
		vector<vector<score_t> > const& D = sharedConstData.dist;

		score_t sumX, sumY, sumZ, sumW;
		index_t threshold = (cntX + cntY + cntZ + cntW) / 2;
		if (cntX <= threshold) sumX = listSum(X, iTaxon);
		if (cntY <= threshold) sumY = listSum(Y, iTaxon);
		if (cntZ <= threshold) sumZ = listSum(Z, iTaxon);
		if (cntW <= threshold) sumW = listSum(W, iTaxon);
		if (cntX > threshold) sumX = rowSum[iTaxon] - sumY - sumZ - sumW;
		if (cntY > threshold) sumY = rowSum[iTaxon] - sumX - sumZ - sumW;
		if (cntZ > threshold) sumZ = rowSum[iTaxon] - sumX - sumY - sumW;
		if (cntW > threshold) sumW = rowSum[iTaxon] - sumX - sumY - sumZ;

		if (iColor == 0) {
			--cntX;
			remove(iTaxon);
			sxy -= sumY;
			sxz -= sumZ;
			swx -= sumW;
		}
		else if (iColor == 1) {
			--cntY;
			remove(iTaxon);
			sxy -= sumX;
			syz -= sumZ;
			swy -= sumW;
		}
		else if (iColor == 2) {
			--cntZ;
			remove(iTaxon);
			sxz -= sumX;
			syz -= sumY;
			swz -= sumW;
		}
		else if (iColor == 3) {
			--cntW;
			remove(iTaxon);
			swx -= sumX;
			swy -= sumY;
			swz -= sumZ;
		}
		if (jColor == 0) {
			++cntX;
			addNext(X, iTaxon);
			sxy += sumY;
			sxz += sumZ;
			swx += sumW;
		}
		else if (jColor == 1) {
			++cntY;
			addNext(Y, iTaxon);
			sxy += sumX;
			syz += sumZ;
			swy += sumW;
		}
		else if (jColor == 2) {
			++cntZ;
			addNext(Z, iTaxon);
			sxz += sumX;
			syz += sumY;
			swz += sumW;
		}
		else if (jColor == 3) {
			++cntW;
			addNext(W, iTaxon);
			swx += sumX;
			swy += sumY;
			swz += sumZ;
		}
	}

	score_t elementScore(size_t iElement) const noexcept {
		return -(sxy * cntZ * (cntZ - 1) + sxz * cntY * (cntY - 1) + syz * cntX * (cntX - 1)) / 2;
	}
	
	array<score_t, 3> elementQuadripartitionScores(size_t iElement) const noexcept {
		//return { R.xy_zw, R.xz_yw, R.xw_yz };
		return { -(sxy * cntZ * cntW + swz * cntX * cntY), -(sxz * cntY * cntW + swy * cntX * cntZ), -(syz * cntX * cntW + swx * cntY * cntZ) };
	}

	index_t quartetCount() const noexcept {
		return cntX * cntY * cntZ * cntW;
	}
	
	Color(SharedConstData const& data) noexcept : sharedConstData(data), nodes(data.nTaxa + 4), X(data.nTaxa + 0), Y(data.nTaxa + 1), Z(data.nTaxa + 2), W(data.nTaxa + 3), rowSum(data.nTaxa){
		for (index_t i = 0; i < data.nTaxa + 4; ++i) {
			nodes[i].prev = i;
			nodes[i].next = i;
		}
	}
};

class BranchLength {
public:
	static inline string const SHORT_NAME = "len";
	static inline string const FULL_NAME = "Branch Length";
	static inline string const MAIN_KEY = "len";
	static inline vector<string> const KEYS = { "len" };

	template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE C> class Data {
	public:
		using Color = C;
		using score_t = typename Color::score_t;
		using index_t = typename Color::index_t;
		static Data const ZERO;

	private:
		score_t score;

	public:
		template<typename T> static void staticInitialize(common::Random<T>& random, size_t nElmnts) noexcept {}

		Data(score_t score = Color::ZERO) : score(score) {}

		static array<Data, 3> map(Color& color, size_t iElement) noexcept {
			array<score_t, 3> scores = color.elementQuadripartitionScores(iElement);
			index_t cnt = color.quartetCount();
			Data s0(scores[0] / cnt), s1(scores[1] / cnt), s2(scores[2] / cnt);
			return { s0, s1, s2 };
		}

		static Data reduce(Data const& a, Data const& b) noexcept { return Data(a.score + b.score); }

		static void annotate(common::AnnotatedBinaryTree::Node* node, Data const& left_right, Data const& left_outgroup, Data const& right_outgroup) noexcept {
			node->set("len", (2 * left_right.score - left_outgroup.score - right_outgroup.score) / 4);
		}
	};
};
template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> BranchLength::Data<Color> const	BranchLength::Data<Color>::ZERO = BranchLength::Data<Color>();
using BL = quadripartition_support::BranchAnnotation<BranchLength>;

namespace DriverHelper {

	using namespace std;

	template<typename DataClass> DataClass read() {
		DataClass data;
		ifstream fin(ARG.get<string>("input"));
		int nTaxa;
		fin >> nTaxa;
		data.nTaxa = nTaxa;
		data.dist.resize(nTaxa);
		for (int i = 0; i < nTaxa; ++i) {
			data.dist[i].resize(nTaxa);
		}
		for (int i = 0; i < nTaxa; ++i) {
			string TaxonName;
			fin >> TaxonName;
			int iTaxon = common::taxonName2ID[TaxonName];
			if (iTaxon != i) {
				common::LogInfo log(1);
				log.log() << "Error: Taxon names in the distance matrix abnormal. Duplicated? Rooting not supported for now." << std::endl;
				exit(-1);
			}
			for (int j = 0; j < nTaxa; ++j) {
				fin >> data.dist[i][j];
			}

		}
		return data;
	}

};

template<bool> class Driver : public common::LogInfo
{
	using string = std::string;

public:
	using DataClasses = std::variant<typename Color<StepwiseColorDefaultAttributes<double> >::SharedConstData>;
	using BranchAnnotations = branch_annotation::BranchAnnotations<BL>;

	static std::string constexpr DEFAULT_SUPPORT_ANNOTATION() noexcept { return "none"; };
	static std::string constexpr DEFAULT_LENGTH_ANNOTATION() noexcept { return "len"; };

	static std::pair<string, string> programNames() noexcept {
		return { "qdist4p", "Quartet Distance 4-Point (QDist4P)" };
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
		return R"YOHANETYO(# Quartet Distance 4-Point (QDist4P)

)YOHANETYO";
	}

	string input() const noexcept override {
		return R"YOHANETYO(# INPUT

)YOHANETYO";
	}

	string programName() const noexcept override { return "qdist4p"; }

	string exampleInput() const noexcept override { return "dist.tsv"; }
};

};
#endif // !QDIST4P_HPP

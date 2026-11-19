#ifndef QDIST4P_HPP
#define QDIST4P_HPP

#include "stepwise_colorable.hpp"
#include "quadripartition_support.hpp"
#include "quadripartition_support_spr44.hpp"

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
	index_t X, Y, Z;
	index_t cntX = 0, cntY = 0, cntZ = 0;

	vector<Node> nodes;
	vector<score_t> rowSum;

	SharedConstData const& sharedConstData;
	score_t sxy = 0, sxz = 0, syz = 0;

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
		score_t sumX = listSum(X, iTaxon);
		score_t sumY = listSum(Y, iTaxon);
		score_t sumZ = listSum(Z, iTaxon);
		rowSum[iTaxon] = sumX + sumY + sumZ;
		if (iColor == 0) {
			++cntX;
			addNext(X, iTaxon);
			sxy += sumY;
			sxz += sumZ;
		} else if (iColor == 1) {
			++cntY;
			addNext(Y, iTaxon);
			sxy += sumX;
			syz += sumZ;
		} else if (iColor == 2) {
			++cntZ;
			addNext(Z, iTaxon);
			sxz += sumX;
			syz += sumY;
		}
		sumX = listSum(X, iTaxon);
		sumY = listSum(Y, iTaxon);
		sumZ = listSum(Z, iTaxon);
	}

	void elementClearTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept {
		if (iElement != 0) return;
		vector<vector<score_t> > const& D = sharedConstData.dist;

		for (index_t jTaxon = nodes[X].next; jTaxon != X; jTaxon = nodes[jTaxon].next) rowSum[jTaxon] -= D[iTaxon][jTaxon];
		for (index_t jTaxon = nodes[Y].next; jTaxon != Y; jTaxon = nodes[jTaxon].next) rowSum[jTaxon] -= D[iTaxon][jTaxon];
		for (index_t jTaxon = nodes[Z].next; jTaxon != Z; jTaxon = nodes[jTaxon].next) rowSum[jTaxon] -= D[iTaxon][jTaxon];
		score_t sumX = listSum(X, iTaxon);
		score_t sumY = listSum(Y, iTaxon);
		score_t sumZ = listSum(Z, iTaxon);
		rowSum[iTaxon] = 0;
		if (iColor == 0) {
			--cntX;
			remove(iTaxon);
			sxy -= sumY;
			sxz -= sumZ;
		}
		else if (iColor == 1) {
			--cntY;
			remove(iTaxon);
			sxy -= sumX;
			syz -= sumZ;
		}
		else if (iColor == 2) {
			--cntZ;
			remove(iTaxon);
			sxz -= sumX;
			syz -= sumY;
		}
		sumX = listSum(X, iTaxon);
		sumY = listSum(Y, iTaxon);
		sumZ = listSum(Z, iTaxon);
	}
	
	void elementClearAndSetTaxonColor(size_t iElement, size_t iTaxon, size_t iColor, size_t jColor) noexcept {
		if (iElement != 0) return;
		vector<vector<score_t> > const& D = sharedConstData.dist;

		score_t sumX, sumY, sumZ;
		if (cntX <= cntY + cntZ) sumX = listSum(X, iTaxon);
		if (cntY <= cntX + cntZ) sumY = listSum(Y, iTaxon);
		if (cntZ <= cntX + cntY) sumZ = listSum(Z, iTaxon);
		if (cntX > cntY + cntZ) sumX = rowSum[iTaxon] - sumY - sumZ;
		if (cntY > cntX + cntZ) sumY = rowSum[iTaxon] - sumX - sumZ;
		if (cntZ > cntX + cntY) sumZ = rowSum[iTaxon] - sumX - sumY;

		if (iColor == 0) {
			--cntX;
			remove(iTaxon);
			sxy -= sumY;
			sxz -= sumZ;
		}
		else if (iColor == 1) {
			--cntY;
			remove(iTaxon);
			sxy -= sumX;
			syz -= sumZ;
		}
		else if (iColor == 2) {
			--cntZ;
			remove(iTaxon);
			sxz -= sumX;
			syz -= sumY;
		}
		if (jColor == 0) {
			++cntX;
			addNext(X, iTaxon);
			sxy += sumY;
			sxz += sumZ;
		} else if (jColor == 1) {
			++cntY;
			addNext(Y, iTaxon);
			sxy += sumX;
			syz += sumZ;
		} else if (jColor == 2) {
			++cntZ;
			addNext(Z, iTaxon);
			sxz += sumX;
			syz += sumY;
		}
	}

	score_t elementScore(size_t iElement) const noexcept {
		return -(sxy * cntZ * (cntZ - 1) + sxz * cntY * (cntY - 1) + syz * cntX * (cntX - 1)) / 2;
	}
	/*
	array<score_t, 3> elementQuadripartitionScores(size_t iElement) const noexcept {
		//return { element.wx[2] * element.yz[2], element.wx[1] * element.yz[1], element.wx[0] * element.yz[0] };
	}
	*/
	Color(SharedConstData const& data) noexcept : sharedConstData(data), nodes(data.nTaxa + 3), X(data.nTaxa + 0), Y(data.nTaxa + 1), Z(data.nTaxa + 2), rowSum(data.nTaxa){
		for (index_t i = 0; i < data.nTaxa + 3; ++i) {
			nodes[i].prev = i;
			nodes[i].next = i;
		}
	}
};

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
	using BranchAnnotations = branch_annotation::BranchAnnotations<quadripartition_support::BA_QS, quadripartition_support::BA_LBB, quadripartition_support_spr44::BA_SBB>;

	static std::string constexpr DEFAULT_SUPPORT_ANNOTATION() noexcept { return "sbb"; };
	static std::string constexpr DEFAULT_LENGTH_ANNOTATION() noexcept { return "none"; };

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

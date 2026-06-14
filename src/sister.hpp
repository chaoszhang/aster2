#ifndef SISTER_HPP
#define SISTER_HPP

#include "stepwise_colorable.hpp"

namespace sister {

using std::size_t;
using std::views::iota;
using std::vector;
using std::array;
using std::string;

class Color {
public:
	using score_t = long long;
	using index_t = long long;

	static inline bool constexpr IS_ROOTED = false;
	static inline score_t constexpr ZERO = 0LL;
	static inline score_t constexpr EPSILON = 0LL;

	struct SharedConstData {
		using ParentClass = Color;

		struct Element {
			vector<index_t> alleleType; // alleleType[iTaxon] -> allele type of iTaxon, -1 if missing
			index_t nType = 0; // number of allele types
		};

		vector<Element> elements;
		index_t nTaxa = 0;

		size_t nElements() const noexcept { return elements.size(); }
	};

	struct Element {
		vector<array<index_t, 4> > cnts; // cnts[iAlleleType][iColor] -> count of iAlleleType in iColor
		array<index_t, 3> xx{}, yz{}, xxyz{}, wx{};
	};

private:
	SharedConstData const& sharedConstData;
	vector<Element> elements;

	template<bool isSet> inline void elementSetOrClearTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept {
		typename SharedConstData::Element const& cElement = sharedConstData.elements[iElement];
		Element& element = elements[iElement];

		index_t iType = cElement.alleleType[iTaxon];
		if (iType == -1 || cElement.nType <= 1) return;
		
		index_t& cnt = element.cnts[iType][iColor];
		if (cElement.nType == 2) {
			if constexpr (isSet) cnt++;
			else cnt--;
		}
		else {
			element.xx[0] -= element.cnts[iType][0] * (element.cnts[iType][0] - 1);
			element.xx[1] -= element.cnts[iType][1] * (element.cnts[iType][1] - 1);
			element.xx[2] -= element.cnts[iType][2] * (element.cnts[iType][2] - 1);
			element.yz[0] -= element.cnts[iType][1] * element.cnts[iType][2];
			element.yz[1] -= element.cnts[iType][2] * element.cnts[iType][0];
			element.yz[2] -= element.cnts[iType][0] * element.cnts[iType][1];
			element.wx[0] -= element.cnts[iType][0] * element.cnts[iType][3];
			element.wx[1] -= element.cnts[iType][1] * element.cnts[iType][3];
			element.wx[2] -= element.cnts[iType][2] * element.cnts[iType][3];
			element.xxyz[0] -= element.cnts[iType][0] * (element.cnts[iType][0] - 1) * element.cnts[iType][1] * element.cnts[iType][2];
			element.xxyz[1] -= element.cnts[iType][1] * (element.cnts[iType][1] - 1) * element.cnts[iType][2] * element.cnts[iType][0];
			element.xxyz[2] -= element.cnts[iType][2] * (element.cnts[iType][2] - 1) * element.cnts[iType][0] * element.cnts[iType][1];
			if constexpr (isSet) cnt++;
			else cnt--;
			element.xx[0] += element.cnts[iType][0] * (element.cnts[iType][0] - 1);
			element.xx[1] += element.cnts[iType][1] * (element.cnts[iType][1] - 1);
			element.xx[2] += element.cnts[iType][2] * (element.cnts[iType][2] - 1);
			element.yz[0] += element.cnts[iType][1] * element.cnts[iType][2];
			element.yz[1] += element.cnts[iType][2] * element.cnts[iType][0];
			element.yz[2] += element.cnts[iType][0] * element.cnts[iType][1];
			element.wx[0] += element.cnts[iType][0] * element.cnts[iType][3];
			element.wx[1] += element.cnts[iType][1] * element.cnts[iType][3];
			element.wx[2] += element.cnts[iType][2] * element.cnts[iType][3];
			element.xxyz[0] += element.cnts[iType][0] * (element.cnts[iType][0] - 1) * element.cnts[iType][1] * element.cnts[iType][2];
			element.xxyz[1] += element.cnts[iType][1] * (element.cnts[iType][1] - 1) * element.cnts[iType][2] * element.cnts[iType][0];
			element.xxyz[2] += element.cnts[iType][2] * (element.cnts[iType][2] - 1) * element.cnts[iType][0] * element.cnts[iType][1];
		}
	}

public:
	void elementSetTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept {
		elementSetOrClearTaxonColor<true>(iElement, iTaxon, iColor);
	}

	void elementClearTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept {
		elementSetOrClearTaxonColor<false>(iElement, iTaxon, iColor);
	}

	score_t elementScore(size_t iElement) const noexcept {
		typename SharedConstData::Element const& cElement = sharedConstData.elements[iElement];
		Element const& element = elements[iElement];

		if (cElement.nType <= 1) return 0;

		if (cElement.nType == 2) {
			return element.cnts[0][0] * (element.cnts[0][0] - 1) * element.cnts[1][1] * element.cnts[1][2]
				+ element.cnts[0][1] * (element.cnts[0][1] - 1) * element.cnts[1][2] * element.cnts[1][0]
				+ element.cnts[0][2] * (element.cnts[0][2] - 1) * element.cnts[1][0] * element.cnts[1][1]
				+ element.cnts[1][0] * (element.cnts[1][0] - 1) * element.cnts[0][1] * element.cnts[0][2]
				+ element.cnts[1][1] * (element.cnts[1][1] - 1) * element.cnts[0][2] * element.cnts[0][0]
				+ element.cnts[1][2] * (element.cnts[1][2] - 1) * element.cnts[0][0] * element.cnts[0][1]; 
		}
		else {
			return element.xx[0] * element.yz[0] - element.xxyz[0]
				+ element.xx[1] * element.yz[1] - element.xxyz[1]
				+ element.xx[2] * element.yz[2] - element.xxyz[2];
		}
	}

	array<score_t, 3> elementQuadripartitionScores(size_t iElement) const noexcept {
		typename SharedConstData::Element const& cElement = sharedConstData.elements[iElement];
		Element const& element = elements[iElement];

		if (cElement.nType <= 1) return { 0, 0, 0 };

		if (cElement.nType == 2) {
			return { element.cnts[0][0] * element.cnts[0][1] * element.cnts[1][2] * element.cnts[1][3] + element.cnts[1][0] * element.cnts[1][1] * element.cnts[0][2] * element.cnts[0][3],
				element.cnts[0][0] * element.cnts[0][2] * element.cnts[1][1] * element.cnts[1][3] + element.cnts[1][0] * element.cnts[1][2] * element.cnts[0][1] * element.cnts[0][3],
				element.cnts[0][0] * element.cnts[0][3] * element.cnts[1][1] * element.cnts[1][2] + element.cnts[1][0] * element.cnts[1][3] * element.cnts[0][1] * element.cnts[0][2] };
		}
		else {
			return { element.wx[2] * element.yz[2], element.wx[1] * element.yz[1], element.wx[0] * element.yz[0] };
		}
	}

	Color(SharedConstData const& data) noexcept : sharedConstData(data), elements(data.nElements()) {
		for (index_t iElement : iota(index_t(0), index_t(data.nElements()))) {
			typename SharedConstData::Element const& cElement = sharedConstData.elements[iElement];
			Element& element = elements[iElement];
			element.cnts.resize(cElement.nType);
		}
	}
};

namespace DriverHelper {

using namespace std;

typename Color::SharedConstData read() noexcept {
	typename Color::SharedConstData data;
	ifstream fin(ARG.get<string>("input"));
	int nTaxa, nElements;
	fin >> nTaxa >> nElements;
	data.elements.resize(nElements);
	data.nTaxa = nTaxa;
	for (int iTaxon = 0; iTaxon < nTaxa; ++iTaxon) {
		string TaxonName;
		fin >> TaxonName;
		common::taxonName2ID[TaxonName];
		for (int iElement = 0; iElement < nElements; ++iElement) {
			int alleleType;
			fin >> alleleType;
			data.elements[iElement].alleleType.push_back(alleleType);
			if (alleleType >= data.elements[iElement].nType) data.elements[iElement].nType = alleleType + 1;
		}
	}
	return data;
}

};

template<bool> class Driver : public common::LogInfo
{
	using string = std::string;

public:
	using DataClasses = std::variant<Color::SharedConstData>;

	static std::pair<string, string> programNames() noexcept {
		return { "sister", "Synteny-Induced Species Tree EstimatoR" };
	}

	static void addArguments() noexcept {
		//ARG.addArgument('\0', "chunk", "integer", "The maximum number of sites in each local aligment block for parameter estimation", 0, true, true, "10000");
	}

	static DataClasses getStepwiseColorSharedConstData() noexcept {
		return DriverHelper::read();
	}
};

class Documentation : public common::DocumentationBase {
protected:
	string introduction() const noexcept override {
		return R"YOHANETYO(# Synteny-Induced Species Tree EstimatoR (SISTER)
)YOHANETYO";
	}

	string input() const noexcept override {
		return R"YOHANETYO(# INPUT
)YOHANETYO";
	}

	string programName() const noexcept override { return "sister"; }

	string exampleInput() const noexcept override { return "alignment.phylip"; }
};


};

#endif // !SISTER_HPP

#ifndef CASTER_HPP
#define CASTER_HPP

#include "common.hpp"
#include "stepwise_colorable.hpp"
#include "alignment_utilities.hpp"

namespace caster{

using std::size_t;
using std::views::iota;
using std::vector;
using std::array;
	
template<class Attributes> concept STEPWISE_COLOR_ATTRIBUTES = requires
{
    requires std::integral<typename Attributes::score_t> || std::floating_point<typename Attributes::score_t>;
	requires std::integral<typename Attributes::cnt_t> || std::floating_point<typename Attributes::cnt_t>;
	requires std::integral<typename Attributes::cnt2_t> || std::floating_point<typename Attributes::cnt2_t>;
	requires std::integral<typename Attributes::cnt4_t> || std::floating_point<typename Attributes::cnt4_t>;
	requires std::integral<typename Attributes::index_t>;
	{ Attributes::ZERO } -> std::convertible_to<typename Attributes::score_t>;
	{ Attributes::EPSILON } -> std::convertible_to<typename Attributes::score_t>;
};

struct StepwiseColorDefaultAttributes{
	using score_t = double;
	using cnt_taxon_t = unsigned char;
	using cnt_t = unsigned short;
	using cnt2_t = unsigned int;
	using cnt4_t = long long;
	using index_t = long long;
	static inline score_t constexpr ZERO = 0;
	static inline score_t constexpr EPSILON = 1e-3;
};

template<STEPWISE_COLOR_ATTRIBUTES Attributes> class Color{
public:
	using score_t = Attributes::score_t;
	using cnt_taxon_t = Attributes::cnt_taxon_t;
	using cnt_t = Attributes::cnt_t;
	using cnt2_t = Attributes::cnt2_t;
	using cnt4_t = Attributes::cnt4_t;
	using index_t = Attributes::index_t;
	static inline score_t constexpr ZERO = Attributes::ZERO;
	static inline score_t constexpr EPSILON = Attributes::EPSILON;
	
	struct SharedConstData{
		using ParentClass = Color<Attributes>;

		struct Element{
			index_t iGenomePosBegin = 0;
			index_t nPos = 0;
			vector<vector<array<cnt_taxon_t, 4> > > cnts; // cnts[iRow][iPos][iNucleotide] -> count
			vector<index_t> taxon2row; // taxon2row[iTaxon] -> iRow in cnts
			array<score_t, 4> eqFreqs{}; // eqFreqs[iNucleotide]

			bool hasTaxon(size_t iTaxon) const noexcept{
				return iTaxon < taxon2row.size() && taxon2row[iTaxon] != -1;
			}
		};

		vector<Element> elements;

		size_t nElements = 0;
		index_t nGenomePos = 0;
    };

private:
	SharedConstData const* const sharedConstData;
    vector<array<array<cnt_t, 4>, 3> > colorCnts; // colorCnts[iGenomePos][iColor][iNucleotide] -> count

	template<bool isSet> inline void elementSetOrClearTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept{
		typename SharedConstData::Element const& element = sharedConstData->elements[iElement];
		if (!element.hasTaxon(iTaxon)) return;
		index_t iRow = element.taxon2row[iTaxon];
		for (index_t iPos : iota(0, element.nPos)){
			for (index_t iNucleotide : iota(0, 4)) {
				cnt_t& colorCnt = colorCnts[element.iGenomePosBegin + iPos][iColor][iNucleotide];
				cnt_t cnt = element.cnts[iRow][iPos][iNucleotide];
				if constexpr (isSet) colorCnt += cnt;
				else colorCnt -= cnt;
			}
		}
	}
	
	inline static cnt4_t XXYY(cnt4_t x0, cnt4_t x1, cnt4_t x2, cnt4_t y0, cnt4_t y1, cnt4_t y2) noexcept{
		return x0 * (x0 - 1) * y1 * y2 + x1 * (x1 - 1) * y2 * y0 + x2 * (x2 - 1) * y0 * y1
			 + y0 * (y0 - 1) * x1 * x2 + y1 * (y1 - 1) * x2 * x0 + y2 * (y2 - 1) * x0 * x1;
	}

	inline static score_t scorePos(array<array<cnt_t, 4>, 3> const &cnt, array<score_t, 4> const &pi) noexcept{
		//lst = simplify([sABCD(R, R, Y, Y); sABCD(A, A, Y, Y); sABCD(C, C, R, R); sABCD(A, A, C, C)])
		//sol = [Pa^2*Pc^2; -Pc^2*(Pa + Pg)^2; -Pa^2*(Pc + Pt)^2; (Pa + Pg)^2*(Pc + Pt)^2]
		//lst2 = simplify([sABCD(Y, Y, R, R); sABCD(Y, Y, A, A); sABCD(R, R, C, C); sABCD(C, C, A, A)])

		// (Pa^2+Pg^2)*(Pc^2+Pt^2)*sABCD(R, R, Y, Y)
		// -Pr^2*(Pc^2+Pt^2)*sABCD(A, A, Y, Y) -Pr^2*(Pc^2+Pt^2)*sABCD(G, G, Y, Y) -(Pa^2+Pg^2)*Py^2*sABCD(C, C, R, R) -(Pa^2+Pg^2)*Py^2*sABCD(T, T, R, R)
		// +Pr^2*Py^2*sABCD(A, A, C, C) +Pr^2*Py^2*sABCD(G, G, C, C) +Pr^2*Py^2*sABCD(A, A, T, T) +Pr^2*Py^2*sABCD(G, G, T, T)
		
		score_t const A = pi[0], C = pi[1], G = pi[2], T = pi[3];
		score_t const R = A + G, Y = C + T, R2 = A * A + G * G, Y2 = C * C + T * T;
		cnt4_t const a0 = cnt[0][0], c0 = cnt[0][1], g0 = cnt[0][2], t0 = cnt[0][3], r0 = a0 + g0, y0 = c0 + t0;
		cnt4_t const a1 = cnt[1][0], c1 = cnt[1][1], g1 = cnt[1][2], t1 = cnt[1][3], r1 = a1 + g1, y1 = c1 + t1;
		cnt4_t const a2 = cnt[2][0], c2 = cnt[2][1], g2 = cnt[2][2], t2 = cnt[2][3], r2 = a2 + g2, y2 = c2 + t2;
		
		cnt4_t const rryy = XXYY(r0, r1, r2, y0, y1, y2);

		cnt4_t const aayy = XXYY(a0, a1, a2, y0, y1, y2);
		cnt4_t const ggyy = XXYY(g0, g1, g2, y0, y1, y2);
		cnt4_t const rrcc = XXYY(r0, r1, r2, c0, c1, c2);
		cnt4_t const rrtt = XXYY(r0, r1, r2, t0, t1, t2);
		
		cnt4_t const aacc = XXYY(a0, a1, a2, c0, c1, c2);
		cnt4_t const aatt = XXYY(a0, a1, a2, t0, t1, t2);
		cnt4_t const ggcc = XXYY(g0, g1, g2, c0, c1, c2);
		cnt4_t const ggtt = XXYY(g0, g1, g2, t0, t1, t2);
		
		return rryy * R2 * Y2 - (aayy + ggyy) * (R * R) * Y2 - (rrcc + rrtt) * R2 * (Y * Y)
			 + (aacc + aatt + ggcc + ggtt) * (R * R) * (Y * Y);
	}
	
public:
	void elementSetTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept{
		elementSetOrClearTaxonColor<true>(iElement, iTaxon, iColor);
	}
	
	void elementClearTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept{
		elementSetOrClearTaxonColor<false>(iElement, iTaxon, iColor);
	}
	
	score_t elementScore(size_t iElement) const noexcept{
		index_t iGenomePosBegin = sharedConstData->elements[iElement].iGenomePosBegin;
		index_t nPos = sharedConstData->elements[iElement].nPos;
		typename SharedConstData::Element const& element = sharedConstData->elements[iElement];

		score_t res = 0;
		for (index_t iPos : iota(0, nPos)){
			res += scorePos(colorCnts[iGenomePosBegin + iPos], element.eqFreqs);
		}
		return res;
	}
	
	Color(SharedConstData const* const data) noexcept: sharedConstData(data), colorCnts(data->nGenomePos){}
};

namespace DriverHelper {

using namespace std;

static array<int, 4> add(const array<int, 4>& a, const array<int, 4>& b) {
    array<int, 4> result;
    for (int j = 0; j < 4; j++) {
        result[j] = a[j] + b[j];
    }
    return result;
}

static int sum(const array<int, 4>& cnt) {
    int result = 0;
    for (int j = 0; j < 4; j++) {
        result += cnt[j];
    }
    return result;
}

template<typename DataClasses> DataClasses read() {
	using DataClass0 = std::variant_alternative_t<0, DataClasses>;
	DataClass0 sharedConstData;

	const string& file = ARG.get<string>("input");
	aligment_utilities::AlignmentParser AP(file, 1), AP2(file, 2);
    while (AP.nextAlignment()) {
        AP2.nextAlignment();
        size_t nSites = AP.getLength();
		size_t chunkMaxSize = ARG.get<size_t>("chunk");
		size_t nChunk = (nSites + chunkMaxSize - 1) / chunkMaxSize;
        vector<vector<size_t> > sites(nChunk);
        vector<array<double, 4> > eqfreq;
		size_t iElementBegin = sharedConstData.elements.size();
		unordered_map<size_t, size_t> taxon2row;
        {
            vector<array<unsigned short, 4> > freq;
            freq.resize(AP.getLength());
            while (AP.nextSeq()) {
                size_t iTaxon = common::taxonName2ID[AP.getName()];
				if (!taxon2row.count(iTaxon)) taxon2row[iTaxon] = taxon2row.size();
                string seq = AP.getSeq();
                for (size_t i = 0; i < seq.size(); i++) {
                    switch (seq[i]) {
						case 'A': freq[i][0]++; break;
						case 'C': freq[i][1]++; break;
						case 'G': freq[i][2]++; break;
						case 'T': freq[i][3]++; break;
                    }
                }
            }
            for (size_t i = 0; i < nChunk; i++) {
				size_t s = i * nSites / nChunk, t = (i + 1) * nSites / nChunk;
				size_t sumFreq[4] = {};
                for (size_t j = s; j < t; j++) {
                    for (int k = 0; k < 4; k++) {
                        sumFreq[k] += freq[j][k];
                    }
#ifdef CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH
                    sites[i].push_back(j);
#else
                    if (freq[j][0] + freq[j][2] >= 2 && freq[j][1] + freq[j][3] >= 2) sites[i].push_back(j);
#endif
                }
                double total = sumFreq[0] + sumFreq[1] + sumFreq[2] + sumFreq[3];
                if (total > 0) eqfreq.push_back({ sumFreq[0] / total, sumFreq[1] / total, sumFreq[2] / total, sumFreq[3] / total });
                else eqfreq.push_back({ 0.25, 0.25, 0.25, 0.25 });
            }
        }
		for (size_t i = 0; i < nChunk; i++) {
			typename DataClass0::Element element;
			element.iGenomePosBegin = sharedConstData.nGenomePos;
			element.nPos = sites[i].size();
			element.cnts.resize(taxon2row.size(), vector<array<typename DataClass0::ParentClass::cnt_taxon_t, 4> >(element.nPos));
			element.taxon2row.resize(common::taxonName2ID.nTaxa(), -1);
			element.eqFreqs = eqfreq[i];
			sharedConstData.elements.push_back(element);
			sharedConstData.nElements++;
			sharedConstData.nGenomePos += element.nPos;
		}
        while (AP2.nextSeq()) {
			size_t iTaxon = common::taxonName2ID[AP2.getName()];
			size_t iRow = taxon2row[iTaxon];
            string seq = AP2.getSeq();
			for (size_t iChunk : iota((size_t) 0, nChunk)) {
				typename DataClass0::Element &element = sharedConstData.elements[iElementBegin + iChunk];
				element.taxon2row[iTaxon] = iRow;
				for (size_t iPos : iota((size_t) 0, sites[iChunk].size())) {
					switch (seq[sites[iChunk][iPos]]) {
						case 'A': element.cnts[iRow][iPos][0]++; break;
						case 'C': element.cnts[iRow][iPos][1]++; break;
						case 'G': element.cnts[iRow][iPos][2]++; break;
						case 'T': element.cnts[iRow][iPos][3]++; break;
					}
				}
            }
        }
    }
	return sharedConstData;
}

};

class Driver : public common::LogInfo
{
	using string = std::string;

public:
	using DataClasses = std::variant<Color<StepwiseColorDefaultAttributes>::SharedConstData>;
	
	static std::pair<string, string> programNames() {
		return { "caster-site", "Coalescence-aware Alignment-based Species Tree EstimatoR (Site)" };
	}

	static void addArguments() {
		ARG.addArgument('\0', "chunk", "integer", "The maximum number of sites in each local aligment block for parameter estimation", 0, true, true, "10000");
	}

	static DataClasses getStepwiseColorSharedConstData(){
		return DriverHelper::read<DataClasses>();
	}
};



};
#endif

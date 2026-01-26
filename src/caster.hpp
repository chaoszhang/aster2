#ifndef CASTER_HPP
#define CASTER_HPP

#include "std_include.hpp"
#include "stepwise_colorable.hpp"

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
};

struct StepwiseColorDefaultAttributes{
	using score_t = double;
	using cnt_taxon_t = unsigned char;
	using cnt_t = unsigned short;
	using cnt2_t = unsigned int;
	using cnt4_t = unsigned long long;
	using index_t = long long;
	static inline score_t constexpr ZERO = 0;
};

template<STEPWISE_COLOR_ATTRIBUTES Attributes> class StepwiseColor{
public:
	using score_t = Attributes::score_t;
	using cnt_taxon_t = Attributes::cnt_taxon_t;
	using cnt_t = Attributes::cnt_t;
	using cnt2_t = Attributes::cnt2_t;
	using cnt4_t = Attributes::cnt4_t;
	using index_t = Attributes::index_t;
	static inline score_t constexpr ZERO = Attributes::ZERO;
	
	struct SharedConstData{
		size_t nElements = 0;
        index_t nGenomePos = 0;
		vector<array<vector<cnt_taxon_t>, 4> > cnts; // cnts[iTaxon][iNucleotide][iGenomePos] -> count
		vector<array<index_t, 2> > elementGenomePosRanges; // elementGenomePosRanges[iElement] -> (iGenomePosBegin, iGenomePosEnd)
		vector<array<score_t, 4> > elementEqFreqs; // cnts[iElement][iNucleotide] -> eqFreq
    };

private:
	SharedConstData const* const sharedConstData;
    vector<array<array<cnt_t, 4>, 3> > colorCnts; // colorCnts[iGenomePos][iColor][iNucleotide] -> count
	
	template<bool isSet> inline void elementSetOrClearTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept{
		for (index_t iNucleotide : iota(0, 4)){
			const array<index_t, 2> &elementRange = sharedConstData->elementGenomePosRanges[iElement];
			for (index_t iGenomePos : iota(elementRange[0], elementRange[1])){
				cnt_t& colorCnt = colorCnts[iGenomePos][iColor][iNucleotide];
				cnt_t cnt = sharedConstData->cnts[iTaxon][iNucleotide][iGenomePos];
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
		array<index_t, 2> const &elementRange = sharedConstData->elementGenomePosRanges[iElement];
		score_t res = 0;
		for (index_t iGenomePos : iota(elementRange[0], elementRange[1])){
			res += scorePos(colorCnts[iGenomePos], sharedConstData->elementEqFreqs[iGenomePos]);
		}
		return res;
	}
	
	StepwiseColor(SharedConstData const* const data) noexcept: sharedConstData(data), colorCnts(data->nGenomePos){}
};

};
#endif

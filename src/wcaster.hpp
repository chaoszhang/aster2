#ifndef WCASTER_HPP
#define WCASTER_HPP

#include "stepwise_colorable.hpp"
#include "quadripartition_support.hpp"
#include "quadripartition_support_spr44.hpp"
#include "alignment_utilities.hpp"

namespace wcaster{

using std::size_t;
using std::views::iota;
using std::vector;
using std::array;
using std::string;
	
namespace DriverHelper {
	template<typename DataClasses> DataClasses read();
};

template<class Attributes> concept STEPWISE_COLOR_ATTRIBUTES = requires
{
    requires std::integral<typename Attributes::score_t> || std::floating_point<typename Attributes::score_t>;
	requires std::integral<typename Attributes::cnt_t> || std::floating_point<typename Attributes::cnt_t>;
	requires std::integral<typename Attributes::cnt4_t> || std::floating_point<typename Attributes::cnt4_t>;
	requires std::integral<typename Attributes::index_t>;
	{ Attributes::ZERO } -> std::convertible_to<typename Attributes::score_t>;
	{ Attributes::EPSILON } -> std::convertible_to<typename Attributes::score_t>;
};

// ══════════ WCASTER-MOD ① BEGIN: cnt_t/cnt4_t fixed to double (single cnt_taxon_type param) ══════════
template<typename cnt_taxon_type = unsigned char> struct StepwiseColorDefaultAttributes {
	using score_t = double;
	using cnt_taxon_t = cnt_taxon_type;
	using cnt_t = double;
	using cnt4_t = double;
	using index_t = long long;
	static inline score_t constexpr ZERO = 0;
	static inline score_t constexpr EPSILON = 1e-3;
};
// ══════════ WCASTER-MOD ① END ══════════

ChangeLog logColor("Color",
	"2026-02-02", "Chao Zhang", "Supporting quadripartiton", "minor",
	"2026-09-04", "Zuizhi Chen", "WCASTER: quartet scoring with per-chunk weighting (colorWeight/colorPairWeight)", "minor");

template<STEPWISE_COLOR_ATTRIBUTES Attributes> class Color{
	using cnt_taxon_t = Attributes::cnt_taxon_t;
	using cnt_t = Attributes::cnt_t;
	using cnt4_t = Attributes::cnt4_t;
	using index_t = Attributes::index_t;

public:
	using score_t = Attributes::score_t;
	static inline bool constexpr IS_ROOTED = false;
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
			// ══════════ WCASTER-MOD ② BEGIN ══════════
			vector<cnt_t> speciesWeights;        // WCASTER: per-species-row sequence-similarity weight
			// ══════════ WCASTER-MOD ② END ══════════

			bool hasTaxon(size_t iTaxon) const noexcept{
				return iTaxon < taxon2row.size() && taxon2row[iTaxon] != -1;
			}
		};

		vector<Element> elements;
		index_t nGenomePos = 0;

		size_t nElements() const noexcept { return elements.size(); }
    };

private:
	SharedConstData const& sharedConstData;
	// ══════════ WCASTER-MOD ③ BEGIN: colorCnts → colorWeight + colorPairWeight ══════════
	vector<array<array<cnt_t, 4>, 4> > colorWeight;     // colorWeight[iGenomePos][iColor][iNucleotide] -> weighted sum Σ w_i·c_i
	vector<array<array<cnt_t, 4>, 4> > colorPairWeight; // colorPairWeight[iGenomePos][iColor][iNucleotide] -> unordered pair sum Σ_{i<j} w_i·w_j (same nucleotide)
	// ══════════ WCASTER-MOD ③ END ══════════

	// ══════════ WCASTER-MOD ④ BEGIN: weighted set/clear maintaining 1st- and 2nd-order terms ══════════
	template<bool isSet> inline void elementSetOrClearTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept{
		typename SharedConstData::Element const& element = sharedConstData.elements[iElement];
		if (!element.hasTaxon(iTaxon)) return;
		index_t iRow = element.taxon2row[iTaxon];
		index_t iPosBegin = element.iGenomePosBegin;

		cnt_t w = element.speciesWeights[iRow]; // WCASTER: per-species-row weight
		for (index_t iPos : iota((index_t)0, element.nPos)){
			index_t gPos = iPosBegin + iPos;
			for (index_t iNuc : iota((index_t)0, (index_t)4)) {
				cnt_taxon_t cnt = element.cnts[iRow][iPos][iNuc];
				if (cnt == 0) continue;
				cnt_t dw = (cnt_t)cnt * w;
				if constexpr (isSet) {
					// add: pair first (uses OLD colorWeight), then weight
					colorPairWeight[gPos][iColor][iNuc] += colorWeight[gPos][iColor][iNuc] * dw;
					colorWeight[gPos][iColor][iNuc] += dw;
				} else {
					// remove: weight first, then pair (uses NEW colorWeight)
					colorWeight[gPos][iColor][iNuc] -= dw;
					colorPairWeight[gPos][iColor][iNuc] -= colorWeight[gPos][iColor][iNuc] * dw;
				}
			}
		}
	}
	// ══════════ WCASTER-MOD ④ END ══════════
	
	// ══════════ WCASTER-MOD ⑤ BEGIN: weighted XXYY over the 3 colours (naming follows wtrial)
	// pw(): unordered pair-weight of nucleotide set {n1,n2} within one colour.
	//   n1==n2 (same nucleotide): cp[col][n1] (colourPairWeight tracks same-nuc unordered pairs)
	//   n1!=n2 (mixed set, e.g. A+G): AA + GG pairs plus cross products cw[n1]*cw[n2]
	inline static cnt_t pw(array<array<cnt_t, 4>, 4> const& cw, array<array<cnt_t, 4>, 4> const& cp, int col, int n1, int n2) noexcept{
		if (n1 == n2) return cp[col][n1];
		return cp[col][n1] + cp[col][n2] + cw[col][n1] * cw[col][n2];
	}

	// Integer CASTER uses x*(x-1) (ORDERED pairs); weighted analog = 2 x (unordered pair-weight).
	// x00,x11,x22 (resp. y00,y11,y22) = unordered pair-weight of the X (resp. Y) letter-set in colour 0/1/2.
	inline static cnt_t XXYY(cnt_t x0, cnt_t x1, cnt_t x2, cnt_t x00, cnt_t x11, cnt_t x22,
	                             cnt_t y0, cnt_t y1, cnt_t y2, cnt_t y00, cnt_t y11, cnt_t y22) noexcept{
		return 2 * (x00 * y1 * y2 + x11 * y2 * y0 + x22 * y0 * y1) +
		       2 * (y00 * x1 * x2 + y11 * x2 * x0 + y22 * x0 * x1);
	}
	// ══════════ WCASTER-MOD ⑤ END ══════════

	// ══════════ WCASTER-MOD ⑥ BEGIN: weighted scorePos over 3 colour groups (symmetric) ══════════
	inline static score_t scorePos(array<array<cnt_t, 4>, 4> const& cw, // colorWeight (weighted sums)
	                               array<array<cnt_t, 4>, 4> const& cp, // colorPairWeight (pair sums)
	                               array<score_t, 4> const& pi) noexcept{
		// unordered pair-weights per colour: purine (A+G) and pyrimidine (C+T) pairs
		cnt_t const r00 = pw(cw, cp, 0, 0, 2), r11 = pw(cw, cp, 1, 0, 2), r22 = pw(cw, cp, 2, 0, 2);
		cnt_t const y00 = pw(cw, cp, 0, 1, 3), y11 = pw(cw, cp, 1, 1, 3), y22 = pw(cw, cp, 2, 1, 3);
		// single-nucleotide unordered pairs per colour
		cnt_t const a00 = cp[0][0], a11 = cp[1][0], a22 = cp[2][0];
		cnt_t const g00 = cp[0][2], g11 = cp[1][2], g22 = cp[2][2];
		cnt_t const c00 = cp[0][1], c11 = cp[1][1], c22 = cp[2][1];
		cnt_t const t00 = cp[0][3], t11 = cp[1][3], t22 = cp[2][3];

		cnt_t const a0 = cw[0][0], c0 = cw[0][1], g0 = cw[0][2], t0 = cw[0][3];
		cnt_t const a1 = cw[1][0], c1 = cw[1][1], g1 = cw[1][2], t1 = cw[1][3];
		cnt_t const a2 = cw[2][0], c2 = cw[2][1], g2 = cw[2][2], t2 = cw[2][3];

		score_t const A = pi[0], C = pi[1], G = pi[2], T = pi[3];
		score_t const R = A + G, Y = C + T, R2 = A * A + G * G, Y2 = C * C + T * T;
		cnt_t const r0 = a0 + g0, y0 = c0 + t0;
		cnt_t const r1 = a1 + g1, y1 = c1 + t1;
		cnt_t const r2 = a2 + g2, y2 = c2 + t2;

		cnt_t const rryy = XXYY(r0, r1, r2, r00, r11, r22, y0, y1, y2, y00, y11, y22);
		cnt_t const aayy = XXYY(a0, a1, a2, a00, a11, a22, y0, y1, y2, y00, y11, y22);
		cnt_t const ggyy = XXYY(g0, g1, g2, g00, g11, g22, y0, y1, y2, y00, y11, y22);
		cnt_t const rrcc = XXYY(r0, r1, r2, r00, r11, r22, c0, c1, c2, c00, c11, c22);
		cnt_t const rrtt = XXYY(r0, r1, r2, r00, r11, r22, t0, t1, t2, t00, t11, t22);
		cnt_t const aacc = XXYY(a0, a1, a2, a00, a11, a22, c0, c1, c2, c00, c11, c22);
		cnt_t const aatt = XXYY(a0, a1, a2, a00, a11, a22, t0, t1, t2, t00, t11, t22);
		cnt_t const ggcc = XXYY(g0, g1, g2, g00, g11, g22, c0, c1, c2, c00, c11, c22);
		cnt_t const ggtt = XXYY(g0, g1, g2, g00, g11, g22, t0, t1, t2, t00, t11, t22);

		return (score_t)(rryy * R2 * Y2 - (aayy + ggyy) * (R * R) * Y2 - (rrcc + rrtt) * R2 * (Y * Y)
		     + (aacc + aatt + ggcc + ggtt) * (R * R) * (Y * Y));
	}
	// ══════════ WCASTER-MOD ⑥ END ══════════
	
	// ══════════ WCASTER-MOD ⑦ BEGIN: quadXXYY/quadPos — formula unchanged; cnt_t=cnt4_t=double now,
	// quadripartition takes one individual per colour (single sums, no pair terms) → fed by colorWeight. ══════════
	inline static cnt4_t quadXXYY(cnt4_t x0, cnt4_t x1, cnt4_t x2, cnt4_t x3, cnt4_t y0, cnt4_t y1, cnt4_t y2, cnt4_t y3) noexcept {
		return x0 * x1 * y2 * y3 + y0 * y1 * x2 * x3;
	}

	inline static score_t quadPos(array<cnt_t, 4> const& cnt0, array<cnt_t, 4> const& cnt1,
		array<cnt_t, 4> const& cnt2, array<cnt_t, 4> const& cnt3, array<score_t, 4> const& pi) noexcept {
		//lst = simplify([sABCD(R, R, Y, Y); sABCD(A, A, Y, Y); sABCD(C, C, R, R); sABCD(A, A, C, C)])
		//sol = [Pa^2*Pc^2; -Pc^2*(Pa + Pg)^2; -Pa^2*(Pc + Pt)^2; (Pa + Pg)^2*(Pc + Pt)^2]
		//lst2 = simplify([sABCD(Y, Y, R, R); sABCD(Y, Y, A, A); sABCD(R, R, C, C); sABCD(C, C, A, A)])

		// (Pa^2+Pg^2)*(Pc^2+Pt^2)*sABCD(R, R, Y, Y)
		// -Pr^2*(Pc^2+Pt^2)*sABCD(A, A, Y, Y) -Pr^2*(Pc^2+Pt^2)*sABCD(G, G, Y, Y) -(Pa^2+Pg^2)*Py^2*sABCD(C, C, R, R) -(Pa^2+Pg^2)*Py^2*sABCD(T, T, R, R)
		// +Pr^2*Py^2*sABCD(A, A, C, C) +Pr^2*Py^2*sABCD(G, G, C, C) +Pr^2*Py^2*sABCD(A, A, T, T) +Pr^2*Py^2*sABCD(G, G, T, T)

		score_t const A = pi[0], C = pi[1], G = pi[2], T = pi[3];
		score_t const R = A + G, Y = C + T, R2 = A * A + G * G, Y2 = C * C + T * T;
		cnt4_t const a0 = cnt0[0], c0 = cnt0[1], g0 = cnt0[2], t0 = cnt0[3], r0 = a0 + g0, y0 = c0 + t0;
		cnt4_t const a1 = cnt1[0], c1 = cnt1[1], g1 = cnt1[2], t1 = cnt1[3], r1 = a1 + g1, y1 = c1 + t1;
		cnt4_t const a2 = cnt2[0], c2 = cnt2[1], g2 = cnt2[2], t2 = cnt2[3], r2 = a2 + g2, y2 = c2 + t2;
		cnt4_t const a3 = cnt3[0], c3 = cnt3[1], g3 = cnt3[2], t3 = cnt3[3], r3 = a3 + g3, y3 = c3 + t3;

		cnt4_t const rryy = quadXXYY(r0, r1, r2, r3, y0, y1, y2, y3);

		cnt4_t const aayy = quadXXYY(a0, a1, a2, a3, y0, y1, y2, y3);
		cnt4_t const ggyy = quadXXYY(g0, g1, g2, g3, y0, y1, y2, y3);
		cnt4_t const rrcc = quadXXYY(r0, r1, r2, r3, c0, c1, c2, c3);
		cnt4_t const rrtt = quadXXYY(r0, r1, r2, r3, t0, t1, t2, t3);

		cnt4_t const aacc = quadXXYY(a0, a1, a2, a3, c0, c1, c2, c3);
		cnt4_t const aatt = quadXXYY(a0, a1, a2, a3, t0, t1, t2, t3);
		cnt4_t const ggcc = quadXXYY(g0, g1, g2, g3, c0, c1, c2, c3);
		cnt4_t const ggtt = quadXXYY(g0, g1, g2, g3, t0, t1, t2, t3);

		return rryy * R2 * Y2 - (aayy + ggyy) * (R * R) * Y2 - (rrcc + rrtt) * R2 * (Y * Y)
			+ (aacc + aatt + ggcc + ggtt) * (R * R) * (Y * Y);
	}

	inline static array<score_t, 3> quadPos(array<array<cnt_t, 4>, 4> const& cnt, array<score_t, 4> const& pi) noexcept {
		return { quadPos(cnt[0], cnt[1], cnt[2], cnt[3], pi),
				quadPos(cnt[0], cnt[2], cnt[1], cnt[3], pi),
				quadPos(cnt[0], cnt[3], cnt[1], cnt[2], pi) };
	}

public:
	void elementSetTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept{
		elementSetOrClearTaxonColor<true>(iElement, iTaxon, iColor);
	}
	
	void elementClearTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept{
		elementSetOrClearTaxonColor<false>(iElement, iTaxon, iColor);
	}
	
	// ══════════ WCASTER-MOD ⑧ BEGIN: read colorWeight/colorPairWeight ══════════
	score_t elementScore(size_t iElement) const noexcept{
		index_t iGenomePosBegin = sharedConstData.elements[iElement].iGenomePosBegin;
		index_t nPos = sharedConstData.elements[iElement].nPos;
		typename SharedConstData::Element const& element = sharedConstData.elements[iElement];

		score_t res = 0;
		for (index_t iPos : iota((index_t)0, nPos)){
			index_t gPos = iGenomePosBegin + iPos;
			res += scorePos(colorWeight[gPos], colorPairWeight[gPos], element.eqFreqs);
		}
		return res;
	}

	//Edge NNI score — quadPos uses weighted single sums only (no pair terms needed)
	array<score_t, 3> elementQuadripartitionScores(size_t iElement) const noexcept {
		index_t iGenomePosBegin = sharedConstData.elements[iElement].iGenomePosBegin;
		index_t nPos = sharedConstData.elements[iElement].nPos;
		typename SharedConstData::Element const& element = sharedConstData.elements[iElement];

		array<score_t, 3> res = {0, 0, 0};
		for (index_t iPos : iota((index_t)0, nPos)) {
			array<score_t, 3> part = quadPos(colorWeight[iGenomePosBegin + iPos], element.eqFreqs);
			for (index_t i : iota((index_t)0, (index_t)3)) res[i] += part[i];
		}
		return res;
	}
	// ══════════ WCASTER-MOD ⑧ END ══════════

	Color(SharedConstData const& data) noexcept : sharedConstData(data),
		colorWeight(data.nGenomePos), colorPairWeight(data.nGenomePos) {}

	template<typename DataClasses> friend DataClasses DriverHelper::read();
};

ChangeLog logDriverHelper("DriverHelper",
	"2026-02-04", "Chao Zhang", "Little code refactoring, no functional change", "patch",
	"2026-09-04", "Zuizhi Chen", "WCASTER: per-chunk weight computation in read()", "minor",
	"2026-11-13", "Zuizhi Chen", "Merge per-alignment weighting into the per-chunk tool as --alignment-weight", "minor");

namespace DriverHelper {

using namespace std;

template<typename T, typename T2> array<T, 4>& operator+=(array<T, 4>& a, const array<T2, 4>& b) {
	for (int j = 0; j < 4; j++) {
		a[j] += b[j];
	}
	return a;
}

template<typename T> T sum(const array<T, 4>& cnt) {
	T result = 0;
	for (int j = 0; j < 4; j++) {
		result += cnt[j];
	}
	return result;
}

// WCASTER: sequence-similarity weight of one species against the reference,
// measured on the half-open site range [iBegin, iEnd).
// Gap-containing sites are ignored. The weight is 0 when no site can be compared
// (this also keeps the per-chunk behaviour when a chunk is entirely gaps).
inline double sequenceSimilarityWeight(string const& refSeq, string const& seq, size_t iBegin, size_t iEnd) noexcept {
	size_t hamming = 0, nonGap = 0;
	for (size_t iPos = iBegin; iPos < iEnd; iPos++) {
		if (refSeq[iPos] == '-' || seq[iPos] == '-') continue;
		nonGap++;
		if (refSeq[iPos] != seq[iPos]) hamming++;
	}
	if (nonGap == 0) return 0.0;
	double const similarity = 1.0 - (double)hamming / (double)nonGap;
	return (similarity < 0.25) ? 0.0 : (similarity - 0.25) / 0.75;
}

// ══════════ WCASTER-MOD ⑨ BEGIN: fasta2ref input + per-chunk/per-alignment weight computation (照 wtrial) ══════════
template<typename DataClass> DataClass read() {
	using cnt_taxon_t = DataClass::ParentClass::cnt_taxon_t;
	using cnt_t = DataClass::ParentClass::cnt_t;

	common::LogInfo log(1);
	log.log() << "Parsing fasta2ref file and reading input..." << std::endl;
	DataClass sharedConstData;

	string fasta2refFile = ARG.get<string>("input");
	ifstream finF2R(fasta2refFile);
	if (!finF2R.is_open()) {
		throw std::logic_error("Cannot open fasta2ref file: " + fasta2refFile);
	}

	string fasta2refDir;
	{
		size_t pos = fasta2refFile.find_last_of("/\\");
		if (pos != string::npos) fasta2refDir = fasta2refFile.substr(0, pos + 1);
	}

	vector<string> fastaFiles;
	vector<size_t> fileRefTaxonIds;

	string fastaPath, refSpeciesName;
	while (finF2R >> fastaPath >> refSpeciesName) {
		if (fastaPath[0] != '/' && fastaPath[0] != '\\' && (fastaPath.size() < 2 || fastaPath[1] != ':')) {
			fastaPath = fasta2refDir + fastaPath;
		}
		fastaFiles.push_back(fastaPath);
		size_t refId = common::taxonName2ID[refSpeciesName];
		fileRefTaxonIds.push_back(refId);
	}
	finF2R.close();

	if (fastaFiles.empty()) {
		throw std::logic_error("No fasta files found in fasta2ref file: " + fasta2refFile);
	}

	log.log() << "#Fasta files: " << fastaFiles.size() << std::endl;

	//WCASTER: create temporary fasta list file for AlignmentParser
	string tempListFile = "/tmp/wcaster_fasta_list_" + std::to_string(std::hash<std::thread::id>()(std::this_thread::get_id())) + ".txt";
	{
		ofstream fout(tempListFile);
		for (const string& f : fastaFiles) fout << f << "\n";
	}

	aligment_utilities::AlignmentParser AP(tempListFile, 2), AP2(tempListFile, 3);
	size_t iFile = 0;

	// WCASTER: --alignment-weight selects one weight per alignment; otherwise weights are per chunk.
	bool const alignmentWeight = ARG.has("alignment-weight");
	size_t const chunkMaxSize = ARG.get<size_t>("chunk");

	while (AP.nextAlignment()) {
		size_t fileRefTaxonId = fileRefTaxonIds[iFile++];
		size_t nSites = AP.getLength();
		size_t nChunk = (nSites + chunkMaxSize - 1) / chunkMaxSize;
		vector<vector<size_t> > sites(nChunk);
		vector<array<double, 4> > eqfreq;
		size_t iElementBegin = sharedConstData.elements.size();
		unordered_map<size_t, size_t> taxon2row;

		string refSeq; // WCASTER: store ref sequence for weight computation (weight anchor)

		{
			size_t nTotalSpeciesmen = 0;
			unordered_map<size_t, size_t> nSpeciesmen;
			vector<array<unsigned short, 4> > freq;
			freq.resize(AP.getLength());
			while (AP.nextSeq()) {
				size_t iTaxon = common::taxonName2ID[AP.getName()];
				nTotalSpeciesmen++;
				nSpeciesmen[iTaxon]++;

				if (!taxon2row.count(iTaxon)) taxon2row[iTaxon] = taxon2row.size();
				string seq = AP.getSeq();
				if (iTaxon == fileRefTaxonId && refSeq.empty()) refSeq = seq; // WCASTER: capture ref
				for (size_t i = 0; i < seq.size(); i++) {
					switch (seq[i]) {
						case 'A': freq[i][0]++; break;
						case 'C': freq[i][1]++; break;
						case 'G': freq[i][2]++; break;
						case 'T': freq[i][3]++; break;
					}
				}
			}

			size_t maxSpeciesman = 0;
			for (auto const& element : nSpeciesmen) {
				maxSpeciesman = std::max(maxSpeciesman, element.second);
			}

			if (std::same_as<cnt_taxon_t, bool> && maxSpeciesman >= 2) {
				log.log() << "Seems there is more than one haploid genome per taxon and thus bool type cannot be used..." << std::endl;
				throw(std::logic_error("Incompatible data structure"));
			}
			if (std::same_as<cnt_taxon_t, unsigned char> && maxSpeciesman >= 256) {
				log.log() << "Seems there are more than 255 haploid genomes per taxon (which is fishy) and thus unsigned char type cannot be used..." << std::endl;
				throw(std::logic_error("Incompatible data structure"));
			}
			if (std::same_as<cnt_taxon_t, unsigned short> && maxSpeciesman >= 65536) {
				common::LogInfo err(-100);
				err.log() << "Seems there are more than 65535 haploid genomes per taxon (which is astonishing)! Please ask the author for a specially made version..." << std::endl;
				exit(-1);
			}

			for (size_t i = 0; i < nChunk; i++) {
					size_t s = i * nSites / nChunk, t = (i + 1) * nSites / nChunk;
					array<size_t, 4> sumFreq = {};
					for (size_t j = s; j < t; j++) {
						sumFreq += freq[j];
						#ifdef CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH
							sites[i].push_back(j);
						#else
							if (freq[j][0] + freq[j][2] >= 2 && freq[j][1] + freq[j][3] >= 2) sites[i].push_back(j);
						#endif
					}
					double total = sum(sumFreq);
					if (total > 0) eqfreq.push_back({ sumFreq[0] / total, sumFreq[1] / total, sumFreq[2] / total, sumFreq[3] / total });
					else eqfreq.push_back({ 0.25, 0.25, 0.25, 0.25 });
			}
		}

		AP2.nextAlignment();
		for (size_t i = 0; i < nChunk; i++) {
			typename DataClass::Element element;
			element.iGenomePosBegin = sharedConstData.nGenomePos;
			element.nPos = sites[i].size();
			element.cnts.resize(taxon2row.size(), vector<array<typename DataClass::ParentClass::cnt_taxon_t, 4> >(element.nPos));
			element.taxon2row.resize(common::taxonName2ID.nTaxa(), -1);
			element.eqFreqs = eqfreq[i];
			element.speciesWeights.resize(taxon2row.size(), 1.0); // WCASTER: init all weights to 1.0
			sharedConstData.elements.push_back(element);
			sharedConstData.nGenomePos += element.nPos;
		}

		while (AP2.nextSeq()) {
			size_t iTaxon = common::taxonName2ID[AP2.getName()];
			size_t iRow = taxon2row[iTaxon];
			string seq = AP2.getSeq();

			// WCASTER: weights are computed on all sites (including conserved ones), not only diverse sites.
			//   default: one weight per chunk (local alignment quality)
			//   --alignment-weight: one weight for the whole alignment, shared by all chunks
			if (iTaxon != fileRefTaxonId && !refSeq.empty()) {
				if (alignmentWeight) {
					double const w = sequenceSimilarityWeight(refSeq, seq, 0, nSites);
					for (size_t iChunk : iota((size_t) 0, nChunk))
						sharedConstData.elements[iElementBegin + iChunk].speciesWeights[iRow] = (cnt_t) w;
				} else {
					for (size_t iChunk : iota((size_t) 0, nChunk)) {
						size_t s = iChunk * nSites / nChunk, t = (iChunk + 1) * nSites / nChunk;
						double const w = sequenceSimilarityWeight(refSeq, seq, s, t);
						sharedConstData.elements[iElementBegin + iChunk].speciesWeights[iRow] = (cnt_t) w;
					}
				}
			}

			for (size_t iChunk : iota((size_t) 0, nChunk)) {
				typename DataClass::Element &element = sharedConstData.elements[iElementBegin + iChunk];
				element.taxon2row[iTaxon] = iRow;
				for (size_t iPos : iota((size_t) 0, sites[iChunk].size())) {
					switch (seq[sites[iChunk][iPos]]) {
						case 'A': element.cnts[iRow][iPos][0] += 1; break;
						case 'C': element.cnts[iRow][iPos][1] += 1; break;
						case 'G': element.cnts[iRow][iPos][2] += 1; break;
						case 'T': element.cnts[iRow][iPos][3] += 1; break;
					}
				}
			}
		}
	}

	std::remove(tempListFile.c_str());

	return sharedConstData;
}
// ══════════ WCASTER-MOD ⑨ END ══════════

};

ChangeLog logDriver("Driver",
	"2026-02-01", "Chao Zhang", "Change prgramName to caster", "patch",
	"2026-02-08", "Chao Zhang", "Adding more type support", "patch",
	"2026-09-04", "Zuizhi Chen", "WCASTER: add wcaster driver", "minor",
	"2026-11-13", "Zuizhi Chen", "Merge per-alignment weighting into the per-chunk tool as --alignment-weight", "minor");

template<bool> class Driver : public common::LogInfo
{
	using string = std::string;

public:
	// ══════════ WCASTER-MOD ⑪ BEGIN: 3 DataClass variants (cnt_t = double fixed) ══════════
	using DataClasses = std::variant<typename Color<StepwiseColorDefaultAttributes<bool> >::SharedConstData, typename Color<StepwiseColorDefaultAttributes<unsigned char> >::SharedConstData, typename Color<StepwiseColorDefaultAttributes<unsigned short> >::SharedConstData>;
	// ASTER2 branch annotations exposed by the driver as -u/--support and -l/--length:
	// quartet score (qs), local block bootstrap (lbb), SPR44 block bootstrap (sbb).
	using BranchAnnotations = branch_annotation::BranchAnnotations<quadripartition_support::BA_QS, quadripartition_support::BA_LBB, quadripartition_support_spr44::BA_SBB>;

	static std::string constexpr DEFAULT_SUPPORT_ANNOTATION() noexcept { return "sbb"; };
	static std::string constexpr DEFAULT_LENGTH_ANNOTATION() noexcept { return "none"; };
	
	static std::pair<string, string> programNames() noexcept {
		return { "wcaster", "Weighted CASTER: per-chunk (default) or per-alignment sequence-similarity weighting" };
	}

	static void addArguments() noexcept {
		ARG.addArgument('\0', "chunk", "integer", "The maximum number of sites in each local aligment block for parameter estimation", 0, true, true, "1000");
		ARG.addArgument('\0', "alignment-weight", "flag", "Use one sequence-similarity weight per alignment instead of one per chunk", 1, true);
	}

	static DataClasses getStepwiseColorSharedConstData() noexcept {
		try { return DriverHelper::read<std::variant_alternative_t<0, DataClasses> >(); } catch (...) {}
		try { return DriverHelper::read<std::variant_alternative_t<1, DataClasses> >(); } catch (...) {}
		return DriverHelper::read<std::variant_alternative_t<2, DataClasses> >();
	}
	// ══════════ WCASTER-MOD ⑪ END ══════════
};

class Documentation : public common::DocumentationBase {
protected:
	string introduction() const noexcept override {
		return R"YOHANETYO(# Weighted CASTER (WCASTER): CASTER with sequence-similarity weighting

WCASTER keeps CASTER's symmetric quartet scoring but weights sites by their sequence similarity to a reference species. Each alignment in the `fasta2ref` input has its own reference, which serves only as the weighting anchor: unlike TRIAL, quartets do not need to contain the reference.

## Weighting modes
* Per chunk (default): each alignment is split into chunks (elements), and each species receives one weight per chunk based on local sequence similarity. This down-weights poorly conserved or likely-misaligned regions at a fine scale.
* Per alignment (`--alignment-weight`): each species receives one weight for the whole alignment, shared by all of its chunks. This is the older per-alignment behaviour.

The weight is estimated from ungapped sites as `similarity = 1 - (#mismatch / #non-gap sites)` and mapped to `[0, 1]`:
```
weight = 0                          if similarity < 0.25
       = (similarity - 0.25) / 0.75 otherwise
```
The reference species itself has weight 1. The default chunk size is 1000 in both modes.

## Relationship to CASTER and TRIAL
CASTER scores all quartets without weighting. TRIAL restricts scoring to quartets containing the reference and prioritizes the reference during search. WCASTER keeps CASTER's scoring (the reference is only the weight anchor), so it reduces to CASTER when all weights are 1.

## Background: CASTER

[<img src="../misc/CASTER.png" width="500"/>](../misc/CASTER.png)

Genome-wide data have the promise of dramatically improving phylogenetic inferences. Yet, inferring the true phylogeny remains a challenge, mainly because the evolutionary histories of different genomic regions differ. The traditional concatenation approach ignores such differences, resulting in both theoretical and empirical shortcomings. In response, many discordance-aware inference methods have been developed. Yet, all have their own weaknesses. Many methods rely on short recombination-free genomic segments to build gene trees and thus suffer from a lack of signals for gene tree reconstruction, resulting in poor species tree. Some methods wrongly assume that the rate of evolution is uniform across the species tree. Yet, others lack enough scalability to analyze phylogenomic data.

We introduce a new site-based species tree inference method that seeks to address these challenges without reconstructing gene trees. Our method, called CASTER (Coalescence-aware Alignment-based Species Tree EstimatoR), has two flavors: CASTER-site and CASTER-pair. The first is based on patterns in individual sites and the second is based on pairs of sites.

CASTER has several outstanding features:
1. CASTER introduces two new optimization objectives based on genomic site patterns of four species; we show that optimizing these objectives produces two estimators: CASTER-site is statistically consistent under MSC+F84 model while allowing mutation rate to change across sites and across species tree branches; CASTER-pair is statistically consistent under MSC+LM model under further assumptions.
2. CASTER comes with a scalable algorithm to optimize the objectives summed over all species quartets. Remarkably, its time complexity is linear to the number of sites and at most quasi-quadratic with respect to the number of species.
3. CASTER can handle multiple samples per species, and CASTER-site specifically can work with allele frequencies of unphased multiploid.
4. CASTER is extremenly memory efficent, requiring <1 byte per SNP per sample

Under extensive simulation of genome-wide data, including recombination, we show that both CASTER-site and CASTER-pair out-perform concatenation using RAxML-ng, as well as discordance-aware methods SVDQuartets and wASTRAL in terms of both accuracy and running time. Noticeably, CASTER-site is 60�C150X fASTER2 than the alternative methods. It reconstructs an Avian tree of 51 species from aligned genomes with 254 million SNPs in only 3.5 hours on an 8-core desktop machine with 32 GB memory. It can also reconstruct a species tree of 201 species with approximately 2 billion SNPs using a server of 256 GB memory.

Our results suggest that CASTER-site and CASTER-pair can fulfill the need for large-scale phylogenomic inferences.

## Publication

Chao Zhang, Rasmus Nielsen, Siavash Mirarab, CASTER: Direct species tree inference from whole-genome alignments. Science (2025) https://www.science.org/doi/10.1126/science.adk9688

## Notice

Since CASTER-site and CASTER-pair assume different models, please run both and choose the result that makes more sense if you can.
)YOHANETYO";
	}

	string input() const noexcept override {
		return R"YOHANETYO(# INPUT

WCASTER uses a `fasta2ref` file to specify the mapping between fasta alignment files and their reference species. The reference is used only to compute sequence-similarity weights.

## fasta2ref file format
Each line contains a fasta file path followed by the reference species name, separated by whitespace:
```
/path/to/alignment1.fasta    Felis_catus
/path/to/alignment2.fasta    Otocolobus_manul
```

## Usage
```
bin/wcaster -i example/fasta2ref.txt -o output.tre
bin/wcaster -i example/fasta2ref.txt --alignment-weight -o output_alignment_weighted.tre
```
)YOHANETYO";
	}

	string programName() const noexcept override { return "wcaster"; }

	string exampleInput() const noexcept override { return "fasta2ref.txt"; }
};

};
#endif

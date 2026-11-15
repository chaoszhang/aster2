#ifndef TRIAL_HPP
#define TRIAL_HPP

#include "stepwise_colorable.hpp"
#include "quadripartition_support.hpp"
#include "quadripartition_support_spr44.hpp"
#include "alignment_utilities.hpp"

namespace trial{

using std::size_t;
using std::views::iota;
using std::vector;
using std::array;
using std::string;
using std::unordered_map;

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

template<typename cnt_taxon_type = unsigned char, typename cnt_type = unsigned short> struct StepwiseColorDefaultAttributes {
	using score_t = double;
	using cnt_taxon_t = cnt_taxon_type;
	using cnt_t = cnt_type;
	using cnt4_t = std::conditional_t<std::same_as<cnt_type, unsigned char>, unsigned int, unsigned long long>;
	using index_t = long long;
	static inline score_t constexpr ZERO = 0;
	static inline score_t constexpr EPSILON = 1e-3;
};

ChangeLog logColor("Color",
	"2026-02-02", "Chao Zhang", "Supporting quadripartiton", "minor",
	"2026-06-12", "Zuizhi Chen", "TRIAL: only score quartets containing reference species", "minor",
	"2026-06-15", "Zuizhi Chen", "TRIAL: adapt to TAXON_ORDER_PRIORITIZING interface", "minor",
	"2026-06-16", "Zuizhi Chen", "TRIAL: refactor scoring to direct calculation with refCnt/refColor", "minor",
	"2026-07-06", "Zuizhi Chen", "TRIAL: fix segfault when N_ref >= 5 (missing refColor guard in scoring)", "patch");

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
			size_t iReferenceTaxonId = (size_t)-1; //reference species this element was aligned to

			bool hasTaxon(size_t iTaxon) const noexcept{
				return iTaxon < taxon2row.size() && taxon2row[iTaxon] != -1;
			}
		};

		vector<Element> elements;
		index_t nGenomePos = 0;
		//all reference species taxon IDs used for priority ordering
		vector<size_t> priorityTaxa;


		size_t nElements() const noexcept { return elements.size(); }
    };

	//static method for TAXON_ORDER_PRIORITIZING concept
	static void taxonOrderPrioritizing(std::vector<size_t>& taxonOrder) noexcept {
		if (s_priorityTaxa.empty()) return;
		std::vector<size_t> priority, others;
		for (size_t t : taxonOrder) {
			if (std::find(s_priorityTaxa.begin(), s_priorityTaxa.end(), t) != s_priorityTaxa.end())
				priority.push_back(t);
			else
				others.push_back(t);
		}
		size_t i = 0;
		for (size_t t : priority) taxonOrder[i++] = t;
		for (size_t t : others) taxonOrder[i++] = t;
	}

	//static cache of priority taxa, set during read()
	static inline std::vector<size_t> s_priorityTaxa;

private:
	SharedConstData const& sharedConstData;
    vector<array<array<cnt_t, 4>, 4> > colorCnts; // colorCnts[iGenomePos][iColor][iNucleotide] -> count (non-ref only)
	vector<array<cnt_t, 4> > refCnt;     // refCnt[iGenomePos][iNucleotide] — reference species counts at each position
	vector<size_t> refColor;             // refColor[iGenomePos] — current color of ref at each position, (size_t)-1 if not set

	template<bool isSet> inline void elementSetOrClearTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept{
		typename SharedConstData::Element const& element = sharedConstData.elements[iElement];
		if (!element.hasTaxon(iTaxon)) return;
		index_t iRow = element.taxon2row[iTaxon];
		index_t iPosBegin = element.iGenomePosBegin;

		//ref species is excluded from colorCnts; only tracked via refColor
		if (iTaxon == element.iReferenceTaxonId) {
			for (index_t iPos : iota((index_t)0, element.nPos))
				refColor[iPosBegin + iPos] = (isSet) ? iColor : (size_t)-1;
			return;
		}

		for (index_t iPos : iota((index_t)0, element.nPos)){
			for (index_t iNucleotide : iota((index_t)0, (index_t)4)) {
				cnt_t& colorCnt = colorCnts[iPosBegin + iPos][iColor][iNucleotide];
				cnt_t cnt = element.cnts[iRow][iPos][iNucleotide];
				if constexpr (isSet) colorCnt += cnt;
				else colorCnt -= cnt;
			}
		}
	}
	
	inline static cnt4_t XXYY(cnt4_t xR, cnt4_t x0, cnt4_t x1, cnt4_t x2, cnt4_t yR, cnt4_t y0, cnt4_t y1, cnt4_t y2) noexcept{
		return xR * x0 * y1 * y2 * 2 + yR * y0 * x1 * x2 * 2 +
		       xR * (xR - 1) * y1 * y2 + yR * (yR - 1) * x1 * x2 +
		       xR * x1 * y2 * (y2 - 1) + yR * y1 * x2 * (x2 - 1) +
		       xR * x2 * y1 * (y1 - 1) + yR * y2 * x1 * (x1 - 1);
	}

	//scorePos with ref in virtual color 0
	inline static score_t scorePos(array<array<cnt_t, 4>, 4> const &cnt, array<cnt_t, 4> const &rCnt, array<score_t, 4> const &pi) noexcept{
		cnt4_t const aR = rCnt[0], cR = rCnt[1], gR = rCnt[2], tR = rCnt[3];
		cnt4_t const a0 = cnt[0][0], c0 = cnt[0][1], g0 = cnt[0][2], t0 = cnt[0][3];
		cnt4_t const a1 = cnt[1][0], c1 = cnt[1][1], g1 = cnt[1][2], t1 = cnt[1][3];
		cnt4_t const a2 = cnt[2][0], c2 = cnt[2][1], g2 = cnt[2][2], t2 = cnt[2][3];

		score_t const A = pi[0], C = pi[1], G = pi[2], T = pi[3];
		score_t const R = A + G, Y = C + T, R2 = A * A + G * G, Y2 = C * C + T * T;
		cnt4_t const r0 = a0 + g0, y0 = c0 + t0;
		cnt4_t const r1 = a1 + g1, y1 = c1 + t1;
		cnt4_t const r2 = a2 + g2, y2 = c2 + t2;
		cnt4_t const rR = aR + gR, yR = cR + tR;

		cnt4_t const rryy = XXYY(rR, r0, r1, r2, yR, y0, y1, y2);
		cnt4_t const aayy = XXYY(aR, a0, a1, a2, yR, y0, y1, y2);
		cnt4_t const ggyy = XXYY(gR, g0, g1, g2, yR, y0, y1, y2);
		cnt4_t const rrcc = XXYY(rR, r0, r1, r2, cR, c0, c1, c2);
		cnt4_t const rrtt = XXYY(rR, r0, r1, r2, tR, t0, t1, t2);
		cnt4_t const aacc = XXYY(aR, a0, a1, a2, cR, c0, c1, c2);
		cnt4_t const aatt = XXYY(aR, a0, a1, a2, tR, t0, t1, t2);
		cnt4_t const ggcc = XXYY(gR, g0, g1, g2, cR, c0, c1, c2);
		cnt4_t const ggtt = XXYY(gR, g0, g1, g2, tR, t0, t1, t2);

		return rryy * R2 * Y2 - (aayy + ggyy) * (R * R) * Y2 - (rrcc + rrtt) * R2 * (Y * Y)
		     + (aacc + aatt + ggcc + ggtt) * (R * R) * (Y * Y);
	}
	
	inline static cnt4_t quadXXYY(cnt4_t x0, cnt4_t x1, cnt4_t x2, cnt4_t x3, cnt4_t y0, cnt4_t y1, cnt4_t y2, cnt4_t y3) noexcept {
		return x0 * x1 * y2 * y3 + y0 * y1 * x2 * x3;
	}

	static score_t quadPosSingle(array<cnt_t, 4> const& cnt0, array<cnt_t, 4> const& cnt1,
		array<cnt_t, 4> const& cnt2, array<cnt_t, 4> const& cnt3, array<score_t, 4> const& pi) noexcept {

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
		return {quadPosSingle(cnt[0], cnt[1], cnt[2], cnt[3], pi),
		        quadPosSingle(cnt[0], cnt[2], cnt[1], cnt[3], pi),
		        quadPosSingle(cnt[0], cnt[3], cnt[1], cnt[2], pi)};
	}

public:
	void elementSetTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept{
		elementSetOrClearTaxonColor<true>(iElement, iTaxon, iColor);
	}
	
	void elementClearTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept{
		elementSetOrClearTaxonColor<false>(iElement, iTaxon, iColor);
	}
	
	//Tripartition score: scorePos with ref in virtual color 0;
	score_t elementScore(size_t iElement) const noexcept{
		if (refColor[sharedConstData.elements[iElement].iGenomePosBegin] == (size_t)-1) return 0;
		index_t iGenomePosBegin = sharedConstData.elements[iElement].iGenomePosBegin;
		index_t nPos = sharedConstData.elements[iElement].nPos;
		typename SharedConstData::Element const& element = sharedConstData.elements[iElement];

		score_t res = 0;
		for (index_t iPos : iota((index_t)0, nPos)){
			index_t gPos = iGenomePosBegin + iPos;
			size_t C = refColor[gPos];
			if (C == (size_t)-1) continue;
			array<array<cnt_t, 4>, 4> c = colorCnts[gPos];
			if (C != 0) std::swap(c[0], c[C]);
			res += scorePos(c, refCnt[gPos], element.eqFreqs);
		}
		return res;
	}

	//Edge NNI score
	array<score_t, 3> elementQuadripartitionScores(size_t iElement) const noexcept {
		if (refColor[sharedConstData.elements[iElement].iGenomePosBegin] == (size_t)-1) return {0, 0, 0};
		index_t iGenomePosBegin = sharedConstData.elements[iElement].iGenomePosBegin;
		index_t nPos = sharedConstData.elements[iElement].nPos;
		typename SharedConstData::Element const& element = sharedConstData.elements[iElement];

		array<score_t, 3> res = {0, 0, 0};
		for (index_t iPos : iota((index_t)0, nPos)) {
			index_t gPos = iGenomePosBegin + iPos;
			size_t C = refColor[gPos];
			if (C == (size_t)-1) continue;
			array<array<cnt_t, 4>, 4> c = colorCnts[gPos];
			c[C] = refCnt[gPos];
			array<score_t, 3> part = quadPos(c, element.eqFreqs);
			for (index_t i : iota((index_t)0, (index_t)3)) res[i] += part[i];
		}
		return res;
	}

	//TRIAL: constructor precomputes refCnt from elements;
	Color(SharedConstData const& data) noexcept : sharedConstData(data), colorCnts(data.nGenomePos), refCnt(data.nGenomePos), refColor(data.nGenomePos, (size_t)-1) {
		for (auto const& element : sharedConstData.elements) {
			if (!element.hasTaxon(element.iReferenceTaxonId)) continue;
			index_t refRow = element.taxon2row[element.iReferenceTaxonId];
			for (index_t iPos : iota((index_t)0, element.nPos)) {
				index_t gPos = element.iGenomePosBegin + iPos;
				for (index_t nuc : iota((index_t)0, (index_t)4))
					refCnt[gPos][nuc] = element.cnts[refRow][iPos][nuc];
			}
		}
	}

	template<typename DataClasses> friend DataClasses DriverHelper::read();
};

ChangeLog logDriverHelper("DriverHelper",
	"2026-02-04", "Chao Zhang", "Little code refactoring, no functional change", "patch",
	"2026-06-12", "Zuizhi Chen", "TRIAL: add fasta2ref parsing, multi-fasta input, and multi-ref support", "minor");

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

//TRIAL: fasta2ref parsing, multi-fasta input, multi-ref support
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
	vector<size_t>& priorityTaxa = sharedConstData.priorityTaxa;
	unordered_set<size_t> seenRefs;

	string fastaPath, refSpeciesName;
	while (finF2R >> fastaPath >> refSpeciesName) {
		if (fastaPath[0] != '/' && fastaPath[0] != '\\' && (fastaPath.size() < 2 || fastaPath[1] != ':')) {
			fastaPath = fasta2refDir + fastaPath;
		}
		fastaFiles.push_back(fastaPath);
		size_t refId = common::taxonName2ID[refSpeciesName];
		fileRefTaxonIds.push_back(refId);
		if (!seenRefs.count(refId)) {
			seenRefs.insert(refId);
			priorityTaxa.push_back(refId);
		}
	}
	finF2R.close();

	if (fastaFiles.empty()) {
		throw std::logic_error("No fasta files found in fasta2ref file: " + fasta2refFile);
	}

	log.log() << "#Reference species (distinct): " << priorityTaxa.size() << std::endl;
	for (size_t rId : priorityTaxa) log.log() << "  " << common::taxonName2ID[rId] << std::endl;
	log.log() << "#Fasta files: " << fastaFiles.size() << std::endl;

	//TRIAL: create temporary fasta list file for AlignmentParser
	string tempListFile = "/tmp/trial_fasta_list_" + std::to_string(std::hash<std::thread::id>()(std::this_thread::get_id())) + ".txt";
	{
		ofstream fout(tempListFile);
		for (const string& f : fastaFiles) fout << f << "\n";
	}

	aligment_utilities::AlignmentParser AP(tempListFile, 2), AP2(tempListFile, 3);
	size_t iFile = 0;

	while (AP.nextAlignment()) {
		size_t fileRefTaxonId = fileRefTaxonIds[iFile++];
		size_t nSites = AP.getLength();
		size_t chunkMaxSize = ARG.get<size_t>("chunk");
		size_t nChunk = (nSites + chunkMaxSize - 1) / chunkMaxSize;
		vector<vector<size_t> > sites(nChunk);
		vector<array<double, 4> > eqfreq;
		size_t iElementBegin = sharedConstData.elements.size();
		unordered_map<size_t, size_t> taxon2row;
			
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
			for (auto const& elem : nSpeciesmen) {
				maxSpeciesman = std::max(maxSpeciesman, elem.second);
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
			if (std::same_as<cnt_t, unsigned char> && nTotalSpeciesmen >= 256) {
				log.log() << "Seems there are more than 255 haploid genomes in total and thus unsigned char type cannot be used..." << std::endl;
				throw(std::logic_error("Incompatible data structure"));
			}
			if (std::same_as<cnt_t, unsigned short> && nTotalSpeciesmen >= 65536) {
				common::LogInfo err(-100);
				err.log() << "Seems there are more than 65535 haploid genomes in total (which is fishy)! Please ask the author for a specially made version..." << std::endl;
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
			element.iReferenceTaxonId = fileRefTaxonId; //tag element with its reference species
			sharedConstData.elements.push_back(element);
			sharedConstData.nGenomePos += element.nPos;
		}

		while (AP2.nextSeq()) {
			size_t iTaxon = common::taxonName2ID[AP2.getName()];
			size_t iRow = taxon2row[iTaxon];
			string seq = AP2.getSeq();
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

	std::remove(tempListFile.c_str()); //clean up temp file

	DataClass::ParentClass::s_priorityTaxa = sharedConstData.priorityTaxa;

	return sharedConstData;
}

};

ChangeLog logDriver("Driver",
	"2026-02-01", "Chao Zhang", "Change prgramName to caster", "patch",
	"2026-02-08", "Chao Zhang", "Adding more type support", "patch",
	"2026-06-12", "Zuizhi Chen", "TRIAL: use -i for fasta2ref input", "minor");

template<bool> class Driver : public common::LogInfo
{
	using string = std::string;

public:
	using DataClasses = std::variant<typename Color<StepwiseColorDefaultAttributes<bool, unsigned char> >::SharedConstData, typename Color<StepwiseColorDefaultAttributes<unsigned char, unsigned char> >::SharedConstData, typename Color<StepwiseColorDefaultAttributes<bool, unsigned short> >::SharedConstData, typename Color<StepwiseColorDefaultAttributes<unsigned char, unsigned short> >::SharedConstData, typename Color<StepwiseColorDefaultAttributes<unsigned short, unsigned short> >::SharedConstData>;
	// ASTER2 branch annotations exposed by the driver as -u/--support and -l/--length:
	// quartet score (qs), local block bootstrap (lbb), SPR44 block bootstrap (sbb).
	using BranchAnnotations = branch_annotation::BranchAnnotations<quadripartition_support::BA_QS, quadripartition_support::BA_LBB, quadripartition_support_spr44::BA_SBB>;

	static std::string constexpr DEFAULT_SUPPORT_ANNOTATION() noexcept { return "sbb"; };
	static std::string constexpr DEFAULT_LENGTH_ANNOTATION() noexcept { return "none"; };

	static std::pair<string, string> programNames() noexcept {
		return { "trial", "TRIAL: Triplet Algorithm-based CASTER" };
	}

	static void addArguments() noexcept {
		ARG.addArgument('\0', "chunk", "integer", "The maximum number of sites in each local aligment block for parameter estimation", 0, true, true, "10000");
	}

	static DataClasses getStepwiseColorSharedConstData() noexcept {
		try { return DriverHelper::read<std::variant_alternative_t<0, DataClasses> >(); } catch (...) {}
		try { return DriverHelper::read<std::variant_alternative_t<1, DataClasses> >(); } catch (...) {}
		try { return DriverHelper::read<std::variant_alternative_t<2, DataClasses> >(); } catch (...) {}
		try { return DriverHelper::read<std::variant_alternative_t<3, DataClasses> >(); } catch (...) {}
		return DriverHelper::read<std::variant_alternative_t<4, DataClasses> >();
	}
};

class Documentation : public common::DocumentationBase {
protected:
	string introduction() const noexcept override {
		return R"YOHANETYO(# TRIAL: Triplet Algorithm-based CASTER

TRIAL (Triplet Algorithm-based CASTER) is a CASTER variant designed to reduce the effect of homology errors in whole-genome alignments. In a multiple sequence alignment, sequences that are more closely related to the reference genome are generally aligned more accurately, while more distant sequences accumulate misalignments (for example, paralogous regions caused by gene duplication and loss). TRIAL therefore scores only quartets that contain a designated reference species (R), reducing the number of potentially misaligned non-reference sequences in each scored quartet from four to three. Multiple reference genomes can be combined to cover different sets of close relatives and to reuse existing public reference assemblies.

## How it works
1. Multiple fasta alignment files are provided via a `fasta2ref` file; each alignment has its own reference species.
2. Only quartets containing the reference species contribute to the tree score.
3. The reference species is placed first during heuristic search.
4. Multiple references provide cross-coverage: different references can cover different close relatives, and no single reference has to cover all taxa.

## Why reference anchoring reduces noise
CASTER scores all quartets in the alignment. When the alignment contains homology errors, every quartet can be affected. TRIAL excludes quartets without the reference, so each scored quartet contains at most three non-reference sequences that may carry alignment errors. Using several references spreads this risk across different sets of close relatives and makes use of existing public reference genomes.

The underlying quartet model and statistical consistency follow CASTER. See the CASTER publication for details.

)YOHANETYO";
	}

	string input() const noexcept override {
		return R"YOHANETYO(# INPUT

TRIAL uses a `fasta2ref` file to specify the mapping between fasta alignment files and their reference species.

## fasta2ref file format
Each line contains a fasta file path followed by the reference species name, separated by whitespace:
```
/path/to/alignment1.fasta    Felis_catus
/path/to/alignment2.fasta    Otocolobus_manul
```

## Usage
```
bin/trial -i example/fasta2ref.txt -o output.tre
```
)YOHANETYO";
	}

	string programName() const noexcept override { return "trial"; }

	string exampleInput() const noexcept override { return "fasta2ref.txt"; }
};

};
#endif

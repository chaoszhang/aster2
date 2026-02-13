#include "common.hpp"

#ifndef ALIGNMENT_UTILITIES
#define ALIGNMENT_UTILITIES

namespace aligment_utilities {

using namespace std;

ChangeLog logAlignmentParser("AlignmentParser",
    "2026-02-13", "Chao Zhang", "Speeding up preprocess", "patch");

class AlignmentParser : common::LogInfo {
    bool isFasta = false, isFastaList = false, isPhylip = false;
    bool formatAmbiguity = false, formatNA = false;
    bool ambiguitySecond = true;

    queue<string> fastaQueue;
    string bufferName, bufferSeq, nextLine;
    size_t length;
    int phylipNspecies;
    bool firstFastaSeq;
    ifstream fin;
    vector<char> file_buffer;
    char NA[256] = {};

    static size_t constexpr BUFFER_SIZE = 1024 * 1024;

public:
    AlignmentParser(string fileName, int verbose = common::LogInfo::DEFAULT_VERBOSE, string fileFormat = "auto", string seqFormat = "NA"): LogInfo(verbose), file_buffer(BUFFER_SIZE) {
        fin.rdbuf()->pubsetbuf(file_buffer.data(), BUFFER_SIZE);

        ifstream ftemp(fileName);
        string temp;
        if (!getline(ftemp, temp) || temp.length() == 0) {
            cerr << "Input file path '" << fileName << "' seems empty!\n";
			throw invalid_argument("Empty input file");
        }
        if (fileFormat == "phylip" || (fileFormat == "auto" && seemsPhylip(temp))) initPhylip(fileName);
        else if (fileFormat == "fasta" || (fileFormat == "auto" && seemsFasta(temp))) initFasta(fileName);
        else initFastaList(fileName);
        if (seqFormat == "ambiguity") formatAmbiguity = true;
        else formatNA = true;

        for (size_t i : std::views::iota(0, 256)) NA[i] = '-';
        NA['A'] = NA['a'] = 'A';
        NA['C'] = NA['c'] = 'C';
        NA['G'] = NA['g'] = 'G';
        NA['T'] = NA['t'] = NA['U'] = NA['u'] = 'T';
    }

    bool nextAlignment() {
        ambiguitySecond = true;
        if (isFasta) return nextAlignmentFasta();
        else if (isPhylip) return nextAlignmentPhylip();
        else return false;
    }

    size_t getLength() {
        return length;
    }

    bool nextSeq() {
        if (formatAmbiguity) {
            ambiguitySecond = (!ambiguitySecond);
            if (ambiguitySecond) return true;
        }

        if (isFasta) return nextSeqFasta();
        else if (isPhylip) return nextSeqPhylip();
        else return false;
    }

    string getName(){
		size_t len = bufferName.length();
        if (len > 0 && bufferName[len - 1] == '\r') {
            bufferName = bufferName.substr(0, len - 1);
            LogInfo vlog(verbose + 1);
			vlog << "Warning: Carriage return (\\r) detected in sequence name '" << bufferName << "'. Removed, but may case bugs elsewhere.\n";
        }
        return bufferName;
    }

    string getSeq() {
		size_t len = bufferSeq.length();
        if (len > 0 && bufferSeq[len - 1] == '\r') {
            bufferSeq = bufferSeq.substr(0, len - 1);
			log() << "Warning: Carriage return (\\r) detected in sequence of '" << bufferName << "'. Removed, but may case bugs elsewhere.\n";
        }
        if (bufferSeq.length() != length) {
            cerr << "Error: Sequence length of " << bufferName << " is not consistent!\n";
			throw invalid_argument("Inconsistent sequence length");
        }
        if (formatNA) return toFormatNA(bufferSeq);
        if (formatAmbiguity) return toAmbiguity(bufferSeq);
        return bufferSeq;
    }

private:
    static string getFastaName(const string& fasta) {
        size_t i = 0;
        string res;
        if (i == fasta.size() || fasta[i] != '>') {
            cerr << "Error in parsing line '" << fasta << "' in FASTA format.\n";
			throw invalid_argument("FASTA format error");
        }
        i++;
        while (i < fasta.size() && fasta[i] == ' ') i++;
        while (i < fasta.size() && fasta[i] != ' ') {
            res += fasta[i];
            i++;
        }
        if (res.size() == 0) {
            cerr << "Error in parsing line '" << fasta << "' in FASTA format.\n";
            throw invalid_argument("FASTA format error");
        }
        return res;
    }

    static string removeSpace(string& seq) {
        string res = std::move(seq);
        size_t i = 0;
        for (char c : res) {
            if (c != ' ' && c != '\t' && c != '\r') res[i++] = c;
        }
        res.resize(i);
        return res;
    }

    static bool seemsPhylip(const string& line) {
        size_t i = 0;
        while (i < line.length() && (line[i] == ' ' || line[i] == '\t')) i++;
        if (i == line.length() || line[i] < '0' || line[i] > '9') return false;
        while (i < line.length() && line[i] >= '0' && line[i] <= '9') i++;
        if (i == line.length() || (line[i] != ' ' && line[i] != '\t')) return false;
        while (i < line.length() && (line[i] == ' ' || line[i] == '\t')) i++;
        if (i == line.length() || line[i] < '0' || line[i] > '9') return false;
        while (i < line.length() && line[i] >= '0' && line[i] <= '9') i++;
        while (i < line.length() && (line[i] == ' ' || line[i] == '\t')) i++;
        return i == line.length();
    }

    static bool seemsFasta(const string& line) {
        return line[0] == '>';
    }

    void initPhylip(const string& fileName) {
        isPhylip = true;
        fin.open(fileName);
    }

    void initFasta(const string& fileName) {
        isFasta = true;
        fastaQueue.push(fileName);
    }

    void initFastaList(const string& fileName) {
        isFasta = true;
        ifstream ftemp(fileName);
        string temp;
        while (getline(ftemp, temp)) {
            fastaQueue.push(temp);
        }
    }

    bool parseSeqFasta() {
        bufferName = "";
        bufferSeq = "";
        string line;
        if (nextLine.length()) {
            line = nextLine;
            nextLine = "";
        }
        else if (!getline(fin, line)) return false;

        bufferName = getFastaName(line);
        while (getline(fin, line)) {
            if (line[0] == '>') {
                nextLine = line;
                break;
            }
            bufferSeq += removeSpace(line);
        }

        return true;
    }

    bool nextAlignmentFasta() {
        if (fastaQueue.empty()) return false;
        string fileName = fastaQueue.front();
        fastaQueue.pop();
        log() << "Processing " << fileName << " ... \n";
        fin.close();
        fin.open(fileName);
        if (!parseSeqFasta()) {
            cerr << "Error: FASTA file empty!\n";
			throw invalid_argument("Empty FASTA file");
        }
        length = bufferSeq.size();
        firstFastaSeq = true;
        return true;
    }

    bool nextSeqFasta() {
        if (firstFastaSeq) {
            firstFastaSeq = false;
            return true;
        }
        return parseSeqFasta();
    }

    bool nextAlignmentPhylip() {
        if (!(fin >> phylipNspecies)) return false;
        fin >> length;
        log() << "Reading alignment of " << phylipNspecies << " species of length " << length << " ...\n";
        return true;
    }

    bool nextSeqPhylip() {
        if (phylipNspecies == 0) return false;
        phylipNspecies--;
        fin >> bufferName >> bufferSeq;
        bufferSeq = removeSpace(bufferSeq);
        return true;
    }

    string toFormatNA(string& seq) {  
        string res = std::move(seq);
        for (char &c : res) c = NA[(unsigned char) c];
        return res;
    }

    string toAmbiguity(const string& seq) {
        string res;
        if (!ambiguitySecond) {
            for (char c : seq) {
                switch (c) {
                case 'A': case 'a': case 'M': case 'm': case 'R': case 'r': case 'W': case 'w': res += 'A'; break;
                case 'C': case 'c': case 'S': case 's': case 'Y': case 'y': res += 'C'; break;
                case 'G': case 'g': case 'K': case 'k': res += 'G'; break;
                case 'T': case 't': case 'U': case 'u': res += 'T'; break;
                default: res += '-';
                }
            }
        }
        else {
            for (char c : seq) {
                switch (c) {
                case 'A': case 'a': res += 'A'; break;
                case 'C': case 'c': case 'M': case 'm': res += 'C'; break;
                case 'G': case 'g': case 'R': case 'r': case 'S': case 's': res += 'G'; break;
                case 'T': case 't': case 'U': case 'u': case 'W': case 'w': case 'Y': case 'y': case 'K': case 'k': res += 'T'; break;
                default: res += '-';
                }
            }
        }
        return res;
    }
};

};

#endif

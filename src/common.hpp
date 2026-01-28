#ifndef COMMON_HPP
#define COMMON_HPP

#include "std_include.hpp"

namespace common{

using std::size_t;
using std::string;
using std::any;
using std::vector;
using std::unordered_map;
using std::unique_ptr;
using namespace std::string_literals;

class ChangeLog {
	struct Change {
		string code;
		string date;
		string author;
		string description;
		string level;

		Change(string const& code, string const& date, string const& author, string const& description, string const& level) noexcept : code(code), date(date), author(author), description(description), level(level) {}

		bool operator<(const Change& other) const noexcept {
			return date < other.date;
		}
	};

	struct Shared {
		vector<Change> all;
		unordered_map<string, vector<Change> > byCode;
		string globalVersion;
	};

public:

	static Shared shared;

	template<typename... Args> ChangeLog(string const& code, string const& date, string const& author, string const& description, string const& level, Args... args) {
		if (date.size() != 10) {
			std::cerr << "Error: Invalid date format for change log entry (YYYY-MM-DD): " << date << std::endl;
			throw;
		}
		if (level != "minor" && level != "patch") {
			std::cerr << "Error: Invalid level format for change log entry (minor or patch): " << level << std::endl;
			throw;
		}

		Change change(code, date, author, description, level);
		shared.all.push_back(change);
		shared.byCode[code].push_back(change);
		if constexpr (sizeof...(args) > 0) {
			ChangeLog(code, std::forward<Args>(args)...);
		}
	}

	template<bool isGlobal> static vector<string> versioning(vector<Change>& changes) {
		string PREFIX = (isGlobal) ? "v2."s : "v"s;
		std::sort(changes.begin(), changes.end());
		string lastDate = "1994-12-07"s, updates, version;
		bool hasMinor = false, hasPatch = false;
		vector<string> blocks;
		size_t minorVersion = 0, patchVersion = 0;
		for (const Change& change : changes) {
			if (change.date > lastDate) {
				blocks.push_back(version + updates);
				version = ""s;
				updates = " - "s + change.date + ":\n"s;
				lastDate = change.date;
				hasMinor = false;
				hasPatch = false;
			}
			if (change.level == "minor"s) {
				if (!hasMinor) minorVersion++;
				patchVersion = 0;
				hasMinor = true;
				hasPatch = true;
			}
			else if (change.level == "patch"s) {
				if (!hasPatch) patchVersion++;
				hasPatch = true;
			}
			version = PREFIX + std::to_string(minorVersion) + "."s + std::to_string(patchVersion);
			updates += "  * "s;
			if (isGlobal) updates += change.code + ": "s;
			updates += change.author + " - "s + change.description + "\n"s;
		}
		blocks.push_back(version + updates);
		if (isGlobal) shared.globalVersion = (blocks.size() > 1) ? version : "v2.0.0";
		return blocks;
	}

	template<bool isGlobal> static void displayAll(std::ostream& out, vector<Change>& changes) {
		for (const string& block : versioning<isGlobal>(changes) | std::views::reverse) out << block;
	}

	string static getGlobalVersion() {
		versioning<true>(shared.all);
		return shared.globalVersion;
	}
};
ChangeLog::Shared ChangeLog::shared;
/*
ChangeLog log("ChangeLog",
"2026-01-25", "Chao Zhang", "New function", "minor",
"2026-01-26", "Chao Zhang", "Initial patch", "patch");
*/

class LogInfo {
	struct Shared {
		unique_ptr<std::ostream> fout;
		int logVerbose = 0, errVerbose = 0;
	};

	static Shared shared;

public:
	static inline int constexpr DEFAULT_VERBOSE = 99;

	LogInfo const& log = *this;
	int verbose;

	LogInfo(int verbose = DEFAULT_VERBOSE) : verbose(verbose) {}

	template<int verbose_delta> LogInfo vlog() const {
		return LogInfo(verbose + verbose_delta);
	}

	LogInfo const& operator<<(auto const& t) const requires requires { {std::cerr << t}; } {
		if (verbose <= shared.errVerbose) std::cerr << t;
		if (shared.fout && verbose <= shared.logVerbose) (*shared.fout) << t;
		return log;
	}

	LogInfo const& operator<<(std::ostream&t(std::ostream&)) const {
		if (verbose <= shared.errVerbose) std::cerr << t;
		if (shared.fout && verbose <= shared.logVerbose) (*shared.fout) << t;
		return log;
	}

	static void setVerbose(std::ostream* fout, int logVerbose, int errVerbose) {
		shared.errVerbose = errVerbose;
		shared.logVerbose = logVerbose;
		if (fout) shared.fout.reset(fout);
	}
};
LogInfo::Shared LogInfo::shared;

template<typename T> concept LOG_OR_OSTREAM = std::same_as<T, const LogInfo> || std::derived_from<T, std::ostream>;

class Attributes: public LogInfo{
public:
	template<typename T> using String = string;
	template<typename T> using StringVector = vector<string>;

private:
	template<typename T, typename... Args> std::ostream &displayAttributesHelper(std::ostream &out, string const &attr, String<Args> const &... args) const{
		if (has(attr)) out << ";" << attr << "=" << std::any_cast<T>(get(attr));
		if constexpr (sizeof...(Args) > 0) displayAttributesHelper<Args...>(out, args...);
		return out;
	}
	
	template<typename T, typename... Args> std::ostream &displayAttributesHelper(std::ostream &out, vector<string> const &attrs, StringVector<Args> const &... args) const{
		for (const string &attr: attrs){
			if (has(attr)) out << ";" << attr << "=" << std::any_cast<T>(get(attr));
		}
		if constexpr (sizeof...(Args) > 0) displayAttributesHelper<Args...>(out, args...);
		return out;
	}
	
	unordered_map<string, any> attributes;
	
public:
	bool has(string const &attr) const noexcept{
		return attr.size() > 0 && attributes.contains(attr);
	}
	
	template<typename T> bool has(string const &attr) const noexcept{
		return has(attr) && attributes.at(attr).type() == typeid(T);
	}
	
	template<typename T> void set(string const &attr, T const &value) noexcept{
		attributes[attr] = value;
	}
	
	any get(string const &attr) const{
		return attributes.at(attr);
	}
	
	template<typename T> T get(string const &attr) const{
		return std::any_cast<T>(get(attr));
	}
	
	any& operator[](string const &attr) noexcept{
		return attributes[attr];
	}
	
	void erase(string const &attr) noexcept{
		attributes.erase(attr);
	}
	
	unordered_map<string, any> const &getAll() const noexcept{
		return attributes;
	}
	
	Attributes() noexcept{}
	
	template<std::ranges::range Range> requires std::convertible_to<std::ranges::range_value_t<Range>, string> Attributes(Attributes const &other, Range const &range) noexcept{
		for (string const &attr: range){
			if (other.has(attr)) attributes[attr] = other.get(attr);
		}
	}
	
	template<std::ranges::range Range> requires std::convertible_to<std::ranges::range_value_t<Range>, std::pair<string, any> > Attributes(Range const &range) noexcept{
		for (std::pair<string, any> const &attr: range){
			attributes[attr.first] = attr.second;
		}
	}
	
	template<typename S, typename T, typename... Args> requires std::convertible_to<S, string> Attributes(S attr, T value, Args... args) noexcept: Attributes(std::forward<Args>(args)...){
		set(attr, value);
	}
	
	template<typename T, typename... Args> std::ostream &displayAttributes(std::ostream &out, string const &attr, String<Args> const &... args) const{
		if (has(attr)) {
			out << attr << "=" << std::any_cast<T>(get(attr));
			if constexpr (sizeof...(Args) > 0) displayAttributesHelper<Args...>(out, args...);
		}
		else if constexpr (sizeof...(Args) > 0) displayAttributes<Args...>(out, args...);
		return out;
	}
	
	template<typename T, typename... Args> std::ostream &displayAttributes(std::ostream &out, vector<string> const &attrs, StringVector<Args> const &... args) const{
		bool first = true;
		for (const string &attr: attrs){
			if (has(attr)) {
				if (first) first = false;
				else out << ";";
				out << attr << "=" << std::any_cast<T>(get(attr));
			}
		}
		if constexpr (sizeof...(Args) > 0) {
			if (first) displayAttributes<Args...>(out, args...);
			else displayAttributesHelper<Args...>(out, args...);
		}
		return out;
	}
};

class InputParser : public Attributes {
	struct Argument {
		char shortcut;
		string name, type, description, defaultValue;
		bool optional, hasDefaultValue;
		int priority;
		std::shared_ptr<vector<string> > defaultEquivalence;

		bool operator<(const Argument& other) const noexcept {
			if (shortcut != 0 && other.shortcut == 0) return true;
			if (shortcut == 0 && other.shortcut != 0) return false;
			if (priority != other.priority) return priority > other.priority;
			if (optional != other.optional) return !optional;
			return name < other.name;
		}

		Argument(){}

		Argument(char shortcut, string const& name, string const& type, string const& description, int priority = 0, bool optional = false, bool hasDefaultValue = false, string defaultValue = "", vector<string>* defaultEquivalence = nullptr) noexcept : shortcut(shortcut), name(name), type(type), description(description), priority(priority), optional(optional), hasDefaultValue(hasDefaultValue), defaultValue(defaultValue), defaultEquivalence(defaultEquivalence) {}
	};

	vector<Argument> arguments;
	unordered_map<string, Argument> nameToArgument;
	unordered_map<char, Argument> shortcutToArgument;
	string argv0;

	Argument const& parseArg(string const& arg) const {
		if (arg.rfind("--", 0) == 0) {
			string name = arg.substr(2);
			if (name == "help") {
				displayHelp(std::cout);
				exit(0);
			}
			if (nameToArgument.count(name) == 0) {
				displayHelp(std::cerr);
				std::cerr << "Error: Unknown argument " << arg << std::endl;
				throw std::invalid_argument("Unknown argument");
			}
			return nameToArgument.at(name);
		}
		else if (arg.rfind("-", 0) == 0 && arg.length() > 1) {
			char shortcut = arg[1];
			if (shortcut == 'h') {
				displayHelp(std::cout);
				exit(0);
			}
			if (shortcutToArgument.count(shortcut) == 0) {
				displayHelp(std::cerr);
				std::cerr << "Error: Unknown argument " << arg << std::endl;
				throw std::invalid_argument("Unknown argument");
			}
			return shortcutToArgument.at(shortcut);
		}
		else {
			displayHelp(std::cerr);
			std::cerr << "Error: Invalid argument format near: " << arg << std::endl;
			throw std::invalid_argument("Invalid argument format");
		}
	}

	template<bool needParseArg = true> void parse(const vector<string>& argv) {
		size_t argc = argv.size();
		for (int i = 0; i < argc; ++i) {
			if (needParseArg && i + 1 == argc && argv[i].size() > 0 && argv[i][0] != '-') {
				set("input", argv[i]);
				log << "Warning: " << argv[i] << " is interpreted as: -i " << argv[i] << std::endl;
				break;
			}
			Argument const& argument = (needParseArg) ? parseArg(argv[i]) : nameToArgument.at(argv[i]);
			if (argument.type == "flag") {
				set(argument.name, true);
				if (argument.defaultEquivalence) parse<false>(*argument.defaultEquivalence);
			}
			else if (i + 1 < argc) {
				string value = argv[++i];
				if (argument.type == "integer") {
					size_t realValue;
					try {
						realValue = std::stoull(value);
					}
					catch (const std::exception& e) {
						displayHelp(std::cerr);
						std::cerr << "Error: Invalid integer value for argument " << argument.name << ": " << value << std::endl;
						throw std::invalid_argument("Invalid integer value");
					}
					set(argument.name, realValue);
				}
				else if (argument.type == "numeric") {
					double realValue;
					try {
						realValue = std::stod(value);
					}
					catch (const std::exception& e) {
						displayHelp(std::cerr);
						std::cerr << "Error: Invalid numeric value for argument " << argument.name << ": " << value << std::endl;
						throw std::invalid_argument("Invalid numeric value");
					}
					set(argument.name, realValue);
				}
				else if (argument.type == "string") {
					set(argument.name, value);
				}
				else {
					std::cerr << "Error: Unknown argument type for argument " << argument.name << ". Please fill a bug report! Many thanks!" << std::endl;
					throw std::invalid_argument("Unknown argument type");
				}
			}
			else {
				displayHelp(std::cerr);
				std::cerr << "Error: Missing value for argument " << argument.name << std::endl;
				throw std::invalid_argument("Missing value for argument");
			}
		}
	}

public:
	InputParser() { verbose = 0; }

	template<typename... Args> void addArgument(Args... args) noexcept requires requires { {Argument{ std::forward<Args>(args)... } }; } {
		Argument argument(std::forward<Args>(args)...);
		arguments.push_back(argument);
		nameToArgument[argument.name] = argument;
		if (argument.shortcut > 0) shortcutToArgument[argument.shortcut] = argument;
	}
	
	void displayHelp(std::ostream& out) const {
		out << get<string>("FULL_NAME") << " " << ChangeLog::getGlobalVersion() << std::endl;
		out << "Available arguments:" << std::endl;
		for (Argument const& argument : arguments) {
			out << "  ";
			if (argument.shortcut > 0) out << "-" << argument.shortcut << ", ";
			else out << "    ";
			out << "--" << argument.name << " <" << argument.type << "> : " << argument.description;
			if (argument.optional) out << " (optional";
			else out << " (required";
			if (argument.hasDefaultValue) out << ", default=" << argument.defaultValue;
			out << ")" << std::endl;
		}
	}

	void sortArguments() {
		std::sort(arguments.begin(), arguments.end());
		vector<Argument*> argumentPtrs;
		for (Argument& argument : arguments) argumentPtrs.push_back(&argument);
		vector<Argument> newArguments;
		std::sort(argumentPtrs.begin(), argumentPtrs.end(), [](Argument* a, Argument* b) -> bool { return *a < *b; });
		for (Argument* argumentPtr : argumentPtrs) newArguments.push_back(*argumentPtr);
		arguments = newArguments;
	}

	void parse(int argc, char* argv[]) {
		if (argc == 1) {
			displayHelp(std::cout);
			exit(0);
		}
		argv0 = argv[0];
		vector<string> args;
		for (int i = 1; i < argc; ++i) args.push_back(string(argv[i]));
		parse(args);
		sortArguments();
		for (Argument const& argument : arguments) {
			if (!has(argument.name)) {
				if (argument.optional && argument.hasDefaultValue) {
					if (argument.type == "integer") {
						try {
							size_t defaultValue = std::stoull(argument.defaultValue);
							set(argument.name, defaultValue);
						}
						catch (const std::exception& e) {
							std::cerr << "Error: Invalid default integer value for argument " << argument.name << ": " << argument.defaultValue << ". Please fill a bug report! Many thanks!" << std::endl;
							throw std::invalid_argument("Invalid default integer value");
						}
					}
					else if (argument.type == "numeric") {
						try {
							double defaultValue = std::stod(argument.defaultValue);
							set(argument.name, defaultValue);
						}
						catch (const std::exception& e) {
							std::cerr << "Error: Invalid default numeric value for argument " << argument.name << ": " << argument.defaultValue << ". Please fill a bug report! Many thanks!" << std::endl;
							throw std::invalid_argument("Invalid default numeric value");
						}
					}
					else if (argument.type == "string") set(argument.name, argument.defaultValue);
					else {
						displayHelp(std::cerr);
						std::cerr << "Error: Unknown argument type for argument " << argument.name << ". Please fill a bug report! Many thanks!" << std::endl;
						throw std::invalid_argument("Unknown argument type");
					}
				}
				else if (!argument.optional) {
					displayHelp(std::cerr);
					std::cerr << "Error: Missing required argument --" << argument.name << std::endl;
					throw std::invalid_argument("Missing required argument");
				}
			}
		}
	}

	void print() {
		log << get<string>("FULL_NAME") << " " << ChangeLog::getGlobalVersion() << std::endl;
		log << argv0;
		for (Argument const& argument : arguments) {
			if (has(argument.name)) {
				log << " --" << argument.name;
				if (argument.type == "integer") log << " " << get<size_t>(argument.name);
				else if (argument.type == "numeric") log << " " << get<double>(argument.name);
				else if (argument.type == "string") log << " " << get<string>(argument.name);
				else if (argument.type != "flag") {
					std::cerr << "Error: Unknown argument type for argument " << argument.name << ". Please fill a bug report! Many thanks!" << std::endl;
					throw std::invalid_argument("Unknown argument type");
				}
			}
		}
		log << std::endl;
	}
};

};
common::InputParser	ARG;

namespace common{

class TaxonName2ID {
	vector<string> names;
	unordered_map<string, size_t> name2id;
	unordered_map<string, string> realName;
	bool mapParsed = false;

public:
	void parseMap() {
		if (mapParsed) return;
		mapParsed = true;
		if (!ARG.has("mapping")) return;
		std::ifstream fin(ARG.get<string>("mapping"));
		if (!fin.is_open()) {
			std::cerr << "Error: Unable to open mapping file: " << ARG.get<string>("mapping") << std::endl;
			throw std::invalid_argument("Unable to open mapping file");
		}
		string geneName, speciesName;
		while (fin >> geneName) {
			fin >> speciesName;
			realName[geneName] = speciesName;
		}
	}

	size_t operator[](string const& name) {
		parseMap();
		string real_name = realName.contains(name) ? realName.at(name) : name;
		if (!name2id.contains(real_name)) {
			name2id[real_name] = names.size();
			names.push_back(real_name);
		}
		return name2id.at(real_name);
	}

	string operator[](size_t iTaxon) const noexcept{
		return names[iTaxon];
	}
	size_t nTaxa() const noexcept {
		return names.size();
	}
} taxonName2ID;

template<template<typename> typename T, typename... Args> concept ATTRIBUTES_DISPLAYABLE = requires(Attributes attrs, std::ostream &out, T<Args> const &... args){
	{attrs.template displayAttributes<Args...>(out, std::forward<T<Args> const &>(args)...)} -> std::same_as<std::ostream &>;
};

class AnnotatedBinaryTree : public Attributes{
public:
	static inline string const NUM_LEAVES = "NUM_LEAVES";
	static inline string const LEAF_ID = "LEAF_ID";
	static inline string const TRIPARTITION_SCORE = "TRIPARTITION_SCORE";

	class Node : public Attributes{
		unique_ptr<Node> lc, rc;
		Node* p = nullptr;
	
		Node* placeAbove(unique_ptr<Node> &newChild) noexcept{
			unique_ptr<Node> newParent(new Node());
			Node* oldP = p;
			Node* newC = newChild.get();
			newChild->p = newParent.get();
			newParent->p = oldP;
			p = newParent.get();
			newParent->rc.swap(newChild);
			// newParent: lc=null rc=theNewChild; newChild=null
			
			if (oldP->lc.get() == this) {
				oldP->lc.swap(newParent->lc);
				// oldP: lc=null rc=sister; newParent: lc=this rc=theNewChild; newChild=null
				oldP->lc.swap(newParent);
				// oldP: lc=realNewParent rc=sister; newParent=null; newChild=null
			}
			else {
				oldP->rc.swap(newParent->lc);
				oldP->rc.swap(newParent);
			}
			return newC;
		}
	
		template<template<typename> typename T, typename... Args> requires ATTRIBUTES_DISPLAYABLE<T, Args...> std::ostream &displaySubtree(std::ostream &out, T<Args> const &... args) const{
			if (!isLeaf()) {
				out << "(";
				lc->displaySubtree<T, Args...>(out, std::forward<T<Args> const &>(args)...);
				out << ",";
				rc->displaySubtree<T, Args...>(out, std::forward<T<Args> const &>(args)...);
				out << ")";
			}
			out << "'[";
			displayAttributes<Args...>(out, std::forward<T<Args> const &>(args)...);
			out << "]'";
			return out;
		}
	
		template<typename... Args> requires std::constructible_from<Attributes, Args...> Node(Args... args) noexcept : Attributes(std::forward<Args>(args)...) {}

	public:

		bool isLeaf() const noexcept{ return !lc; }
		
		bool isRoot() const noexcept{ return p->p == nullptr; }
		
		Node* parent() const noexcept{ return (isRoot()) ? nullptr : p; }

		Node* leftChild() const noexcept{ return lc.get(); }

		Node* rightChild() const noexcept{ return rc.get(); }
		
		void swapChildren() noexcept { lc.swap(rc); }
		
		template<typename... Args> Node* emplaceAbove(Args... args) noexcept requires requires{ {Node{ std::forward<Args>(args)... } } noexcept; } {
			unique_ptr<Node> newChild(new Node(std::forward<Args>(args)...));
			return placeAbove(newChild);
		}
		
		Node* regraftAbove(AnnotatedBinaryTree &tree) noexcept{
			if (tree.empty()) return nullptr;
			unique_ptr<Node> newChild;
			newChild.swap(tree.rootRef());
			return placeAbove(newChild);
		}
		
		AnnotatedBinaryTree pruneAbove() noexcept{
			AnnotatedBinaryTree tree;
			
			if (isRoot()){
				tree.rootRef().swap(p->lc);
				p = tree.dummy_root.get();
				return tree;
			}
			
			unique_ptr<Node> cur, sister, parent;
			
			if (p->lc.get() == this) {
				cur.swap(p->lc);
				sister.swap(p->rc);
			}
			else {
				cur.swap(p->rc);
				sister.swap(p->lc);
			}
			// p: lc->null rc->null; parent=null
			
			if (p->p->lc.get() == p){
				parent.swap(p->p->lc);
				sister.swap(p->p->lc);
			}
			else{
				parent.swap(p->p->rc);
				sister.swap(p->p->rc);
			}
			// parent: lc->null rc->null; sister=null
			
			tree.rootRef().swap(cur);
			p = tree.dummy_root.get();
			return tree;
		}
		
		template<typename... Args> std::ostream &displaySubtree(std::ostream &out, String<Args> const &... args) const{
			return displaySubtree<String, Args...>(out, std::forward<String<Args> const &>(args)...);
		}
		
		template<typename... Args> std::ostream &displaySubtree(std::ostream &out, StringVector<Args> const &... args) const{
			return displaySubtree<StringVector, Args...>(out, std::forward<StringVector<Args> const &>(args)...);
		}

		template<typename Internal, typename Length, LOG_OR_OSTREAM Ostream> Ostream& displaySimpleNewickSubtree(Ostream &out, string const &internalAttr, string const &lengthAttr) const{
			if (!isLeaf()) {
				out << "(";
				lc->displaySimpleNewickSubtree<Internal, Length>(out, internalAttr, lengthAttr);
				out << ",";
				rc->displaySimpleNewickSubtree<Internal, Length>(out, internalAttr, lengthAttr);
				out << ")";
				if (has<Internal>(internalAttr)) out << get<Internal>(internalAttr);
			}
			else {
				if (has<size_t>(LEAF_ID)) out << taxonName2ID[get<size_t>(LEAF_ID)];
			}
			if (has<Length>(lengthAttr)) out << ":" << get<Length>(lengthAttr);
			return out;
		}
		
		size_t makeLeafHeavyByLeafCount() noexcept{
			size_t leafcnt;
			if (isLeaf()) leafcnt = 1;
			else{
				size_t lcLeafcnt = lc->makeLeafHeavyByLeafCount();
				size_t rcLeafcnt = rc->makeLeafHeavyByLeafCount();
				leafcnt = lcLeafcnt + rcLeafcnt;
				if (lcLeafcnt < rcLeafcnt) swapChildren();
			}
			set(NUM_LEAVES, leafcnt);
			return leafcnt;
		}
		
		void subtreeLeaves(vector<Node*> &leaves) noexcept{
			if (isLeaf()) leaves.push_back(this);
			else{
				lc->subtreeLeaves(leaves);
				rc->subtreeLeaves(leaves);
			}
		}
		
		friend class AnnotatedBinaryTree;
	};
	
private:
	unique_ptr<Node> dummy_root;
	
	template<template<typename> typename T, typename... Args> requires ATTRIBUTES_DISPLAYABLE<T, Args...> std::ostream &display(std::ostream &out, T<Args> const &... args) const{
		out << "[";
		displayAttributes<Args...>(out, std::forward<T<Args> const &>(args)...);
		out << "]";
		if (!empty()) root()->displaySubtree<T, Args...>(out, std::forward<T<Args> const &>(args)...);
		out << ";";
		return out;
	}
	
	unique_ptr<Node> &rootRef() const{
		return dummy_root->lc;
	}

public:
	template<typename... Args> requires std::constructible_from<Attributes, Args...> AnnotatedBinaryTree(Args... args) noexcept: Attributes(std::forward<Args>(args)...), dummy_root(new Node){}
	
	bool empty() const noexcept{ return !rootRef();}

	Node* root() const noexcept{ return rootRef().get(); }
	
	template<typename... Args> Node* emplaceRoot(Args... args) noexcept requires requires{ {Node{ std::forward<Args>(args)... } } noexcept; } {
		if (empty()){
			rootRef().reset(new Node(std::forward<Args>(args)...));
			root()->p = dummy_root.get();
		}
		else root()->emplaceAbove(std::forward<Args>(args)...);
		return root();
	}

	Node* regraftRoot(AnnotatedBinaryTree &tree) noexcept{
		if (empty()){
			rootRef().swap(tree.rootRef());
			root()->p = dummy_root.get();
		}
		else root()->regraftAbove(tree);
		return root();
	}
	
	template<typename... Args> std::ostream &display(std::ostream &out, String<Args> const &... args) const{
		return display<String, Args...>(out, std::forward<String<Args> const &>(args)...);
	}
	
	template<typename... Args> std::ostream &display(std::ostream &out, StringVector<Args> const &... args) const{
		return display<StringVector, Args...>(out, std::forward<StringVector<Args> const &>(args)...);
	}

	template<typename Internal, typename Length, LOG_OR_OSTREAM Ostream> Ostream& displaySimpleNewick(Ostream& out, string const& internalAttr, string const& lengthAttr) const {
		if (!empty()) root()->displaySimpleNewickSubtree<Internal, Length>(out, internalAttr, lengthAttr);
		out << ";";
		return out;
	}
	
	void makeLeafHeavyByLeafCount() noexcept{ if (!empty()) root()->makeLeafHeavyByLeafCount(); }
	
	vector<Node*> leaves() noexcept{
		vector<Node*> res;
		if (!empty()) root()->subtreeLeaves(res);
		return res;
	}
};

template<class G> class Random {
public:
	using Generator = G;

	Generator generator;

	template<typename... Args> Random(Args... args) requires requires { {Generator{ std::forward<Args>(args)... }}; } : generator{ std::forward<Args>(args)... } {}

	vector<size_t> randomTaxonOrder(size_t nTaxa) {
		vector<size_t> indices;
		for (size_t iTaxa : std::views::iota((size_t) 0, nTaxa)) indices.push_back(iTaxa);
		std::shuffle(indices.begin(), indices.end(), generator);
		return indices;
	}
};

};

using common::ChangeLog;

#endif
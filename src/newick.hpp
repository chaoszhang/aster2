#ifndef NEWICK_HPP
#define NEWICK_HPP

#include <cstddef>
#include <exception>
#include <stdexcept>

#include "common.hpp"

namespace newick {

	ChangeLog logNewick("Newick",
		"2026-10-15", "Chao Zhang", "AI generated Newick parser, printing error messages when throwing exceptions", "patch");

	struct Newick {
		struct Node {
			std::string name;
			std::vector<Node*> children;
			Node* lc = nullptr;
			Node* rc = nullptr;

			bool isLeaf() const noexcept { return children.empty() && !lc; }
		};

	private:
		std::vector<std::unique_ptr<Node> > nodes_;
		Node* root_ = nullptr;
		bool resolved_ = false;

		static bool isSpace(char c) noexcept {
			return c == ' ' || c == '\t' || c == '\n' || c == '\r';
		}

		Node* newNode() {
			nodes_.push_back(std::unique_ptr<Node>(new Node()));
			return nodes_.back().get();
		}

		static std::string readName(std::string const& s, std::size_t& i) {
			std::string name;
			if (i < s.size() && (s[i] == '\'' || s[i] == '"')) {
				char const quote = s[i++];
				while (i < s.size() && s[i] != quote) name += s[i++];
				if (i < s.size()) ++i;
			}
			else {
				while (i < s.size()) {
					char const c = s[i];
					if (c == ';' || c == ')' || c == ',' || c == '(' ||
						c == ':' || c == '[' || isSpace(c)) break;
					name += c;
					++i;
				}
			}
			while (i < s.size() && isSpace(s[i])) ++i;
			return name;
		}

		static void skipTrailing(std::string const& s, std::size_t& i) {
			for (;;) {
				if (i < s.size() && s[i] == ':') {
					++i;
					while (i < s.size() && s[i] != ',' && s[i] != ')' &&
						s[i] != ';' && s[i] != '[') ++i;
				}
				else if (i < s.size() && s[i] == '[') {
					int depth = 0;
					while (i < s.size()) {
						if (s[i] == '[') ++depth;
						else if (s[i] == ']') { --depth; if (depth == 0) { ++i; break; } }
						++i;
					}
				}
				else if (i < s.size() && isSpace(s[i])) ++i;
				else break;
			}
		}

		Node* parseSubtree(std::string const& s, std::size_t& i) {
			while (i < s.size() && isSpace(s[i])) ++i;
			Node* const node = newNode();
			if (i < s.size() && s[i] == '(') {
				++i;
				for (;;) {
					node->children.push_back(parseSubtree(s, i));
					while (i < s.size() && isSpace(s[i])) ++i;
					if (i < s.size() && s[i] == ',') { ++i; continue; }
					break;
				}
				if (i < s.size() && s[i] == ')') ++i;
				else throw std::invalid_argument("Newick: unbalanced parentheses");
				node->name = readName(s, i);
				skipTrailing(s, i);
			}
			else {
				node->name = readName(s, i);
				skipTrailing(s, i);
			}
			return node;
		}

		Node* buildBalanced(std::vector<Node*> const& kids, std::size_t a, std::size_t b) {
			if (b - a == 1) return kids[a];
			std::size_t const h = a + (b - a) / 2;
			Node* const n = newNode();
			n->lc = buildBalanced(kids, a, h);
			n->rc = buildBalanced(kids, h, b);
			return n;
		}

		void resolveNode(Node* n) {
			if (n->children.empty()) return;
			for (Node* const c : n->children) resolveNode(c);
			std::vector<Node*> const kids = std::move(n->children);
			n->children.clear();
			if (kids.size() == 1) {
				throw std::invalid_argument("Newick: unary node is not supported");
			}
			if (kids.size() == 2) {
				n->lc = kids[0];
				n->rc = kids[1];
			}
			else {
				Node* const r = buildBalanced(kids, 0, kids.size());
				n->lc = r->lc;
				n->rc = r->rc;
			}
		}

		void assignImpl(std::string const& text) {
			root_ = nullptr;
			nodes_.clear();
			resolved_ = false;
			std::size_t i = 0;
			if (text.empty()) return;
			root_ = parseSubtree(text, i);
		}

	public:
		Newick() noexcept = default;

		explicit Newick(std::string const& text) try {
			assignImpl(text);
		}
		catch (std::exception const& e) {
			std::cerr << "[newick] Newick(text): " << e.what() << std::endl;
			throw;
		}
		catch (...) {
			std::cerr << "[newick] Newick(text): unknown exception" << std::endl;
			throw;
		}

		void assign(std::string const& text) try {
			assignImpl(text);
		}
		catch (std::exception const& e) {
			std::cerr << "[newick] assign: " << e.what() << std::endl;
			throw;
		}
		catch (...) {
			std::cerr << "[newick] assign: unknown exception" << std::endl;
			throw;
		}

		static std::vector<Newick> readFile(std::string const& path) try {
			std::ifstream fin(path);
			if (!fin.is_open()) {
				throw std::invalid_argument("Newick: unable to open file: " + path);
			}
			std::vector<Newick> trees;
			std::string line;
			while (std::getline(fin, line)) {
				if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;
				trees.push_back(Newick(line));
				trees.back().arbitrarilyResolve();
			}
			return trees;
		}
		catch (std::exception const& e) {
			std::cerr << "[newick] readFile: " << e.what() << std::endl;
			throw;
		}
		catch (...) {
			std::cerr << "[newick] readFile: unknown exception" << std::endl;
			throw;
		}

		bool empty() const noexcept { return !root_; }
		Node* root() const noexcept { return root_; }
		bool resolved() const noexcept { return resolved_; }

		void arbitrarilyResolve()
			try {
				if (!root_) return;
				resolveNode(root_);
				resolved_ = true;
			}
		catch (std::exception const& e) {
			std::cerr << "[newick] arbitrarilyResolve: " << e.what() << std::endl;
			throw;
		}
		catch (...) {
			std::cerr << "[newick] arbitrarilyResolve: unknown exception" << std::endl;
			throw;
		}

		std::size_t nLeaves() const noexcept {
			std::size_t n = 0;
			std::vector<Node const*> stack;
			if (root_) stack.push_back(root_);
			while (!stack.empty()) {
				Node const* const cur = stack.back();
				stack.pop_back();
				if (cur->isLeaf()) ++n;
				else {
					if (cur->lc) stack.push_back(cur->lc);
					if (cur->rc) stack.push_back(cur->rc);
					for (Node* const c : cur->children) stack.push_back(c);
				}
			}
			return n;
		}
	};

	namespace detail {

		inline common::AnnotatedBinaryTree::BinaryCladogram* buildCladogram(
			Newick::Node const* const src,
			std::vector<std::unique_ptr<common::AnnotatedBinaryTree::BinaryCladogram> >& arena) {
			arena.push_back(std::make_unique<common::AnnotatedBinaryTree::BinaryCladogram>());
			common::AnnotatedBinaryTree::BinaryCladogram* const clado = arena.back().get();
			if (src->lc == nullptr) {
				clado->name = src->name;
			}
			else {
				clado->lc = buildCladogram(src->lc, arena);
				clado->rc = buildCladogram(src->rc, arena);
			}
			return clado;
		}

	}

	inline std::vector<common::AnnotatedBinaryTree> readGuideTrees(std::string const& path) try {
		std::vector<common::AnnotatedBinaryTree> trees;
		if (path.empty()) return trees;
		std::vector<Newick> const newicks = Newick::readFile(path);
		trees.reserve(newicks.size());
		std::vector<std::unique_ptr<common::AnnotatedBinaryTree::BinaryCladogram> > arena;
		for (Newick const& nwk : newicks) {
			arena.clear();
			trees.emplace_back(*detail::buildCladogram(nwk.root(), arena));
		}
		return trees;
	}
	catch (std::exception const& e) {
		std::cerr << "[newick] readGuideTrees: " << e.what() << std::endl;
		throw;
	}
	catch (...) {
		std::cerr << "[newick] readGuideTrees: unknown exception" << std::endl;
		throw;
	}

}

#endif // !NEWICK_HPP

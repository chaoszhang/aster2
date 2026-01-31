#ifndef CONSTRAINED_DP_ALGORITHM_HPP
#define CONSTRAINED_DP_ALGORITHM_HPP

#include "common.hpp"

namespace constrained_dp_algorithm {

using std::size_t;
using std::array;
using std::vector;
using std::unordered_map;
using std::views::iota;

template<typename score_t, typename Random> class ConstrainedDP : public common::LogInfo {
	struct NodeHash {
		array<size_t, 8> arr;

		NodeHash() noexcept {
			arr.fill(0);
		}

		NodeHash(typename Random::Generator& generator) noexcept {
			std::uniform_int_distribution<size_t> dist;
			for (size_t& x : arr) {
				size_t a = dist(generator);
				size_t b = dist(generator);
				size_t c = dist(generator);
				x = (a ^ b) + c;
			}
		}

		bool operator==(const NodeHash& other) const noexcept {
			return arr == other.arr;
		}

		NodeHash operator+(const NodeHash& other) const noexcept {
			NodeHash result = *this;
			for (size_t i = 0; i < arr.size(); ++i) {
				result.arr[i] += other.arr[i];
			}
			return result;
		}

		NodeHash operator-(const NodeHash& other) const noexcept {
			NodeHash result = *this;
			for (size_t i = 0; i < arr.size(); ++i) {
				result.arr[i] -= other.arr[i];
			}
			return result;
		}

		NodeHash& operator+=(const NodeHash& other) noexcept {
			for (size_t i = 0; i < arr.size(); ++i) {
				arr[i] += other.arr[i];
			}
			return *this;
		}

		NodeHash& operator-=(const NodeHash& other) noexcept {
			for (size_t i = 0; i < arr.size(); ++i) {
				arr[i] -= other.arr[i];
			}
			return *this;
		}

		struct HashFunction {
			size_t operator()(const NodeHash& nodeHash) const noexcept {
				return nodeHash.arr[0];
			}
		};
	};

	struct Node {
		NodeHash hash, bestChild;
		vector<std::pair<NodeHash, score_t> > children; // pair of (child hash, score)
		bool isLeaf = false;
		size_t iTaxon, dpRound = 0;
		score_t bestScore;
		
		Node() noexcept {};
		Node(NodeHash const& hash) noexcept : hash(hash) {}
	};

	unordered_map<NodeHash, Node, typename NodeHash::HashFunction> hash2node;
	vector<NodeHash> leafHashes;
	NodeHash totalHash;
	size_t dpRound = 0;

	using Tree = common::AnnotatedBinaryTree;


	void addNode(const NodeHash& hash) noexcept {
		if (!hash2node.contains(hash)) {
			hash2node[hash] = Node(hash);
		}
	}

	NodeHash addSubtree(Tree::Node* treeNode) noexcept {
		if (treeNode->isLeaf()) return leafHashes[treeNode->get<size_t>(Tree::LEAF_ID)];
		NodeHash leftHash = addSubtree(treeNode->leftChild());
		NodeHash rightHash = addSubtree(treeNode->rightChild());
		NodeHash hash = leftHash + rightHash;
		addNode(hash);
		score_t score = treeNode->get<score_t>(Tree::TRIPARTITION_SCORE);
		hash2node[hash].children.emplace_back(leftHash, score);
		if (hash != totalHash) {
			addNode(totalHash - leftHash);
			addNode(totalHash - rightHash);
			hash2node[totalHash - leftHash].children.emplace_back(rightHash, score);
			hash2node[totalHash - rightHash].children.emplace_back(leftHash, score);
			hash2node[totalHash].children.emplace_back(hash, score);
		}
		return hash;
	}

	score_t computeDP(NodeHash const& hash, score_t const& ZERO) noexcept {
		Node& node = hash2node.at(hash);
		if (node.dpRound == dpRound) return node.bestScore;
		node.dpRound = dpRound;
		if (node.isLeaf) return node.bestScore = ZERO;
		//if (node.children.empty()) throw std::logic_error("No children for non-leaf node in DP computation.");
		bool first = true;
		score_t bestScore;
		NodeHash bestChild;
		for (auto const& [childHash, score] : node.children) {
			score_t childScore = computeDP(childHash, ZERO) + computeDP(hash - childHash, ZERO) + score;
			if (first || childScore > bestScore) {
				bestScore = childScore;
				bestChild = childHash;
				first = false;
			}
		}
		node.bestChild = bestChild;
		return node.bestScore = bestScore;
	}

	Tree reconstructSubtree(NodeHash const& hash) noexcept {
		if (hash2node.at(hash).isLeaf) {
			Tree tree;
			tree.emplaceRoot(Tree::LEAF_ID, hash2node.at(hash).iTaxon);
			return tree;
		}
		NodeHash lcHash = hash2node.at(hash).bestChild;
		NodeHash rcHash = hash - lcHash;
		Tree lc = reconstructSubtree(lcHash);
		Tree rc = reconstructSubtree(rcHash);
		lc.root()->regraftAbove(rc);
		return lc;
	}

public:
	ConstrainedDP(Random& random, size_t nTaxa, int verbose = common::LogInfo::DEFAULT_VERBOSE) noexcept : LogInfo(verbose) {
		for (size_t iTaxon : iota((size_t)0, nTaxa)) {
			NodeHash h(random.generator);
			leafHashes.push_back(h);
			addNode(h);
			hash2node[h].isLeaf = true;
			hash2node[h].iTaxon = iTaxon;
			totalHash += h;
		}
		addNode(totalHash);
		for (size_t iTaxon : iota((size_t)0, nTaxa)) {
			addNode(totalHash - leafHashes[iTaxon]);
		}
	}

	void addTree(Tree const& tree) noexcept {
		if (!tree.empty()) {
			addSubtree(tree.root());
		}
		//else throw std::logic_error("Cannot add empty tree to ConstrainedDP.");
	}

	Tree optimalUnrootedTree(score_t const& ZERO, size_t outgroup = 0) noexcept {
		if (dpRound == 0){
			for (NodeHash const& h : leafHashes) {
				hash2node[totalHash].children.emplace_back(h, ZERO);
			}
		}
		dpRound++;
		LogInfo vlog(verbose + 1);
		vlog.log() << "Computing constrained dynamic programming tree (round " << dpRound << ")..." << std::endl;
		computeDP(totalHash, ZERO);
		Tree tree = reconstructSubtree(totalHash - leafHashes[outgroup]);
		tree.emplaceRoot(Tree::LEAF_ID, outgroup);
		return tree;
	}
};

};
#endif
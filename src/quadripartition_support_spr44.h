#ifndef QUADRIPARTION_SUPPORT_SPR44_HPP
#define QUADRIPARTION_SUPPORT_SPR44_HPP

#include "stepwise_colorable.hpp"
#include "spr44.h"
#include "common.hpp"
#include <unordered_map>
#include <tuple>
#include <vector>

namespace quadripartition_support_spr44 {

using std::array;
using std::vector;
using std::string;
using std::pair;
using std::size_t;
using std::unordered_map;
using std::tuple;
using std::make_tuple;

template<typename Support> concept QUADRIPARTITION_SUPPORT_SPR44 = requires(typename Support::Color& color, Support const &support, common::AnnotatedBinaryTree::Node * node, size_t index, common::Random<std::mt19937_64>& random)
{
	requires stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE<typename Support::Color>;
	{ Support::staticInitialize(random, index) } noexcept;
	{ Support::ZERO } noexcept -> std::convertible_to<Support>;
	{ Support::FULL_NAME } noexcept -> std::convertible_to<string>;
	{ Support::map(color, index) } noexcept -> std::convertible_to<array<Support, 3> >;
	{ Support::reduce(support, support) } noexcept -> std::convertible_to<Support>;
	{ Support::annotate(node, support, support, support, string()) } noexcept;
};

template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE C> class QuartetScoreSPR44 {
public:
	using Color = C;
	using score_t = typename Color::score_t;
	static QuartetScoreSPR44 const ZERO;
	static inline string const FULL_NAME = "Quartet Score SPR";

private:
	score_t score;

public:
	template<typename T> static void staticInitialize(common::Random<T>& random, size_t nElmnts) noexcept {}

	QuartetScoreSPR44(score_t score = 0): score(score) {}

	static array<QuartetScoreSPR44, 3> map(Color& color, size_t iElement) noexcept {
		array<score_t, 3> scores = color.elementQuadripartitionScores(iElement);
		QuartetScoreSPR44 s0(scores[0]), s1(scores[1]), s2(scores[2]);
		return {s0, s1, s2};
	}

	static QuartetScoreSPR44 reduce(QuartetScoreSPR44 const& a, QuartetScoreSPR44 const& b) noexcept { return QuartetScoreSPR44( a.score + b.score ); }

	static void annotate(common::AnnotatedBinaryTree::Node* node, QuartetScoreSPR44 const& left, QuartetScoreSPR44 const& right, QuartetScoreSPR44 const& outgroup) noexcept {
		node->set(common::AnnotatedBinaryTree::QUADRIPARTITION_SCORE, left.score);
		node->set(common::AnnotatedBinaryTree::QUADRIPARTITION_ALTERNATIVE_1_SCORE, right.score);
		node->set(common::AnnotatedBinaryTree::QUADRIPARTITION_ALTERNATIVE_2_SCORE, outgroup.score);
		node->set(common::AnnotatedBinaryTree::NORMALIZED_QUADRIPARTITION_SCORE, left.score / (left.score + right.score + outgroup.score));
	}
};
template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> QuartetScoreSPR44<Color> const QuartetScoreSPR44<Color>::ZERO = QuartetScoreSPR44<Color>();

template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE C> class SPR44BlockBootstrap : public common::LogInfo {
public:
	using Color = C;
	using score_t = typename Color::score_t;
	static SPR44BlockBootstrap const ZERO;
	static inline string const FULL_NAME = "SPR44 Block Bootstrap";
	static inline size_t constexpr N_BLOCKS = 50;
	static inline size_t constexpr N_BOOTSTRAPS = 1000;
	static size_t nElements;
	static vector<vector<size_t> > bootstrapSPR44;
	using Tree = common::AnnotatedBinaryTree;

private:
	vector<pair<size_t, score_t> > scores;

public:
	template<typename T> static void staticInitialize(common::Random<T>& random, size_t nElmnts) noexcept {
		nElements = nElmnts;
		std::uniform_int_distribution<size_t> uniform(0, N_BLOCKS - 1);
		for (size_t iBootstrap : std::views::iota((size_t)0, N_BOOTSTRAPS)) {
			bootstrapSPR44.emplace_back();
			for (size_t iBlock : std::views::iota((size_t)0, N_BLOCKS)) {
				bootstrapSPR44.back().emplace_back(uniform(random.generator));
			}
		}
	}

	SPR44BlockBootstrap() {
		for (size_t iBlock : std::views::iota((size_t)0, N_BLOCKS)) {
			scores.emplace_back(iBlock, 0);
		}
	}

	SPR44BlockBootstrap(size_t iElement, score_t score){
		size_t iBlock = iElement * N_BLOCKS / nElements;
		scores.push_back({iBlock, score});
	}

	static array<SPR44BlockBootstrap, 3> map(Color& color, size_t iElement) noexcept {
		array<score_t, 3> scores = color.elementQuadripartitionScores(iElement);
		SPR44BlockBootstrap s0(iElement, scores[0]), s1(iElement, scores[1]), s2(iElement, scores[2]);
		return { s0, s1, s2 };
	}

	static SPR44BlockBootstrap reduce(SPR44BlockBootstrap const& a, SPR44BlockBootstrap const& b) noexcept {
		if (a.scores.size() == 1 && b.scores.size() == 1 && a.scores[0].first == b.scores[0].first) {
			return SPR44BlockBootstrap(a.scores[0].first, a.scores[0].second + b.scores[0].second);
		}
		SPR44BlockBootstrap res;
		if (a.scores.size() == 1) res.scores[a.scores[0].first].second += a.scores[0].second;
		else {
			for (size_t iBlock : std::views::iota((size_t)0, N_BLOCKS)) res.scores[iBlock].second += a.scores[iBlock].second;
		}
		if (b.scores.size() == 1) res.scores[b.scores[0].first].second += b.scores[0].second;
		else {
			for (size_t iBlock : std::views::iota((size_t)0, N_BLOCKS)) res.scores[iBlock].second += b.scores[iBlock].second;
		}
		return res;
	}

	static void annotate(common::AnnotatedBinaryTree::Node* node, SPR44BlockBootstrap const& left, SPR44BlockBootstrap const& right, SPR44BlockBootstrap const& outgroup, string key1) noexcept {
		for (size_t iBootstrap : std::views::iota((size_t)0, N_BOOTSTRAPS)) {
			typename Color::score_t s0 = 0, s1 = 0, s2 = 0;
			for (size_t iBlock : bootstrapSPR44[iBootstrap]) {
				s0 += left.scores[iBlock].second;
				s1 += right.scores[iBlock].second;
				s2 += outgroup.scores[iBlock].second;
			}
			std::string key = key1 + std::to_string(iBootstrap);
			std::array<typename Color::score_t, 3> score_arr = { s0, s1, s2 };
			node->set(key, score_arr);
		}
	}
};
template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> SPR44BlockBootstrap<Color> const SPR44BlockBootstrap<Color>::ZERO = SPR44BlockBootstrap<Color>();
template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> vector<vector<size_t> > SPR44BlockBootstrap<Color>::bootstrapSPR44 = vector<vector<size_t> >();
template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> size_t SPR44BlockBootstrap<Color>::nElements = 0;

template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE C> struct ProcedureAttributesSPR44 {
	using Color = C;
	using score_t = Color::score_t;
	template<typename T> using Scheduler = thread_pool::SimpleScheduler<T>;
	template<QUADRIPARTITION_SUPPORT_SPR44 Support> using QuadripartitionScore = spr44_algorithm::StepwiseColorQuadripartitionScore_SPR44<spr44_algorithm::StepwiseColorQuadripartitionScoreDefaultAttributes_SPR44<Color, Support> > ;
	using Random = common::Random<std::mt19937_64>;
};

template<typename Attributes> class ProcedureSPR44 : common::LogInfo {
public:
	using Log = LogInfo;
	using score_t = Attributes::score_t;
	using Tree = common::AnnotatedBinaryTree;
	using Random = Attributes::Random;
	using Color = Attributes::Color;
	using Data = Color::SharedConstData;
	template<typename T> using Scheduler = thread_pool::SimpleScheduler<T>;
	template<typename T> using ThreadPool3 = thread_pool::ThreadPool<std::array<T, 3>, Scheduler>;
	template<QUADRIPARTITION_SUPPORT_SPR44 Support> using QuadripartitionScore = typename Attributes::QuadripartitionScore<Support>;

	template<typename Support> static void annotate(Color& color, Data const& data, Tree& tree, size_t nThreads, int verbose) {
		Log log(verbose);
		log.log() << "Annotating " << Support::FULL_NAME << " ..." << std::endl;
		
		Random random;
		Support::staticInitialize(random, data.nElements());
		ThreadPool3<Support> threadpool(nThreads, 0, data.nElements());
		QuadripartitionScore<Support> qs(color, threadpool, tree, verbose);
		qs.labelTree(SPR44BlockBootstrap<Color>::N_BOOTSTRAPS);
	}

	static void annotate(vector<string> const& supports, Data const& data, Tree& tree, size_t nThreads, int verbose) {
		Color color(data);
		return annotate<SPR44BlockBootstrap<Color> >(color, data, tree, nThreads, verbose);
	}
};

};

#endif


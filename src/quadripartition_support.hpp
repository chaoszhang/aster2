#ifndef QUADRIPARTION_SUPPORT_HPP
#define QUADRIPARTION_SUPPORT_HPP

#include "stepwise_colorable.hpp"
#include "nni_algorithm.hpp"

namespace quadripartition_support {

using std::array;
using std::vector;
using std::string;
using std::pair;
using std::size_t;

template<typename Support> concept QUADRIPARTITION_SUPPORT = requires(typename Support::Color& color, Support const &support, common::AnnotatedBinaryTree::Node * node, size_t index, common::Random<std::mt19937_64>& random)
{
	requires stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE<typename Support::Color>;
	{ Support::staticInitialize(random, index) } noexcept;
	{ Support::ZERO } noexcept -> std::convertible_to<Support>;
	{ Support::FULL_NAME } noexcept -> std::convertible_to<string>;
	{ Support::map(color, index) } noexcept -> std::convertible_to<array<Support, 3> >;
	{ Support::reduce(support, support) } noexcept -> std::convertible_to<Support>;
	{ Support::annotate(node, support, support, support) } noexcept;
};

template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE C> class QuartetScore {
public:
	using Color = C;
	using score_t = typename Color::score_t;
	static QuartetScore const ZERO;
	static inline string const FULL_NAME = "Quartet Score";

private:
	score_t score;

public:
	template<typename T> static void staticInitialize(common::Random<T>& random, size_t nElmnts) noexcept {}

	QuartetScore(score_t score = 0): score(score) {}

	static array<QuartetScore, 3> map(Color& color, size_t iElement) noexcept {
		array<score_t, 3> scores = color.elementQuadripartitionScores(iElement);
		QuartetScore s0(scores[0]), s1(scores[1]), s2(scores[2]);
		return {s0, s1, s2};
	}

	static QuartetScore reduce(QuartetScore const& a, QuartetScore const& b) noexcept { return QuartetScore( a.score + b.score ); }

	static void annotate(common::AnnotatedBinaryTree::Node* node, QuartetScore const& left_right, QuartetScore const& left_outgroup, QuartetScore const& right_outgroup) noexcept {
		node->set(common::AnnotatedBinaryTree::QUADRIPARTITION_SCORE, left_right.score);
		node->set(common::AnnotatedBinaryTree::QUADRIPARTITION_ALTERNATIVE_1_SCORE, left_outgroup.score);
		node->set(common::AnnotatedBinaryTree::QUADRIPARTITION_ALTERNATIVE_2_SCORE, right_outgroup.score);
		node->set(common::AnnotatedBinaryTree::NORMALIZED_QUADRIPARTITION_SCORE, left_right.score / (left_right.score + left_outgroup.score + right_outgroup.score));
	}
};
template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> QuartetScore<Color> const QuartetScore<Color>::ZERO = QuartetScore<Color>();

ChangeLog logLocalBlockBootstrap("LocalBlockBootstrap",
	"2026-02-07", "Chao Zhang", "First version", "minor");

template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE C> class LocalBlockBootstrap {
public:
	using Color = C;
	using score_t = typename Color::score_t;
	static LocalBlockBootstrap const ZERO;
	static inline string const FULL_NAME = "Local Block Bootstrap";
	static inline size_t constexpr N_BLOCKS = 50;
	static inline size_t constexpr N_BOOTSTRAPS = 1000;
	static size_t nElements;
	static vector<vector<size_t> > bootstrap;

private:
	vector<pair<size_t, score_t> > scores;

public:
	template<typename T> static void staticInitialize(common::Random<T>& random, size_t nElmnts) noexcept {
		nElements = nElmnts;
		std::uniform_int_distribution<size_t> uniform(0, N_BLOCKS - 1);
		for (size_t iBootstrap : std::views::iota((size_t)0, N_BOOTSTRAPS)) {
			bootstrap.emplace_back();
			for (size_t iBlock : std::views::iota((size_t)0, N_BLOCKS)) {
				bootstrap.back().emplace_back(uniform(random.generator));
			}
		}
	}

	LocalBlockBootstrap() {
		for (size_t iBlock : std::views::iota((size_t)0, N_BLOCKS)) {
			scores.emplace_back(iBlock, 0);
		}
	}

	LocalBlockBootstrap(size_t iElement, score_t score){
		size_t iBlock = iElement * N_BLOCKS / nElements;
		scores.push_back({iBlock, score});
	}

	static array<LocalBlockBootstrap, 3> map(Color& color, size_t iElement) noexcept {
		array<score_t, 3> scores = color.elementQuadripartitionScores(iElement);
		LocalBlockBootstrap s0(iElement, scores[0]), s1(iElement, scores[1]), s2(iElement, scores[2]);
		return { s0, s1, s2 };
	}

	static LocalBlockBootstrap reduce(LocalBlockBootstrap const& a, LocalBlockBootstrap const& b) noexcept {
		if (a.scores.size() == 1 && b.scores.size() == 1 && a.scores[0].first == b.scores[0].first) {
			return LocalBlockBootstrap(a.scores[0].first, a.scores[0].second + b.scores[0].second);
		}
		LocalBlockBootstrap res;
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

	static void annotate(common::AnnotatedBinaryTree::Node* node, LocalBlockBootstrap const& left_right, LocalBlockBootstrap const& left_outgroup, LocalBlockBootstrap const& right_outgroup) noexcept {
		size_t cnt = 0;
		for (size_t iBootstrap : std::views::iota((size_t)0, N_BOOTSTRAPS)) {
			score_t s0 = 0, s1 = 0, s2 = 0;
			for (size_t iBlock : bootstrap[iBootstrap]) {
				s0 += left_right.scores[iBlock].second;
				s1 += left_outgroup.scores[iBlock].second;
				s2 += right_outgroup.scores[iBlock].second;
			}
			if (s0 > s1 && s0 > s2) cnt++;
		}
		if (cnt == 0) {
			for (size_t iBootstrap : std::views::iota((size_t)0, N_BOOTSTRAPS)) {
				score_t s0 = 0, s1 = 0, s2 = 0;
				for (size_t iBlock : bootstrap[iBootstrap]) {
					s0 += left_right.scores[iBlock].second;
					s1 += left_outgroup.scores[iBlock].second;
					s2 += right_outgroup.scores[iBlock].second;
				}
			}
		}
		node->set("LocalBlockBootstrap", cnt * 100.0 / N_BOOTSTRAPS);
	}
};
template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> LocalBlockBootstrap<Color> const LocalBlockBootstrap<Color>::ZERO = LocalBlockBootstrap<Color>();
template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> vector<vector<size_t> > LocalBlockBootstrap<Color>::bootstrap = vector<vector<size_t> >();
template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> size_t LocalBlockBootstrap<Color>::nElements = 0;

template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE C> struct ProcedureAttributes {
	using Color = C;
	using score_t = Color::score_t;
	template<typename T> using Scheduler = thread_pool::SimpleScheduler<T>;
	template<QUADRIPARTITION_SUPPORT Support> using QuadripartitionScore = nni_algorithm::StepwiseColorQuadripartitionScore<nni_algorithm::StepwiseColorQuadripartitionScoreDefaultAttributes<Color, Support> > ;
	using Random = common::Random<std::mt19937_64>;
};

template<typename Attributes> class Procedure : common::LogInfo {
public:
	using Log = LogInfo;
	using score_t = Attributes::score_t;
	using Tree = common::AnnotatedBinaryTree;
	using Random = Attributes::Random;
	using Color = Attributes::Color;
	using Data = Color::SharedConstData;
	template<typename T> using Scheduler = thread_pool::SimpleScheduler<T>;
	template<typename T> using ThreadPool3 = thread_pool::ThreadPool<std::array<T, 3>, Scheduler>;
	template<QUADRIPARTITION_SUPPORT Support> using QuadripartitionScore = typename Attributes::QuadripartitionScore<Support>;

	template<typename Support> static void annotate(Color& color, Data const& data, Tree& tree, size_t nThreads, int verbose) {
		Log log(verbose);
		log.log() << "Annotating " << Support::FULL_NAME << " ..." << std::endl;
		
		Random random;
		Support::staticInitialize(random, data.nElements);
		ThreadPool3<Support> threadpool(nThreads, 0, data.nElements);
		QuadripartitionScore<Support> qs(color, threadpool, tree, verbose);
		qs.labelTree();
	}

	static void annotate(vector<string> const& supports, Data const& data, Tree& tree, size_t nThreads, int verbose) {
		Color color(&data);
		return annotate<LocalBlockBootstrap<Color> >(color, data, tree, nThreads, verbose);
	}
};

};

#endif
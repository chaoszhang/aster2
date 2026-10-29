#ifndef QUADRIPARTION_SUPPORT_HPP
#define QUADRIPARTION_SUPPORT_HPP

#include "branch_annotation.hpp"

namespace quadripartition_support {

using std::array;
using std::vector;
using std::string;
using std::pair;
using std::size_t;

template<typename Support> concept QUADRIPARTITION_SUPPORT = requires{
	{ Support::SHORT_NAME } noexcept -> std::convertible_to<string>;
	{ Support::FULL_NAME } noexcept -> std::convertible_to<string>;
	{ Support::MAIN_KEY } noexcept -> std::convertible_to<string>;
	{ Support::KEYS } noexcept -> std::convertible_to<vector<string> const>;
};

template<typename Support, typename Color> concept QUADRIPARTITION_SUPPORT_STEPWISE_COLORABLE = requires(Color & color, Support::template Data<Color> const& support, common::AnnotatedBinaryTree::Node * node, size_t index, common::Random<std::mt19937_64>& random)
{
	requires QUADRIPARTITION_SUPPORT<Support>;
	requires stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE<Color>;
	{ Support::template Data<Color>::ZERO } noexcept -> std::convertible_to<typename Support::template Data<Color> >;
	{ Support::template Data<Color>::staticInitialize(random, index) } noexcept;
	{ Support::template Data<Color>::map(color, index) } noexcept -> std::convertible_to<array<typename Support::template Data<Color>, 3> >;
	{ Support::template Data<Color>::reduce(support, support) } noexcept -> std::convertible_to<typename Support::template Data<Color> >;
	{ Support::template Data<Color>::annotate(node, support, support, support) } noexcept;
};

ChangeLog logQUADRIPARTITION_SUPPORT_POSTPROCESS("QUADRIPARTITION_SUPPORT_POSTPROCESS",
	"2026-09-03", "Chao Zhang", "Supporting post-processing", "patch");

template<typename Support> concept QUADRIPARTITION_SUPPORT_POSTPROCESS = requires(common::AnnotatedBinaryTree & tree)
{
	requires QUADRIPARTITION_SUPPORT<Support>;
	{ Support::postProcess(tree) } noexcept;
};

class QuartetScore {
public:
	static inline string const SHORT_NAME = "qs";
	static inline string const FULL_NAME = "Quartet Score";
	static inline string const MAIN_KEY = common::AnnotatedBinaryTree::NORMALIZED_QUADRIPARTITION_SCORE;
	static inline vector<string> const KEYS = { common::AnnotatedBinaryTree::QUADRIPARTITION_SCORE, common::AnnotatedBinaryTree::QUADRIPARTITION_ALTERNATIVE_1_SCORE, common::AnnotatedBinaryTree::QUADRIPARTITION_ALTERNATIVE_2_SCORE, common::AnnotatedBinaryTree::NORMALIZED_QUADRIPARTITION_SCORE };

	template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE C> class Data{
	public:
		using Color = C;
		using score_t = typename Color::score_t;
		static Data const ZERO;

	private:
		score_t score;

	public:
		template<typename T> static void staticInitialize(common::Random<T>& random, size_t nElmnts) noexcept {}

		Data(score_t score = Color::ZERO): score(score) {}

		static array<Data, 3> map(Color& color, size_t iElement) noexcept {
			array<score_t, 3> scores = color.elementQuadripartitionScores(iElement);
			Data s0(scores[0]), s1(scores[1]), s2(scores[2]);
			return {s0, s1, s2};
		}

		static Data reduce(Data const& a, Data const& b) noexcept { return Data( a.score + b.score ); }

		static void annotate(common::AnnotatedBinaryTree::Node* node, Data const& left_right, Data const& left_outgroup, Data const& right_outgroup) noexcept {
			node->set(common::AnnotatedBinaryTree::QUADRIPARTITION_SCORE, left_right.score);
			node->set(common::AnnotatedBinaryTree::QUADRIPARTITION_ALTERNATIVE_1_SCORE, left_outgroup.score);
			node->set(common::AnnotatedBinaryTree::QUADRIPARTITION_ALTERNATIVE_2_SCORE, right_outgroup.score);
			score_t const total = left_right.score + left_outgroup.score + right_outgroup.score;
			double const normalized = (total == 0) ? 0.0 : static_cast<double>(left_right.score) / static_cast<double>(total);
			node->set(common::AnnotatedBinaryTree::NORMALIZED_QUADRIPARTITION_SCORE, normalized);
		}
	};
};
template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> QuartetScore::Data<Color> const QuartetScore::Data<Color>::ZERO = QuartetScore::Data<Color>();

ChangeLog logLocalBlockBootstrap("LocalBlockBootstrap",
	"2026-02-07", "Chao Zhang", "First version", "minor");

class LocalBlockBootstrap {
public:
	static inline string const SHORT_NAME = "lbb";
	static inline string const FULL_NAME = "Local Block Bootstrap";
	static inline string const MAIN_KEY = "LocalBlockBootstrap";
	static inline vector<string> const KEYS = { "LocalBlockBootstrap" };

	template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE C> class Data {
	public:
		using Color = C;
		using score_t = typename Color::score_t;
		static Data const ZERO;

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

		Data() {
			for (size_t iBlock : std::views::iota((size_t)0, N_BLOCKS)) {
				scores.emplace_back(iBlock, 0);
			}
		}

		Data(size_t iElement, score_t score) {
			size_t iBlock = iElement * N_BLOCKS / nElements;
			scores.push_back({ iBlock, score });
		}

		static array<Data, 3> map(Color& color, size_t iElement) noexcept {
			array<score_t, 3> scores = color.elementQuadripartitionScores(iElement);
			Data s0(iElement, scores[0]), s1(iElement, scores[1]), s2(iElement, scores[2]);
			return { s0, s1, s2 };
		}

		static Data reduce(Data const& a, Data const& b) noexcept {
			if (a.scores.size() == 1 && b.scores.size() == 1 && a.scores[0].first == b.scores[0].first) {
				return Data(a.scores[0].first, a.scores[0].second + b.scores[0].second);
			}
			Data res;
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

		static void annotate(common::AnnotatedBinaryTree::Node* node, Data const& left_right, Data const& left_outgroup, Data const& right_outgroup) noexcept {
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
};
template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> LocalBlockBootstrap::Data<Color> const LocalBlockBootstrap::Data<Color>::ZERO = LocalBlockBootstrap::Data<Color>();
template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> size_t LocalBlockBootstrap::Data<Color>::nElements{};
template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> vector<vector<size_t> > LocalBlockBootstrap::Data<Color>::bootstrap{};

template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE C> struct ProcedureAttributes {
	using Color = C;
	using score_t = Color::score_t;
	template<typename T> using Scheduler = thread_pool::SimpleScheduler<T>;
	template<QUADRIPARTITION_SUPPORT Support> using QuadripartitionScore = nni_algorithm::StepwiseColorQuadripartitionScore<nni_algorithm::StepwiseColorQuadripartitionScoreDefaultAttributes<Color, typename Support::template Data<Color> > > ;
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

	template<typename Support, typename Color> static void annotate(Color& color, Data const& data, Tree& tree, size_t nThreads, int verbose) {
		Log log(verbose);
		log.log() << "Annotating " << Support::FULL_NAME << " ..." << std::endl;
		
		Random random;
		Support::template Data<Color>::staticInitialize(random, data.nElements());
		ThreadPool3<typename Support::template Data<Color> > threadpool(nThreads, 0, data.nElements());
		QuadripartitionScore<Support> qs(color, threadpool, tree, verbose);
		qs.labelTree();
		if constexpr (QUADRIPARTITION_SUPPORT_POSTPROCESS<Support>) Support::postProcess(tree);
	}
};

template<QUADRIPARTITION_SUPPORT S> class BranchAnnotation {
public:
	static std::string const shortName() {
		return S::SHORT_NAME;
	}

	static std::string const fullName() {
		return S::FULL_NAME;
	}

	template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> requires QUADRIPARTITION_SUPPORT_STEPWISE_COLORABLE<S, Color> static std::string const annotate(Color::SharedConstData const& data, common::AnnotatedBinaryTree& tree, size_t nThreads, int verbose) {
		Color color(data);
		Procedure<ProcedureAttributes<Color> >::template annotate<S, Color>(color, data, tree, nThreads, verbose);
		return S::MAIN_KEY;
	}

	template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> requires QUADRIPARTITION_SUPPORT_STEPWISE_COLORABLE<S, Color> static std::vector<std::string> annotateAll(Color::SharedConstData const& data, common::AnnotatedBinaryTree& tree, size_t nThreads, int verbose) {
		Color color(data);
		Procedure<ProcedureAttributes<Color> >::template annotate<S, Color>(color, data, tree, nThreads, verbose);
		return S::KEYS;
	}
};

using BA_QS = BranchAnnotation<QuartetScore>;
using BA_LBB = BranchAnnotation<LocalBlockBootstrap>;

};

#endif
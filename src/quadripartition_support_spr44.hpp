#ifndef QUADRIPARTION_SUPPORT_SPR44_HPP
#define QUADRIPARTION_SUPPORT_SPR44_HPP

#include "stepwise_colorable.hpp"
#include "spr44.hpp"
#include "branch_annotation.hpp"
#include "common.hpp"

namespace quadripartition_support_spr44 {

	using std::array;
	using std::vector;
	using std::string;
	using std::pair;
	using std::size_t;
	using std::unordered_map;
	using std::tuple;
	using std::make_tuple;
	using std::shared_ptr;

	ChangeLog logQuadripartitionSupportSPR("logQuadripartitionSupportSPR",
		"2026-10-16", "Yujie Wang", "add quardripartition support for spr", "minor",
		"2026-10-19", "Chao Zhang", "memory efficiency improvements by DeepSeek", "patch");

	template<typename Support> concept QUADRIPARTITION_SUPPORT_SPR44 = requires(typename Support::Color & color, Support const& support, common::AnnotatedBinaryTree::Node * node, size_t index, common::Random<std::mt19937_64>&random)
	{
		requires stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE<typename Support::Color>;
	{ Support::staticInitialize(random, index) } noexcept;
	{ Support::ZERO } noexcept -> std::convertible_to<Support>;
	{ Support::SHORT_NAME } noexcept -> std::convertible_to<string>;
	{ Support::FULL_NAME } noexcept -> std::convertible_to<string>;
	{ Support::MAIN_KEY } noexcept -> std::convertible_to<string>;
	{ Support::KEYS } noexcept -> std::convertible_to<vector<string> const>;
	{ Support::map(color, index) } noexcept -> std::convertible_to<array<Support, 3> >;
	{ Support::reduce(support, support) } noexcept -> std::convertible_to<Support>;
	{ Support::annotate(node, support, support, support, string()) } noexcept;
	};

	template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE C> class SPR44BlockBootstrap : public common::LogInfo {
	public:
		using Color = C;
		using score_t = typename Color::score_t;
		static SPR44BlockBootstrap const ZERO;
		static inline string const SHORT_NAME = "sbb";
		static inline string const FULL_NAME = "SPR44 Block Bootstrap";
		static inline string const MAIN_KEY = "SPR44BlockBootstrap";
		static inline vector<string> const KEYS = { "SPR44BlockBootstrap" };

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

		SPR44BlockBootstrap(size_t iElement, score_t score) {
			size_t iBlock = iElement * N_BLOCKS / nElements;
			scores.push_back({ iBlock, score });
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
			// ---- EVAL FIX (M1) ----------------------------------------------------
			// 旧实现：每个节点挂一个 std::map<std::string, array<double,3>>，键为
			// "PLACE44<切割号><bootstrap 号>"，于是条目数 = 切割数 × 节点数 × N_BOOTSTRAPS
			// （n=201、1000 次 bootstrap 时 1.6e8 条 ≈ 21 GB，且从不释放）。
			// 新实现：每个节点一份按 bootstrap 下标索引的定长数组，切割之间复用/覆盖；
			// 切割内部的消费（markCutForBootstrap）在同一切割着色后立即完成。
			auto& values = spr44_algorithm::bootstrapValues<typename Color::score_t>(node);
			if (values.size() != N_BOOTSTRAPS) {
				values.assign(N_BOOTSTRAPS, std::array<typename Color::score_t, 3>{ 0, 0, 0 });
			}
			(void)key1;
			for (size_t iBootstrap : std::views::iota((size_t)0, N_BOOTSTRAPS)) {
				typename Color::score_t s0 = 0, s1 = 0, s2 = 0;
				for (size_t iBlock : bootstrapSPR44[iBootstrap]) {
					s0 += left.scores[iBlock].second;
					s1 += right.scores[iBlock].second;
					s2 += outgroup.scores[iBlock].second;
				}
				values[iBootstrap] = { s0, s1, s2 };
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
		template<QUADRIPARTITION_SUPPORT_SPR44 Support> using QuadripartitionScore = spr44_algorithm::StepwiseColorQuadripartitionScore_SPR44<spr44_algorithm::StepwiseColorQuadripartitionScoreDefaultAttributes_SPR44<Color, Support> >;
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

	template<template<class> class S> class BranchAnnotationSPR44 {
	public:
		template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> requires QUADRIPARTITION_SUPPORT_SPR44<S<Color> > static std::string const shortName() {
			return S<Color>::SHORT_NAME;
		}

		template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> requires QUADRIPARTITION_SUPPORT_SPR44<S<Color> > static std::string const fullName() {
			return S<Color>::FULL_NAME;
		}

		template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> requires QUADRIPARTITION_SUPPORT_SPR44<S<Color> > static std::string const annotate(Color::SharedConstData const& data, common::AnnotatedBinaryTree& tree, size_t nThreads, int verbose) {
			Color color(data);
			ProcedureSPR44<ProcedureAttributesSPR44<Color> >::template annotate<S<Color> >(color, data, tree, nThreads, verbose);
			return S<Color>::MAIN_KEY;
		}

		template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> requires QUADRIPARTITION_SUPPORT_SPR44<S<Color> > static std::vector<std::string> annotateAll(Color::SharedConstData const& data, common::AnnotatedBinaryTree& tree, size_t nThreads, int verbose) {
			Color color(data);
			ProcedureSPR44<ProcedureAttributesSPR44<Color> >::template annotate<S<Color> >(color, data, tree, nThreads, verbose);
			return S<Color>::KEYS;
		}
	};

	using BA_SBB = BranchAnnotationSPR44<SPR44BlockBootstrap>;

};

#endif
#ifndef BRANCH_ANNOTATION_HPP
#define BRANCH_ANNOTATION_HPP

#include "common.hpp"
#include "stepwise_colorable.hpp"
#include "nni_algorithm.hpp"

namespace branch_annotation {

template<class BA, class Color> concept BRANCH_ANNOTATION = requires(typename Color::SharedConstData const& data, common::AnnotatedBinaryTree & tree, size_t nThreads, int verbose) {
	requires stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE<Color>;
	{ BA::template shortName<Color>() } -> std::convertible_to<std::string const>;
	{ BA::template fullName<Color>() } -> std::convertible_to<std::string const>;
	{ BA::template annotate<Color>(data, tree, nThreads, verbose) } -> std::convertible_to<std::string const>;
	{ BA::template annotateAll<Color>(data, tree, nThreads, verbose) } -> std::convertible_to<std::vector<std::string> const>;
};

template<class BranchAnnotations, class Color> concept BRANCH_ANNOTATIONS = requires(std::string const& name, typename Color::SharedConstData const& data, common::AnnotatedBinaryTree & tree, size_t nThreads, int verbose) {
	requires stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE<Color>;
	{ BranchAnnotations::template shortNames<Color>() } -> std::convertible_to<std::vector<std::string> const>;
	{ BranchAnnotations::template fullNames<Color>() } -> std::convertible_to<std::vector<std::string> const>;
	{ BranchAnnotations::template isValid<Color>(name) } -> std::convertible_to<bool>;
	{ BranchAnnotations::template annotate<Color>(name, data, tree, nThreads, verbose) } -> std::convertible_to<std::string const>;
	{ BranchAnnotations::template annotateAll<Color>(data, tree, nThreads, verbose) } -> std::convertible_to<std::vector<std::string> const>;
};

ChangeLog logLocalBlockBootstrap("BranchAnnotations",
	"2026-10-06", "Chao Zhang", "Refactoring", "patch");

template<class BA, class... Args> class BranchAnnotations {
public:
	template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> requires BRANCH_ANNOTATION<BA, Color> static std::vector<std::string> const shortNames() {
		std::vector<std::string> res;
		if constexpr (sizeof...(Args) > 0) res = BranchAnnotations<Args...>::template shortNames<Color>();
		res.emplace_back(BA::template shortName<Color>());
		return res;
	}

	template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> requires BRANCH_ANNOTATION<BA, Color> static std::vector<std::string> const fullNames() {
		std::vector<std::string> res;
		if constexpr (sizeof...(Args) > 0) res = BranchAnnotations<Args...>::template fullNames<Color>();
		res.emplace_back(BA::template fullName<Color>());
		return res;
	}

	template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> requires BRANCH_ANNOTATION<BA, Color> static bool isValid(std::string const& name) {
		if (name == BA::template shortName<Color>()) return true;
		if constexpr (sizeof...(Args) > 0) {
			return BranchAnnotations<Args...>::template isValid<Color>(name);
		}
		return false;
	}

	template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> requires BRANCH_ANNOTATION<BA, Color> static std::string const annotate(std::string const& name, Color::SharedConstData const& data, common::AnnotatedBinaryTree& tree, size_t nThreads, int verbose) {
		if (name == BA::template shortName<Color>()) return BA::template annotate<Color>(data, tree, nThreads, verbose);
		if constexpr (sizeof...(Args) > 0) return BranchAnnotations<Args...>::template annotate<Color>(name, data, tree, nThreads, verbose);
		return "";
	}

	template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> requires BRANCH_ANNOTATION<BA, Color> static std::vector<std::string> annotateAll(Color::SharedConstData const& data, common::AnnotatedBinaryTree& tree, size_t nThreads, int verbose) {
		std::vector<std::string> res;
		if constexpr (sizeof...(Args) > 0) res = BranchAnnotations<Args...>::template annotateAll<Color>(data, tree, nThreads, verbose);
		for (std::string const& key : BA::template annotateAll<Color>(data, tree, nThreads, verbose)) res.push_back(key);
		return res;
	}
};

};

#endif // !BRANCH_ANNOTATION_HPP

#ifndef BRANCH_ANNOTATION_HPP
#define BRANCH_ANNOTATION_HPP

#include "common.hpp"
#include "stepwise_colorable.hpp"
#include "nni_algorithm.hpp"

namespace branch_annotation {

template<class BA, class Color> concept BRANCH_ANNOTATION = requires(typename Color::SharedConstData const& data, common::AnnotatedBinaryTree & tree, size_t nThreads, int verbose) {
	requires stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE<Color>;
	{ BA::shortName() } -> std::convertible_to<std::string const>;
	{ BA::fullName() } -> std::convertible_to<std::string const>;
	{ BA::template annotate<Color>(data, tree, nThreads, verbose) } -> std::convertible_to<std::string const>;
	{ BA::template annotateAll<Color>(data, tree, nThreads, verbose) } -> std::convertible_to<std::vector<std::string> const>;
};

template<class BranchAnnotations, class Color> concept BRANCH_ANNOTATIONS = requires(std::string const& name, typename Color::SharedConstData const& data, common::AnnotatedBinaryTree & tree, size_t nThreads, int verbose) {
	requires stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE<Color>;
	{ BranchAnnotations::shortNames() } -> std::convertible_to<std::vector<std::string> const>;
	{ BranchAnnotations::fullNames() } -> std::convertible_to<std::vector<std::string> const>;
	{ BranchAnnotations::isValid(name) } -> std::convertible_to<bool>;
	{ BranchAnnotations::template annotate<Color>(name, data, tree, nThreads, verbose) } -> std::convertible_to<std::string const>;
	{ BranchAnnotations::template annotateAll<Color>(data, tree, nThreads, verbose) } -> std::convertible_to<std::vector<std::string> const>;
};

ChangeLog logLocalBlockBootstrap("BranchAnnotations",
	"2026-10-06", "Chao Zhang", "Refactoring", "patch");

template<class BA, class... Args> class BranchAnnotations {
public:
	static std::vector<std::string> const shortNames() {
		std::vector<std::string> res;
		if constexpr (sizeof...(Args) > 0) res = BranchAnnotations<Args...>::shortNames();
		res.emplace_back(BA::shortName());
		return res;
	}

	static std::vector<std::string> const fullNames() {
		std::vector<std::string> res;
		if constexpr (sizeof...(Args) > 0) res = BranchAnnotations<Args...>::fullNames();
		res.emplace_back(BA::fullName());
		return res;
	}

	static bool isValid(std::string const& name) {
		if (name == BA::shortName()) return true;
		if constexpr (sizeof...(Args) > 0) {
			return BranchAnnotations<Args...>::isValid(name);
		}
		return false;
	}

	static std::string const supportOptionText() {
		std::vector<std::string> sNames = shortNames();
		std::vector<std::string> fNames = fullNames();
		std::string res = "none, verbose";
		for (int i = 0; i < sNames.size(); ++i) {
			res += std::string(", ") + sNames[i] + " [" + fNames[i] + "]";
		}
		return res;
	}

	static std::string const lengthOptionText() {
		std::vector<std::string> sNames = shortNames();
		std::vector<std::string> fNames = fullNames();
		std::string res = "none";
		for (int i = 0; i < sNames.size(); ++i) {
			res += std::string(", ") + sNames[i] + " [" + fNames[i] + "]";
		}
		return res;
	}

	template<stepwise_colorable::QUADRIPARTITION_STEPWISE_COLORABLE Color> requires BRANCH_ANNOTATION<BA, Color> static std::string const annotate(std::string const& name, Color::SharedConstData const& data, common::AnnotatedBinaryTree& tree, size_t nThreads, int verbose) {
		if (name == BA::shortName()) return BA::template annotate<Color>(data, tree, nThreads, verbose);
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

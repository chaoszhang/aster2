#ifndef DRIVER_HPP
#define DRIVER_HPP

#include "common.hpp"
#include "branch_annotation.hpp"

namespace driver {

template<template<bool> class Driver> concept DRIVER = requires {
	{ Driver<true>::programNames() } noexcept -> std::convertible_to<std::pair<std::string, std::string> >;
	{ Driver<true>::addArguments() } noexcept;
	{ Driver<true>::getStepwiseColorSharedConstData() } noexcept;
	{ std::visit([]<typename StepwiseColorSharedConstData>(StepwiseColorSharedConstData const& stepwiseColorSharedConstData) {}, Driver<true>::getStepwiseColorSharedConstData()) };
};

template<template<bool> class Driver> concept ANNOTATION_READY = requires {
	requires DRIVER<Driver>;
	requires branch_annotation::BRANCH_ANNOTATIONS<typename Driver<true>::BranchAnnotations, typename std::variant_alternative_t<0, typename Driver<true>::DataClasses>::ParentClass>;
	{ Driver<true>::DEFAULT_SUPPORT_ANNOTATION() } noexcept -> std::convertible_to<std::string const>;
	{ Driver<true>::DEFAULT_LENGTH_ANNOTATION() } noexcept -> std::convertible_to<std::string const>;
};

template<template<bool> class Driver, class Color> requires ANNOTATION_READY<Driver> std::string const annotate(std::string const& name, typename Color::SharedConstData const& data, common::AnnotatedBinaryTree& tree, size_t nThreads, int verbose) {
	return Driver<true>::BranchAnnotations::template annotate<Color>(name, data, tree, nThreads, 0);
}

template<template<bool> class Driver, class Color> requires (!ANNOTATION_READY<Driver>) std::string const annotate(std::string const& name, typename Color::SharedConstData const& data, common::AnnotatedBinaryTree& tree, size_t nThreads, int verbose) {
	return "";
}

template<template<bool> class Driver> requires ANNOTATION_READY<Driver> bool const isValid(std::string const& name) {
	return Driver<true>::BranchAnnotations::isValid(name);
}

template<template<bool> class Driver> requires (!ANNOTATION_READY<Driver>) bool const isValid(std::string const& name) {
	return false;
}

template<template<bool> class Driver> requires ANNOTATION_READY<Driver> std::string const supportOptionText() {
	return Driver<true>::BranchAnnotations::supportOptionText();
}

template<template<bool> class Driver> requires ANNOTATION_READY<Driver> std::string const DEFAULT_SUPPORT_ANNOTATION() {
	return Driver<true>::DEFAULT_SUPPORT_ANNOTATION();
}

template<template<bool> class Driver> requires ANNOTATION_READY<Driver> std::string const lengthOptionText() {
	return Driver<true>::BranchAnnotations::lengthOptionText();
}

template<template<bool> class Driver> requires ANNOTATION_READY<Driver> std::string const DEFAULT_LENGTH_ANNOTATION() {
	return Driver<true>::DEFAULT_LENGTH_ANNOTATION();
}

template<template<bool> class Driver> requires (!ANNOTATION_READY<Driver>) std::string const supportOptionText() {
	return "";
}

template<template<bool> class Driver> requires (!ANNOTATION_READY<Driver>) std::string const DEFAULT_SUPPORT_ANNOTATION() {
	return "";
}

template<template<bool> class Driver> requires (!ANNOTATION_READY<Driver>) std::string const lengthOptionText() {
	return "";
}

template<template<bool> class Driver> requires (!ANNOTATION_READY<Driver>) std::string const DEFAULT_LENGTH_ANNOTATION() {
	return "";
}

};

#endif // !DRIVER_HPP

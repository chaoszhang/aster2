#ifndef DRIVER_HPP
#define DRIVER_HPP

#include "common.hpp"

namespace driver {

template<template<bool> class Driver> concept DRIVER = requires {
	{ Driver<true>::programNames() } noexcept -> std::convertible_to<std::pair<std::string, std::string> >;
	{ Driver<true>::addArguments() } noexcept;
	{ Driver<true>::getStepwiseColorSharedConstData() } noexcept;
	{ std::visit([]<typename StepwiseColorSharedConstData>(StepwiseColorSharedConstData const& stepwiseColorSharedConstData) {}, Driver<true>::getStepwiseColorSharedConstData()) };
};

};

#endif // !DRIVER_HPP

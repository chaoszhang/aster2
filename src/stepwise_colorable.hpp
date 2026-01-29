#ifndef STEPWISE_COLORABLE_HPP
#define STEPWISE_COLORABLE_HPP

#include "std_include.hpp"

namespace stepwise_colorable{

template<class T> concept STEPWISE_COLORABLE = requires(T t, typename T::SharedConstData const *const dataPtr, typename T::SharedConstData data, std::size_t index, typename T::score_t score)
{
	requires std::convertible_to<typename T::SharedConstData::ParentClass, T>;
	{ data.nElements } noexcept -> std::convertible_to<size_t>;
	requires std::integral<typename T::score_t> || std::floating_point<typename T::score_t>;
	{T::ZERO} noexcept -> std::convertible_to<typename T::score_t const>;
	{T::EPSILON} noexcept -> std::convertible_to<typename T::score_t const>;
	{score + score} noexcept -> std::convertible_to<typename T::score_t const>;
	{T{dataPtr}} noexcept;
    {t.elementSetTaxonColor(index, index, index)} noexcept;
	{t.elementClearTaxonColor(index, index, index)} noexcept;
	{t.elementScore(index)} noexcept -> std::same_as<typename T::score_t>;
};

};
#endif
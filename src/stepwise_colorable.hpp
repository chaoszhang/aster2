#ifndef STEPWISE_COLORABLE_HPP
#define STEPWISE_COLORABLE_HPP

#include "common.hpp"

namespace stepwise_colorable{

template<class T> concept STEPWISE_COLORABLE = requires(T t, typename T::SharedConstData const *const dataPtr, typename T::SharedConstData data, std::size_t index, typename T::score_t score)
{
	requires std::convertible_to<typename T::SharedConstData::ParentClass, T>;
	{ data.nElements } noexcept -> std::convertible_to<size_t>;
	requires std::integral<typename T::score_t> || std::floating_point<typename T::score_t>;
	{T::ZERO} noexcept -> std::convertible_to<typename T::score_t const>;
	{T::EPSILON} noexcept -> std::convertible_to<typename T::score_t const>;
	{T{dataPtr}} noexcept;
    {t.elementSetTaxonColor(index, index, index)} noexcept;
	{t.elementClearTaxonColor(index, index, index)} noexcept;
	{t.elementScore(index)} noexcept -> std::same_as<typename T::score_t>;
};

ChangeLog logQUADRIPARTITION_STEPWISE_COLORABLE("QUADRIPARTITION_STEPWISE_COLORABLE",
	"2026-02-02", "Chao Zhang", "Supporting quadripartiton", "patch");

template<class T> concept QUADRIPARTITION_STEPWISE_COLORABLE = requires(T t, std::size_t index)
{
	requires STEPWISE_COLORABLE<T>;
	{ t.elementQuadripartitionScores(index) } noexcept -> std::convertible_to<std::array<typename T::score_t, 3> >; // ( lc+rc|sister+outgroup. lc+outgroup|sister+rc, rc+outgroup|sister+lc )
};

};
#endif
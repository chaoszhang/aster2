#ifndef STEPWISE_COLORABLE_HPP
#define STEPWISE_COLORABLE_HPP

#include "common.hpp"

namespace stepwise_colorable{

ChangeLog logSTEPWISE_COLORABLE("STEPWISE_COLORABLE",
	"2026-03-02", "Chao Zhang", "Signiture change of constructor and data.nElements()", "patch");

template<class T> concept STEPWISE_COLORABLE = requires(T t, typename T::SharedConstData const& data, std::size_t index, typename T::score_t score)
{
	requires std::convertible_to<typename T::SharedConstData::ParentClass, T>;
	{ data.nElements() } noexcept -> std::convertible_to<size_t>;
	requires std::integral<typename T::score_t> || std::floating_point<typename T::score_t>;
	{T::IS_ROOTED} noexcept -> std::convertible_to<bool>;
	{T::ZERO} noexcept -> std::convertible_to<typename T::score_t const>;
	{T::EPSILON} noexcept -> std::convertible_to<typename T::score_t const>;
	{T{data}} noexcept;
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

template<class T> concept ELEMENT_CLEAR_AND_SET_TAXON_COLOR = requires(T t, std::size_t index)
{
	requires STEPWISE_COLORABLE<T>;
	{ t.elementClearAndSetTaxonColor(index, index, index, index) } noexcept;
};

template<class T> concept SET_QUADRIPARTITION_MODE = requires(T t, bool b)
{
	requires STEPWISE_COLORABLE<T>;
	{ t.setQuadripartitionMode(b) } noexcept;
};

};
#endif
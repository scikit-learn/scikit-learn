/*
    pybind11/typing.h: Convenience wrapper classes for basic Python types
    with more explicit annotations.

    Copyright (c) 2023 Dustin Spicuzza <dustin@virtualroadside.com>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include "detail/common.h"
#include "cast.h"
#include "pytypes.h"

#include <algorithm>

#if defined(__cpp_nontype_template_args) && __cpp_nontype_template_args >= 201911L
#    define PYBIND11_TYPING_H_HAS_STRING_LITERAL
#    include <numeric>
#    include <ranges>
#    include <string_view>
#endif

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)
PYBIND11_NAMESPACE_BEGIN(typing)

/*
    The following types can be used to direct pybind11-generated docstrings
    to have have more explicit types (e.g., `list[str]` instead of `list`).
    Just use these in place of existing types.

    There is no additional enforcement of types at runtime.
*/

// Tuple type hint defined in cast.h for use in py::make_tuple to avoid circular includes

template <typename K, typename V>
class Dict : public dict {
    using dict::dict;
};

template <typename T>
class List : public list {
    using list::list;
};

template <typename T>
class Set : public set {
    using set::set;
};

template <typename T>
class Iterable : public iterable {
    using iterable::iterable;
};

template <typename T>
class Iterator : public iterator {
    using iterator::iterator;
};

template <typename Signature>
class Callable;

template <typename Return, typename... Args>
class Callable<Return(Args...)> : public function {
    using function::function;
};

template <typename T>
class Type : public type {
    using type::type;
};

template <typename... Types>
class Union : public object {
    PYBIND11_OBJECT_DEFAULT(Union, object, PyObject_Type)
    using object::object;
};

template <typename T>
class Optional : public object {
    PYBIND11_OBJECT_DEFAULT(Optional, object, PyObject_Type)
    using object::object;
};

template <typename T>
class Final : public object {
    PYBIND11_OBJECT_DEFAULT(Final, object, PyObject_Type)
    using object::object;
};

template <typename T>
class ClassVar : public object {
    PYBIND11_OBJECT_DEFAULT(ClassVar, object, PyObject_Type)
    using object::object;
};

template <typename T>
class TypeGuard : public bool_ {
    using bool_::bool_;
};

template <typename T>
class TypeIs : public bool_ {
    using bool_::bool_;
};

class NoReturn : public none {
    using none::none;
};

class Never : public none {
    using none::none;
};

#if defined(PYBIND11_TYPING_H_HAS_STRING_LITERAL)
template <size_t N>
struct StringLiteral {
    constexpr StringLiteral(const char (&str)[N]) { std::copy_n(str, N, name); }
    char name[N];
};

template <StringLiteral... StrLits>
class Literal : public object {
    PYBIND11_OBJECT_DEFAULT(Literal, object, PyObject_Type)
};

// Example syntax for creating a TypeVar.
// typedef typing::TypeVar<"T"> TypeVarT;
template <StringLiteral>
class TypeVar : public object {
    PYBIND11_OBJECT_DEFAULT(TypeVar, object, PyObject_Type)
    using object::object;
};
#endif

PYBIND11_NAMESPACE_END(typing)

PYBIND11_NAMESPACE_BEGIN(detail)

template <typename... Types>
struct handle_type_name<typing::Tuple<Types...>> {
    static constexpr auto name = const_name("tuple[")
                                 + ::pybind11::detail::concat(make_caster<Types>::name...)
                                 + const_name("]");
};

template <>
struct handle_type_name<typing::Tuple<>> {
    // PEP 484 specifies this syntax for an empty tuple
    static constexpr auto name = const_name("tuple[()]");
};

template <typename T>
struct handle_type_name<typing::Tuple<T, ellipsis>> {
    // PEP 484 specifies this syntax for a variable-length tuple
    static constexpr auto name
        = const_name("tuple[") + make_caster<T>::name + const_name(", ...]");
};

template <typename K, typename V>
struct handle_type_name<typing::Dict<K, V>> {
    static constexpr auto name = const_name("dict[") + make_caster<K>::name + const_name(", ")
                                 + make_caster<V>::name + const_name("]");
};

template <typename T>
struct handle_type_name<typing::List<T>> {
    static constexpr auto name = const_name("list[") + make_caster<T>::name + const_name("]");
};

template <typename T>
struct handle_type_name<typing::Set<T>> {
    static constexpr auto name = const_name("set[") + make_caster<T>::name + const_name("]");
};

template <typename T>
struct handle_type_name<typing::Iterable<T>> {
    static constexpr auto name
        = const_name("collections.abc.Iterable[") + make_caster<T>::name + const_name("]");
};

template <typename T>
struct handle_type_name<typing::Iterator<T>> {
    static constexpr auto name
        = const_name("collections.abc.Iterator[") + make_caster<T>::name + const_name("]");
};

template <typename Return, typename... Args>
struct handle_type_name<typing::Callable<Return(Args...)>> {
    using retval_type = conditional_t<std::is_same<Return, void>::value, void_type, Return>;
    static constexpr auto name
        = const_name("collections.abc.Callable[[")
          + ::pybind11::detail::concat(::pybind11::detail::arg_descr(make_caster<Args>::name)...)
          + const_name("], ") + ::pybind11::detail::return_descr(make_caster<retval_type>::name)
          + const_name("]");
};

template <typename Return>
struct handle_type_name<typing::Callable<Return(ellipsis)>> {
    // PEP 484 specifies this syntax for defining only return types of callables
    using retval_type = conditional_t<std::is_same<Return, void>::value, void_type, Return>;
    static constexpr auto name = const_name("collections.abc.Callable[..., ")
                                 + ::pybind11::detail::return_descr(make_caster<retval_type>::name)
                                 + const_name("]");
};

template <typename T>
struct handle_type_name<typing::Type<T>> {
    static constexpr auto name = const_name("type[") + make_caster<T>::name + const_name("]");
};

template <typename... Types>
struct handle_type_name<typing::Union<Types...>> {
    static constexpr auto name = ::pybind11::detail::union_concat(make_caster<Types>::name...);
};

template <typename T>
struct handle_type_name<typing::Optional<T>> {
    static constexpr auto name = make_caster<T>::name | make_caster<none>::name;
};

template <typename T>
struct handle_type_name<typing::Final<T>> {
    static constexpr auto name = const_name("typing.Final[")
                                 + ::pybind11::detail::return_descr(make_caster<T>::name)
                                 + const_name("]");
};

template <typename T>
struct handle_type_name<typing::ClassVar<T>> {
    static constexpr auto name
        = const_name("typing.ClassVar[") + make_caster<T>::name + const_name("]");
};

template <typename T>
struct handle_type_name<typing::TypeGuard<T>> {
    static constexpr auto name = const_name(PYBIND11_TYPE_GUARD_TYPE_HINT) + const_name("[")
                                 + make_caster<T>::name + const_name("]");
};

template <typename T>
struct handle_type_name<typing::TypeIs<T>> {
    static constexpr auto name = const_name(PYBIND11_TYPE_IS_TYPE_HINT) + const_name("[")
                                 + make_caster<T>::name + const_name("]");
};

template <>
struct handle_type_name<typing::NoReturn> {
    static constexpr auto name = const_name("typing.NoReturn");
};

template <>
struct handle_type_name<typing::Never> {
    static constexpr auto name = const_name(PYBIND11_NEVER_TYPE_HINT);
};

#if defined(PYBIND11_TYPING_H_HAS_STRING_LITERAL)
template <typing::StringLiteral StrLit>
consteval auto sanitize_string_literal() {
    constexpr std::string_view v(StrLit.name);
    constexpr std::string_view special_chars("!@%{}-");
    constexpr auto num_special_chars = std::accumulate(
        special_chars.begin(), special_chars.end(), (size_t) 0, [&v](auto acc, const char &c) {
            return std::move(acc) + std::ranges::count(v, c);
        });
    char result[v.size() + num_special_chars + 1];
    size_t i = 0;
    for (auto c : StrLit.name) {
        if (special_chars.find(c) != std::string_view::npos) {
            result[i++] = '!';
        }
        result[i++] = c;
    }
    return typing::StringLiteral(result);
}

template <typing::StringLiteral... Literals>
struct handle_type_name<typing::Literal<Literals...>> {
    static constexpr auto name
        = const_name("typing.Literal[")
          + pybind11::detail::concat(const_name(sanitize_string_literal<Literals>().name)...)
          + const_name("]");
};
template <typing::StringLiteral StrLit>
struct handle_type_name<typing::TypeVar<StrLit>> {
    static constexpr auto name = const_name(sanitize_string_literal<StrLit>().name);
};
#endif

PYBIND11_NAMESPACE_END(detail)
PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)

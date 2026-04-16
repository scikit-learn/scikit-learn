/*
    pybind11/detail/type_caster_base.h (originally first part of pybind11/cast.h)

    Copyright (c) 2016 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include <pybind11/gil.h>
#include <pybind11/pytypes.h>
#include <pybind11/trampoline_self_life_support.h>

#include "common.h"
#include "cpp_conduit.h"
#include "descr.h"
#include "dynamic_raw_ptr_cast_if_possible.h"
#include "internals.h"
#include "typeid.h"
#include "using_smart_holder.h"
#include "value_and_holder.h"

#include <cstdint>
#include <cstring>
#include <iterator>
#include <new>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>
#include <utility>
#include <vector>

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)
PYBIND11_NAMESPACE_BEGIN(detail)

/// A life support system for temporary objects created by `type_caster::load()`.
/// Adding a patient will keep it alive up until the enclosing function returns.
class loader_life_support {
private:
    // Thread-local top-of-stack for loader_life_support frames (linked via parent).
    // Observation: loader_life_support needs to be thread-local,
    // but we don't need to go to extra effort to keep it
    // per-interpreter (i.e., by putting it in internals) since
    // individual function calls are already isolated to a single
    // interpreter, even though they could potentially call into a
    // different interpreter later in the same call chain.  This
    // saves a significant cost per function call spent in
    // loader_life_support destruction.
    // Note for future C++17 simplification:
    // inline static thread_local loader_life_support *tls_current_frame = nullptr;
    static loader_life_support *&tls_current_frame() {
        static thread_local loader_life_support *frame_ptr = nullptr;
        return frame_ptr;
    }

    loader_life_support *parent = nullptr;
    std::unordered_set<PyObject *> keep_alive;

public:
    /// A new patient frame is created when a function is entered
    loader_life_support() {
        auto &frame = tls_current_frame();
        parent = frame;
        frame = this;
    }

    /// ... and destroyed after it returns
    ~loader_life_support() {
        auto &frame = tls_current_frame();
        if (frame != this) {
            pybind11_fail("loader_life_support: internal error");
        }
        frame = parent;
        for (auto *item : keep_alive) {
            Py_DECREF(item);
        }
    }

    /// This can only be used inside a pybind11-bound function, either by `argument_loader`
    /// at argument preparation time or by `py::cast()` at execution time.
    PYBIND11_NOINLINE static void add_patient(handle h) {
        loader_life_support *frame = tls_current_frame();
        if (!frame) {
            // NOTE: It would be nice to include the stack frames here, as this indicates
            // use of pybind11::cast<> outside the normal call framework, finding such
            // a location is challenging. Developers could consider printing out
            // stack frame addresses here using something like __builtin_frame_address(0)
            throw cast_error("When called outside a bound function, py::cast() cannot "
                             "do Python -> C++ conversions which require the creation "
                             "of temporary values");
        }

        if (frame->keep_alive.insert(h.ptr()).second) {
            Py_INCREF(h.ptr());
        }
    }
};

// Gets the cache entry for the given type, creating it if necessary.  The return value is the pair
// returned by emplace, i.e. an iterator for the entry and a bool set to `true` if the entry was
// just created.
inline std::pair<decltype(internals::registered_types_py)::iterator, bool>
all_type_info_get_cache(PyTypeObject *type);

// Band-aid workaround to fix a subtle but serious bug in a minimalistic fashion. See PR #4762.
inline void all_type_info_add_base_most_derived_first(std::vector<type_info *> &bases,
                                                      type_info *addl_base) {
    for (auto it = bases.begin(); it != bases.end(); it++) {
        type_info *existing_base = *it;
        if (PyType_IsSubtype(addl_base->type, existing_base->type) != 0) {
            bases.insert(it, addl_base);
            return;
        }
    }
    bases.push_back(addl_base);
}

// Populates a just-created cache entry.
PYBIND11_NOINLINE void all_type_info_populate(PyTypeObject *t, std::vector<type_info *> &bases) {
    assert(bases.empty());
    std::vector<PyTypeObject *> check;
    for (handle parent : reinterpret_borrow<tuple>(t->tp_bases)) {
        check.push_back(reinterpret_cast<PyTypeObject *>(parent.ptr()));
    }
    auto const &type_dict = get_internals().registered_types_py;
    for (size_t i = 0; i < check.size(); i++) {
        auto *type = check[i];
        // Ignore Python2 old-style class super type:
        if (!PyType_Check((PyObject *) type)) {
            continue;
        }

        // Check `type` in the current set of registered python types:
        auto it = type_dict.find(type);
        if (it != type_dict.end()) {
            // We found a cache entry for it, so it's either pybind-registered or has pre-computed
            // pybind bases, but we have to make sure we haven't already seen the type(s) before:
            // we want to follow Python/virtual C++ rules that there should only be one instance of
            // a common base.
            for (auto *tinfo : it->second) {
                // NB: Could use a second set here, rather than doing a linear search, but since
                // having a large number of immediate pybind11-registered types seems fairly
                // unlikely, that probably isn't worthwhile.
                bool found = false;
                for (auto *known : bases) {
                    if (known == tinfo) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    all_type_info_add_base_most_derived_first(bases, tinfo);
                }
            }
        } else if (type->tp_bases) {
            // It's some python type, so keep follow its bases classes to look for one or more
            // registered types
            if (i + 1 == check.size()) {
                // When we're at the end, we can pop off the current element to avoid growing
                // `check` when adding just one base (which is typical--i.e. when there is no
                // multiple inheritance)
                check.pop_back();
                i--;
            }
            for (handle parent : reinterpret_borrow<tuple>(type->tp_bases)) {
                check.push_back(reinterpret_cast<PyTypeObject *>(parent.ptr()));
            }
        }
    }
}

/**
 * Extracts vector of type_info pointers of pybind-registered roots of the given Python type.  Will
 * be just 1 pybind type for the Python type of a pybind-registered class, or for any Python-side
 * derived class that uses single inheritance.  Will contain as many types as required for a Python
 * class that uses multiple inheritance to inherit (directly or indirectly) from multiple
 * pybind-registered classes.  Will be empty if neither the type nor any base classes are
 * pybind-registered.
 *
 * The value is cached for the lifetime of the Python type.
 */
inline const std::vector<detail::type_info *> &all_type_info(PyTypeObject *type) {
    return all_type_info_get_cache(type).first->second;
}

/**
 * Gets a single pybind11 type info for a python type.  Returns nullptr if neither the type nor any
 * ancestors are pybind11-registered.  Throws an exception if there are multiple bases--use
 * `all_type_info` instead if you want to support multiple bases.
 */
PYBIND11_NOINLINE detail::type_info *get_type_info(PyTypeObject *type) {
    const auto &bases = all_type_info(type);
    if (bases.empty()) {
        return nullptr;
    }
    if (bases.size() > 1) {
        pybind11_fail(
            "pybind11::detail::get_type_info: type has multiple pybind11-registered bases");
    }
    return bases.front();
}

inline detail::type_info *get_local_type_info_lock_held(const std::type_info &tp) {
    const auto &locals = get_local_internals().registered_types_cpp;
    auto it = locals.find(&tp);
    if (it != locals.end()) {
        return it->second;
    }
    return nullptr;
}

inline detail::type_info *get_local_type_info(const std::type_info &tp) {
    // NB: internals and local_internals share a single mutex
    PYBIND11_LOCK_INTERNALS(get_internals());
    return get_local_type_info_lock_held(tp);
}

inline detail::type_info *get_global_type_info_lock_held(const std::type_info &tp) {
    // This is a two-level lookup. Hopefully we find the type info in
    // registered_types_cpp_fast, but if not we try
    // registered_types_cpp and fill registered_types_cpp_fast for
    // next time.
    detail::type_info *type_info = nullptr;
    auto &internals = get_internals();
#if PYBIND11_INTERNALS_VERSION >= 12
    auto &fast_types = internals.registered_types_cpp_fast;
#endif
    auto &types = internals.registered_types_cpp;
#if PYBIND11_INTERNALS_VERSION >= 12
    auto fast_it = fast_types.find(&tp);
    if (fast_it != fast_types.end()) {
#    ifndef NDEBUG
        auto types_it = types.find(std::type_index(tp));
        assert(types_it != types.end());
        assert(types_it->second == fast_it->second);
#    endif
        return fast_it->second;
    }
#endif // PYBIND11_INTERNALS_VERSION >= 12

    auto it = types.find(std::type_index(tp));
    if (it != types.end()) {
#if PYBIND11_INTERNALS_VERSION >= 12
        // We found the type in the slow map but not the fast one, so
        // some other DSO added it (otherwise it would be in the fast
        // map under &tp) and therefore we must be an alias. Record
        // that.
        it->second->alias_chain.push_front(&tp);
        fast_types.emplace(&tp, it->second);
#endif
        type_info = it->second;
    }
    return type_info;
}

inline detail::type_info *get_global_type_info(const std::type_info &tp) {
    PYBIND11_LOCK_INTERNALS(get_internals());
    return get_global_type_info_lock_held(tp);
}

/// Return the type info for a given C++ type; on lookup failure can either throw or return
/// nullptr.
PYBIND11_NOINLINE detail::type_info *get_type_info(const std::type_info &tp,
                                                   bool throw_if_missing = false) {
    PYBIND11_LOCK_INTERNALS(get_internals());
    if (auto *ltype = get_local_type_info_lock_held(tp)) {
        return ltype;
    }
    if (auto *gtype = get_global_type_info_lock_held(tp)) {
        return gtype;
    }

    if (throw_if_missing) {
        std::string tname = tp.name();
        detail::clean_type_id(tname);
        pybind11_fail("pybind11::detail::get_type_info: unable to find type info for \""
                      + std::move(tname) + '"');
    }
    return nullptr;
}

PYBIND11_NOINLINE handle get_type_handle(const std::type_info &tp, bool throw_if_missing) {
    detail::type_info *type_info = get_type_info(tp, throw_if_missing);
    return handle(type_info ? (reinterpret_cast<PyObject *>(type_info->type)) : nullptr);
}

inline bool try_incref(PyObject *obj) {
    // Tries to increment the reference count of an object if it's not zero.
#if defined(Py_GIL_DISABLED) && PY_VERSION_HEX >= 0x030E00A4
    return PyUnstable_TryIncRef(obj);
#elif defined(Py_GIL_DISABLED)
    // See
    // https://github.com/python/cpython/blob/d05140f9f77d7dfc753dd1e5ac3a5962aaa03eff/Include/internal/pycore_object.h#L761
    uint32_t local = _Py_atomic_load_uint32_relaxed(&obj->ob_ref_local);
    local += 1;
    if (local == 0) {
        // immortal
        return true;
    }
    if (_Py_IsOwnedByCurrentThread(obj)) {
        _Py_atomic_store_uint32_relaxed(&obj->ob_ref_local, local);
#    ifdef Py_REF_DEBUG
        _Py_INCREF_IncRefTotal();
#    endif
        return true;
    }
    Py_ssize_t shared = _Py_atomic_load_ssize_relaxed(&obj->ob_ref_shared);
    for (;;) {
        // If the shared refcount is zero and the object is either merged
        // or may not have weak references, then we cannot incref it.
        if (shared == 0 || shared == _Py_REF_MERGED) {
            return false;
        }

        if (_Py_atomic_compare_exchange_ssize(
                &obj->ob_ref_shared, &shared, shared + (1 << _Py_REF_SHARED_SHIFT))) {
#    ifdef Py_REF_DEBUG
            _Py_INCREF_IncRefTotal();
#    endif
            return true;
        }
    }
#else
    assert(Py_REFCNT(obj) > 0);
    Py_INCREF(obj);
    return true;
#endif
}

// Searches the inheritance graph for a registered Python instance, using all_type_info().
PYBIND11_NOINLINE handle find_registered_python_instance(void *src,
                                                         const detail::type_info *tinfo) {
    return with_instance_map(src, [&](instance_map &instances) {
        auto it_instances = instances.equal_range(src);
        for (auto it_i = it_instances.first; it_i != it_instances.second; ++it_i) {
            for (auto *instance_type : detail::all_type_info(Py_TYPE(it_i->second))) {
                if (instance_type && same_type(*instance_type->cpptype, *tinfo->cpptype)) {
                    auto *wrapper = reinterpret_cast<PyObject *>(it_i->second);
                    if (try_incref(wrapper)) {
                        return handle(wrapper);
                    }
                }
            }
        }
        return handle();
    });
}

// Container for accessing and iterating over an instance's values/holders
struct values_and_holders {
private:
    instance *inst;
    using type_vec = std::vector<detail::type_info *>;
    const type_vec &tinfo;

public:
    explicit values_and_holders(instance *inst)
        : inst{inst}, tinfo(all_type_info(Py_TYPE(inst))) {}

    explicit values_and_holders(PyObject *obj)
        : inst{nullptr}, tinfo(all_type_info(Py_TYPE(obj))) {
        if (!tinfo.empty()) {
            inst = reinterpret_cast<instance *>(obj);
        }
    }

    struct iterator {
    private:
        instance *inst = nullptr;
        const type_vec *types = nullptr;
        value_and_holder curr;
        friend struct values_and_holders;
        iterator(instance *inst, const type_vec *tinfo) : inst{inst}, types{tinfo} {
            if (inst != nullptr) {
                assert(!types->empty());
                curr = value_and_holder(
                    inst /* instance */,
                    (*types)[0] /* type info */,
                    0, /* vpos: (non-simple types only): the first vptr comes first */
                    0 /* index */);
            }
        }
        // Past-the-end iterator:
        explicit iterator(size_t end) : curr(end) {}

    public:
        bool operator==(const iterator &other) const { return curr.index == other.curr.index; }
        bool operator!=(const iterator &other) const { return curr.index != other.curr.index; }
        iterator &operator++() {
            if (!inst->simple_layout) {
                curr.vh += 1 + (*types)[curr.index]->holder_size_in_ptrs;
            }
            ++curr.index;
            curr.type = curr.index < types->size() ? (*types)[curr.index] : nullptr;
            return *this;
        }
        value_and_holder &operator*() { return curr; }
        value_and_holder *operator->() { return &curr; }
    };

    iterator begin() { return iterator(inst, &tinfo); }
    iterator end() { return iterator(tinfo.size()); }

    iterator find(const type_info *find_type) {
        auto it = begin(), endit = end();
        while (it != endit && it->type != find_type) {
            ++it;
        }
        return it;
    }

    size_t size() { return tinfo.size(); }

    // Band-aid workaround to fix a subtle but serious bug in a minimalistic fashion. See PR #4762.
    bool is_redundant_value_and_holder(const value_and_holder &vh) {
        for (size_t i = 0; i < vh.index; i++) {
            if (PyType_IsSubtype(tinfo[i]->type, tinfo[vh.index]->type) != 0) {
                return true;
            }
        }
        return false;
    }
};

/**
 * Extracts C++ value and holder pointer references from an instance (which may contain multiple
 * values/holders for python-side multiple inheritance) that match the given type.  Throws an error
 * if the given type (or ValueType, if omitted) is not a pybind11 base of the given instance.  If
 * `find_type` is omitted (or explicitly specified as nullptr) the first value/holder are returned,
 * regardless of type (and the resulting .type will be nullptr).
 *
 * The returned object should be short-lived: in particular, it must not outlive the called-upon
 * instance.
 */
PYBIND11_NOINLINE value_and_holder
instance::get_value_and_holder(const type_info *find_type /*= nullptr default in common.h*/,
                               bool throw_if_missing /*= true in common.h*/) {
    // Optimize common case:
    if (!find_type || Py_TYPE(this) == find_type->type) {
        return value_and_holder(this, find_type, 0, 0);
    }

    detail::values_and_holders vhs(this);
    auto it = vhs.find(find_type);
    if (it != vhs.end()) {
        return *it;
    }

    if (!throw_if_missing) {
        return value_and_holder();
    }

#if defined(PYBIND11_DETAILED_ERROR_MESSAGES)
    pybind11_fail("pybind11::detail::instance::get_value_and_holder: `"
                  + get_fully_qualified_tp_name(find_type->type)
                  + "' is not a pybind11 base of the given `"
                  + get_fully_qualified_tp_name(Py_TYPE(this)) + "' instance");
#else
    pybind11_fail(
        "pybind11::detail::instance::get_value_and_holder: "
        "type is not a pybind11 base of the given instance "
        "(#define PYBIND11_DETAILED_ERROR_MESSAGES or compile in debug mode for type details)");
#endif
}

PYBIND11_NOINLINE void instance::allocate_layout() {
    const auto &tinfo = all_type_info(Py_TYPE(this));

    const size_t n_types = tinfo.size();

    if (n_types == 0) {
        pybind11_fail(
            "instance allocation failed: new instance has no pybind11-registered base types");
    }

    simple_layout
        = n_types == 1 && tinfo.front()->holder_size_in_ptrs <= instance_simple_holder_in_ptrs();

    // Simple path: no python-side multiple inheritance, and a small-enough holder
    if (simple_layout) {
        simple_value_holder[0] = nullptr;
        simple_holder_constructed = false;
        simple_instance_registered = false;
    } else { // multiple base types or a too-large holder
        // Allocate space to hold: [v1*][h1][v2*][h2]...[bb...] where [vN*] is a value pointer,
        // [hN] is the (uninitialized) holder instance for value N, and [bb...] is a set of bool
        // values that tracks whether each associated holder has been initialized.  Each [block] is
        // padded, if necessary, to an integer multiple of sizeof(void *).
        size_t space = 0;
        for (auto *t : tinfo) {
            space += 1;                      // value pointer
            space += t->holder_size_in_ptrs; // holder instance
        }
        size_t flags_at = space;
        space += size_in_ptrs(n_types); // status bytes (holder_constructed and
                                        // instance_registered)

        // Allocate space for flags, values, and holders, and initialize it to 0 (flags and values,
        // in particular, need to be 0).  Use Python's memory allocation
        // functions: Python is using pymalloc, which is designed to be
        // efficient for small allocations like the one we're doing here;
        // for larger allocations they are just wrappers around malloc.
        // TODO: is this still true for pure Python 3.6?
        nonsimple.values_and_holders = static_cast<void **>(PyMem_Calloc(space, sizeof(void *)));
        if (!nonsimple.values_and_holders) {
            throw std::bad_alloc();
        }
        nonsimple.status
            = reinterpret_cast<std::uint8_t *>(&nonsimple.values_and_holders[flags_at]);
    }
    owned = true;
}

// NOLINTNEXTLINE(readability-make-member-function-const)
PYBIND11_NOINLINE void instance::deallocate_layout() {
    if (!simple_layout) {
        PyMem_Free(reinterpret_cast<void *>(nonsimple.values_and_holders));
    }
}

PYBIND11_NOINLINE bool isinstance_generic(handle obj, const std::type_info &tp) {
    handle type = detail::get_type_handle(tp, false);
    if (!type) {
        return false;
    }
    return isinstance(obj, type);
}

PYBIND11_NOINLINE handle get_object_handle(const void *ptr, const detail::type_info *type) {
    return with_instance_map(ptr, [&](instance_map &instances) {
        auto range = instances.equal_range(ptr);
        for (auto it = range.first; it != range.second; ++it) {
            for (const auto &vh : values_and_holders(it->second)) {
                if (vh.type == type) {
                    return handle(reinterpret_cast<PyObject *>(it->second));
                }
            }
        }
        return handle();
    });
}

// Information about how type_caster_generic::cast() can obtain its source object
struct cast_sources {
    // A type-erased pointer and the type it points to
    struct raw_source {
        const void *cppobj;
        const std::type_info *cpptype;
    };

    // A C++ pointer and the Python type info we will convert it to;
    // we expect that cppobj points to something of type tinfo->cpptype
    struct resolved_source {
        const void *cppobj;
        const type_info *tinfo;
    };

    // Use the given pointer with its compile-time type, possibly downcast
    // via polymorphic_type_hook()
    template <typename itype>
    explicit cast_sources(const itype *ptr);

    // Use the given pointer and type
    // NOLINTNEXTLINE(google-explicit-constructor)
    cast_sources(const raw_source &orig) : original(orig) { result = resolve(); }

    // Use the given object and pybind11 type info. NB: if tinfo is null,
    // this does not provide enough information to use a foreign type or
    // to render a useful error message
    cast_sources(const void *obj, const detail::type_info *tinfo)
        : original{obj, tinfo ? tinfo->cpptype : nullptr}, result{obj, tinfo} {}

    // The object passed to cast(), with its static type.
    // original.type must not be null if resolve() will be called.
    // original.obj may be null if we're converting nullptr to a Python None
    raw_source original;

    // A more-derived version of `original` provided by a
    // polymorphic_type_hook. downcast.type may be null if this is not
    // a relevant concept for the current cast.
    raw_source downcast{};

    // The source to use for this cast, and the corresponding pybind11
    // type_info. If the type_info is null, then pybind11 doesn't know
    // about this type.
    resolved_source result;

    // Returns true if the cast will use a pybind11 type that uses
    // a smart holder.
    bool creates_smart_holder() const {
        return result.tinfo != nullptr
               && result.tinfo->holder_enum_v == detail::holder_enum_t::smart_holder;
    }

private:
    resolved_source resolve() {
        if (downcast.cpptype) {
            if (same_type(*original.cpptype, *downcast.cpptype)) {
                downcast.cpptype = nullptr;
            } else if (const auto *tpi = get_type_info(*downcast.cpptype)) {
                return {downcast.cppobj, tpi};
            }
        }
        if (const auto *tpi = get_type_info(*original.cpptype)) {
            return {original.cppobj, tpi};
        }
        return {nullptr, nullptr};
    }
};

// Forward declarations
void keep_alive_impl(handle nurse, handle patient);
inline PyObject *make_new_instance(PyTypeObject *type);

PYBIND11_WARNING_PUSH
PYBIND11_WARNING_DISABLE_GCC("-Wredundant-decls")

// PYBIND11:REMINDER: Needs refactoring of existing pybind11 code.
inline bool deregister_instance(instance *self, void *valptr, const type_info *tinfo);

PYBIND11_WARNING_POP

PYBIND11_NAMESPACE_BEGIN(smart_holder_type_caster_support)

struct value_and_holder_helper {
    value_and_holder loaded_v_h;

    bool have_holder() const {
        return loaded_v_h.vh != nullptr && loaded_v_h.holder_constructed();
    }

    smart_holder &holder() const { return loaded_v_h.holder<smart_holder>(); }

    void throw_if_uninitialized_or_disowned_holder(const char *typeid_name) const {
        static const std::string missing_value_msg = "Missing value for wrapped C++ type `";
        if (!holder().is_populated) {
            throw value_error(missing_value_msg + clean_type_id(typeid_name)
                              + "`: Python instance is uninitialized.");
        }
        if (!holder().has_pointee()) {
            throw value_error(missing_value_msg + clean_type_id(typeid_name)
                              + "`: Python instance was disowned.");
        }
    }

    void throw_if_uninitialized_or_disowned_holder(const std::type_info &type_info) const {
        throw_if_uninitialized_or_disowned_holder(type_info.name());
    }

    // have_holder() must be true or this function will fail.
    void throw_if_instance_is_currently_owned_by_shared_ptr(const type_info *tinfo) const {
        auto *vptr_gd_ptr = tinfo->get_memory_guarded_delete(holder().vptr);
        if (vptr_gd_ptr != nullptr && !vptr_gd_ptr->released_ptr.expired()) {
            throw value_error("Python instance is currently owned by a std::shared_ptr.");
        }
    }

    void *get_void_ptr_or_nullptr() const {
        if (have_holder()) {
            auto &hld = holder();
            if (hld.is_populated && hld.has_pointee()) {
                return hld.template as_raw_ptr_unowned<void>();
            }
        }
        return nullptr;
    }
};

template <typename T, typename D>
handle smart_holder_from_unique_ptr(std::unique_ptr<T, D> &&src,
                                    return_value_policy policy,
                                    handle parent,
                                    const cast_sources::resolved_source &cs) {
    if (policy == return_value_policy::copy) {
        throw cast_error("return_value_policy::copy is invalid for unique_ptr.");
    }
    if (!src) {
        return none().release();
    }
    // cs.cppobj is the subobject pointer appropriate for tinfo (may differ from src.get()
    // under MI/VI). Use this for Python identity/registration, but keep ownership on T*.
    void *src_raw_void_ptr = const_cast<void *>(cs.cppobj);
    assert(cs.tinfo != nullptr);
    const detail::type_info *tinfo = cs.tinfo;
    if (handle existing_inst = find_registered_python_instance(src_raw_void_ptr, tinfo)) {
        auto *self_life_support = tinfo->get_trampoline_self_life_support(src.get());
        if (self_life_support != nullptr) {
            value_and_holder &v_h = self_life_support->v_h;
            if (v_h.inst != nullptr && v_h.vh != nullptr) {
                auto &holder = v_h.holder<smart_holder>();
                if (!holder.is_disowned) {
                    pybind11_fail("smart_holder_from_unique_ptr: unexpected "
                                  "smart_holder.is_disowned failure.");
                }
                // Critical transfer-of-ownership section. This must stay together.
                self_life_support->deactivate_life_support();
                holder.reclaim_disowned(tinfo->get_memory_guarded_delete);
                (void) src.release();
                // Critical section end.
                return existing_inst;
            }
        }
        throw cast_error("Invalid unique_ptr: another instance owns this pointer already.");
    }

    auto inst = reinterpret_steal<object>(make_new_instance(tinfo->type));
    auto *inst_raw_ptr = reinterpret_cast<instance *>(inst.ptr());
    inst_raw_ptr->owned = true;
    void *&valueptr = values_and_holders(inst_raw_ptr).begin()->value_ptr();
    valueptr = src_raw_void_ptr;

    if (static_cast<void *>(src.get()) == src_raw_void_ptr) {
        // This is a multiple-inheritance situation that is incompatible with the current
        // shared_from_this handling (see PR #3023). Is there a better solution?
        src_raw_void_ptr = nullptr;
    }
    auto smhldr = smart_holder::from_unique_ptr(std::move(src), src_raw_void_ptr);
    tinfo->init_instance(inst_raw_ptr, static_cast<const void *>(&smhldr));

    if (policy == return_value_policy::reference_internal) {
        keep_alive_impl(inst, parent);
    }

    return inst.release();
}

template <typename T, typename D>
handle smart_holder_from_unique_ptr(std::unique_ptr<T const, D> &&src,
                                    return_value_policy policy,
                                    handle parent,
                                    const cast_sources::resolved_source &cs) {
    return smart_holder_from_unique_ptr(
        std::unique_ptr<T, D>(const_cast<T *>(src.release()),
                              std::move(src.get_deleter())), // Const2Mutbl
        policy,
        parent,
        cs);
}

template <typename T>
handle smart_holder_from_shared_ptr(const std::shared_ptr<T> &src,
                                    return_value_policy policy,
                                    handle parent,
                                    const cast_sources::resolved_source &cs) {
    switch (policy) {
        case return_value_policy::automatic:
        case return_value_policy::automatic_reference:
            break;
        case return_value_policy::take_ownership:
            throw cast_error("Invalid return_value_policy for shared_ptr (take_ownership).");
        case return_value_policy::copy:
        case return_value_policy::move:
            break;
        case return_value_policy::reference:
            throw cast_error("Invalid return_value_policy for shared_ptr (reference).");
        case return_value_policy::reference_internal:
            break;
    }
    if (!src) {
        return none().release();
    }

    // cs.cppobj is the subobject pointer appropriate for tinfo (may differ from src.get()
    // under MI/VI). Use this for Python identity/registration, but keep ownership on T*.
    void *src_raw_void_ptr = const_cast<void *>(cs.cppobj);
    assert(cs.tinfo != nullptr);
    const detail::type_info *tinfo = cs.tinfo;
    if (handle existing_inst = find_registered_python_instance(src_raw_void_ptr, tinfo)) {
        // PYBIND11:REMINDER: MISSING: Enforcement of consistency with existing smart_holder.
        // PYBIND11:REMINDER: MISSING: keep_alive.
        return existing_inst;
    }

    auto inst = reinterpret_steal<object>(make_new_instance(tinfo->type));
    auto *inst_raw_ptr = reinterpret_cast<instance *>(inst.ptr());
    inst_raw_ptr->owned = true;
    void *&valueptr = values_and_holders(inst_raw_ptr).begin()->value_ptr();
    valueptr = src_raw_void_ptr;

    auto smhldr = smart_holder::from_shared_ptr(std::shared_ptr<void>(src, src_raw_void_ptr));
    tinfo->init_instance(inst_raw_ptr, static_cast<const void *>(&smhldr));

    if (policy == return_value_policy::reference_internal) {
        keep_alive_impl(inst, parent);
    }

    return inst.release();
}

template <typename T>
handle smart_holder_from_shared_ptr(const std::shared_ptr<T const> &src,
                                    return_value_policy policy,
                                    handle parent,
                                    const cast_sources::resolved_source &cs) {
    return smart_holder_from_shared_ptr(std::const_pointer_cast<T>(src), // Const2Mutbl
                                        policy,
                                        parent,
                                        cs);
}

struct shared_ptr_parent_life_support {
    PyObject *parent;
    explicit shared_ptr_parent_life_support(PyObject *parent) : parent{parent} {
        Py_INCREF(parent);
    }
    // NOLINTNEXTLINE(readability-make-member-function-const)
    void operator()(void *) {
        gil_scoped_acquire gil;
        Py_DECREF(parent);
    }
};

struct shared_ptr_trampoline_self_life_support {
    PyObject *self;
    explicit shared_ptr_trampoline_self_life_support(instance *inst)
        : self{reinterpret_cast<PyObject *>(inst)} {
        gil_scoped_acquire gil;
        Py_INCREF(self);
    }
    // NOLINTNEXTLINE(readability-make-member-function-const)
    void operator()(void *) {
        gil_scoped_acquire gil;
        Py_DECREF(self);
    }
};

template <typename T,
          typename D,
          typename std::enable_if<std::is_default_constructible<D>::value, int>::type = 0>
inline std::unique_ptr<T, D> unique_with_deleter(T *raw_ptr, std::unique_ptr<D> &&deleter) {
    if (deleter == nullptr) {
        return std::unique_ptr<T, D>(raw_ptr);
    }
    return std::unique_ptr<T, D>(raw_ptr, std::move(*deleter));
}

template <typename T,
          typename D,
          typename std::enable_if<!std::is_default_constructible<D>::value, int>::type = 0>
inline std::unique_ptr<T, D> unique_with_deleter(T *raw_ptr, std::unique_ptr<D> &&deleter) {
    if (deleter == nullptr) {
        pybind11_fail("smart_holder_type_casters: deleter is not default constructible and no"
                      " instance available to return.");
    }
    return std::unique_ptr<T, D>(raw_ptr, std::move(*deleter));
}

template <typename T>
struct load_helper : value_and_holder_helper {
    bool was_populated = false;
    bool python_instance_is_alias = false;

    void maybe_set_python_instance_is_alias(handle src) {
        if (was_populated) {
            python_instance_is_alias = reinterpret_cast<instance *>(src.ptr())->is_alias;
        }
    }

    static std::shared_ptr<T> make_shared_ptr_with_responsible_parent(T *raw_ptr, handle parent) {
        return std::shared_ptr<T>(raw_ptr, shared_ptr_parent_life_support(parent.ptr()));
    }

    std::shared_ptr<T> load_as_shared_ptr(const type_info *tinfo,
                                          void *void_raw_ptr,
                                          handle responsible_parent = nullptr,
                                          // to support py::potentially_slicing_weak_ptr
                                          // with minimal added code complexity:
                                          bool force_potentially_slicing_shared_ptr
                                          = false) const {
        if (!have_holder()) {
            return nullptr;
        }
        throw_if_uninitialized_or_disowned_holder(typeid(T));
        smart_holder &hld = holder();
        hld.ensure_is_not_disowned("load_as_shared_ptr");
        if (hld.vptr_is_using_noop_deleter) {
            if (responsible_parent) {
                return make_shared_ptr_with_responsible_parent(static_cast<T *>(void_raw_ptr),
                                                               responsible_parent);
            }
            throw std::runtime_error("Non-owning holder (load_as_shared_ptr).");
        }
        auto *type_raw_ptr = static_cast<T *>(void_raw_ptr);
        if (python_instance_is_alias && !force_potentially_slicing_shared_ptr) {
            auto *vptr_gd_ptr = tinfo->get_memory_guarded_delete(holder().vptr);
            if (vptr_gd_ptr != nullptr) {
                std::shared_ptr<void> released_ptr = vptr_gd_ptr->released_ptr.lock();
                if (released_ptr) {
                    return std::shared_ptr<T>(released_ptr, type_raw_ptr);
                }
                std::shared_ptr<T> to_be_released(
                    type_raw_ptr, shared_ptr_trampoline_self_life_support(loaded_v_h.inst));
                vptr_gd_ptr->released_ptr = to_be_released;
                return to_be_released;
            }
            auto *sptsls_ptr = std::get_deleter<shared_ptr_trampoline_self_life_support>(hld.vptr);
            if (sptsls_ptr != nullptr) {
                // This code is reachable only if there are multiple registered_instances for the
                // same pointee.
                if (reinterpret_cast<PyObject *>(loaded_v_h.inst) == sptsls_ptr->self) {
                    pybind11_fail("smart_holder_type_caster_support load_as_shared_ptr failure: "
                                  "loaded_v_h.inst == sptsls_ptr->self");
                }
            }
            if (sptsls_ptr != nullptr || !memory::type_has_shared_from_this(type_raw_ptr)) {
                return std::shared_ptr<T>(
                    type_raw_ptr, shared_ptr_trampoline_self_life_support(loaded_v_h.inst));
            }
            if (hld.vptr_is_external_shared_ptr) {
                pybind11_fail("smart_holder_type_casters load_as_shared_ptr failure: not "
                              "implemented: trampoline-self-life-support for external shared_ptr "
                              "to type inheriting from std::enable_shared_from_this.");
            }
            pybind11_fail(
                "smart_holder_type_casters: load_as_shared_ptr failure: internal inconsistency.");
        }
        std::shared_ptr<void> void_shd_ptr = hld.template as_shared_ptr<void>();
        return std::shared_ptr<T>(void_shd_ptr, type_raw_ptr);
    }

    template <typename D>
    std::unique_ptr<T, D> load_as_unique_ptr(const type_info *tinfo,
                                             void *raw_void_ptr,
                                             const char *context = "load_as_unique_ptr") {
        if (!have_holder()) {
            return unique_with_deleter<T, D>(nullptr, std::unique_ptr<D>());
        }
        throw_if_uninitialized_or_disowned_holder(typeid(T));
        throw_if_instance_is_currently_owned_by_shared_ptr(tinfo);
        holder().ensure_is_not_disowned(context);
        holder().template ensure_compatible_uqp_del<T, D>(context);
        holder().ensure_use_count_1(context);

        T *raw_type_ptr = static_cast<T *>(raw_void_ptr);

        auto *self_life_support = tinfo->get_trampoline_self_life_support(raw_type_ptr);
        // This is enforced indirectly by a static_assert in the class_ implementation:
        assert(!python_instance_is_alias || self_life_support);

        std::unique_ptr<D> extracted_deleter
            = holder().template extract_deleter<T, D>(context, tinfo->get_memory_guarded_delete);

        // Critical transfer-of-ownership section. This must stay together.
        if (self_life_support != nullptr) {
            holder().disown(tinfo->get_memory_guarded_delete);
        } else {
            holder().release_ownership(tinfo->get_memory_guarded_delete);
        }
        auto result = unique_with_deleter<T, D>(raw_type_ptr, std::move(extracted_deleter));
        if (self_life_support != nullptr) {
            self_life_support->activate_life_support(loaded_v_h);
        } else {
            void *value_void_ptr = loaded_v_h.value_ptr();
            loaded_v_h.value_ptr() = nullptr;
            deregister_instance(loaded_v_h.inst, value_void_ptr, loaded_v_h.type);
        }
        // Critical section end.

        return result;
    }

    // This assumes load_as_shared_ptr succeeded(), and the returned shared_ptr is still alive.
    // The returned unique_ptr is meant to never expire (the behavior is undefined otherwise).
    template <typename D>
    std::unique_ptr<T, D> load_as_const_unique_ptr(const type_info *tinfo,
                                                   T *raw_type_ptr,
                                                   const char *context
                                                   = "load_as_const_unique_ptr") {
        if (!have_holder()) {
            return unique_with_deleter<T, D>(nullptr, std::unique_ptr<D>());
        }
        holder().template ensure_compatible_uqp_del<T, D>(context);
        return unique_with_deleter<T, D>(raw_type_ptr,
                                         std::move(holder().template extract_deleter<T, D>(
                                             context, tinfo->get_memory_guarded_delete)));
    }
};

PYBIND11_NAMESPACE_END(smart_holder_type_caster_support)

class type_caster_generic {
public:
    PYBIND11_NOINLINE explicit type_caster_generic(const std::type_info &type_info)
        : typeinfo(get_type_info(type_info)), cpptype(&type_info) {}

    explicit type_caster_generic(const type_info *typeinfo)
        : typeinfo(typeinfo), cpptype(typeinfo ? typeinfo->cpptype : nullptr) {}

    bool load(handle src, bool convert) { return load_impl<type_caster_generic>(src, convert); }

    static handle cast(const void *src,
                       return_value_policy policy,
                       handle parent,
                       const detail::type_info *tinfo,
                       void *(*copy_constructor)(const void *),
                       void *(*move_constructor)(const void *),
                       const void *existing_holder = nullptr) {
        cast_sources srcs{src, tinfo};
        return cast(srcs, policy, parent, copy_constructor, move_constructor, existing_holder);
    }

    PYBIND11_NOINLINE static handle cast(const cast_sources &srcs,
                                         return_value_policy policy,
                                         handle parent,
                                         void *(*copy_constructor)(const void *),
                                         void *(*move_constructor)(const void *),
                                         const void *existing_holder = nullptr) {
        if (!srcs.result.tinfo) {
            // No pybind11 type info. Raise an exception.
            std::string tname = srcs.downcast.cpptype   ? srcs.downcast.cpptype->name()
                                : srcs.original.cpptype ? srcs.original.cpptype->name()
                                                        : "<unspecified>";
            detail::clean_type_id(tname);
            std::string msg = "Unregistered type : " + tname;
            set_error(PyExc_TypeError, msg.c_str());
            return handle();
        }

        void *src = const_cast<void *>(srcs.result.cppobj);
        if (src == nullptr) {
            return none().release();
        }
        const type_info *tinfo = srcs.result.tinfo;

        if (handle registered_inst = find_registered_python_instance(src, tinfo)) {
            return registered_inst;
        }

        auto inst = reinterpret_steal<object>(make_new_instance(tinfo->type));
        auto *wrapper = reinterpret_cast<instance *>(inst.ptr());
        wrapper->owned = false;
        void *&valueptr = values_and_holders(wrapper).begin()->value_ptr();

        switch (policy) {
            case return_value_policy::automatic:
            case return_value_policy::take_ownership:
                valueptr = src;
                wrapper->owned = true;
                break;

            case return_value_policy::automatic_reference:
            case return_value_policy::reference:
                valueptr = src;
                wrapper->owned = false;
                break;

            case return_value_policy::copy:
                if (copy_constructor) {
                    valueptr = copy_constructor(src);
                } else {
#if defined(PYBIND11_DETAILED_ERROR_MESSAGES)
                    std::string type_name(tinfo->cpptype->name());
                    detail::clean_type_id(type_name);
                    throw cast_error("return_value_policy = copy, but type " + type_name
                                     + " is non-copyable!");
#else
                    throw cast_error("return_value_policy = copy, but type is "
                                     "non-copyable! (#define PYBIND11_DETAILED_ERROR_MESSAGES or "
                                     "compile in debug mode for details)");
#endif
                }
                wrapper->owned = true;
                break;

            case return_value_policy::move:
                if (move_constructor) {
                    valueptr = move_constructor(src);
                } else if (copy_constructor) {
                    valueptr = copy_constructor(src);
                } else {
#if defined(PYBIND11_DETAILED_ERROR_MESSAGES)
                    std::string type_name(tinfo->cpptype->name());
                    detail::clean_type_id(type_name);
                    throw cast_error("return_value_policy = move, but type " + type_name
                                     + " is neither movable nor copyable!");
#else
                    throw cast_error("return_value_policy = move, but type is neither "
                                     "movable nor copyable! "
                                     "(#define PYBIND11_DETAILED_ERROR_MESSAGES or compile in "
                                     "debug mode for details)");
#endif
                }
                wrapper->owned = true;
                break;

            case return_value_policy::reference_internal:
                valueptr = src;
                wrapper->owned = false;
                keep_alive_impl(inst, parent);
                break;

            default:
                throw cast_error("unhandled return_value_policy: should not happen!");
        }

        tinfo->init_instance(wrapper, existing_holder);

        return inst.release();
    }

    // Base methods for generic caster; there are overridden in copyable_holder_caster
    void load_value(value_and_holder &&v_h) {
        if (typeinfo->holder_enum_v == detail::holder_enum_t::smart_holder) {
            smart_holder_type_caster_support::value_and_holder_helper v_h_helper;
            v_h_helper.loaded_v_h = v_h;
            if (v_h_helper.have_holder()) {
                v_h_helper.throw_if_uninitialized_or_disowned_holder(cpptype->name());
                value = v_h_helper.holder().template as_raw_ptr_unowned<void>();
                return;
            }
        }
        auto *&vptr = v_h.value_ptr();
        // Lazy allocation for unallocated values:
        if (vptr == nullptr) {
            const auto *type = v_h.type ? v_h.type : typeinfo;
            if (type->operator_new) {
                vptr = type->operator_new(type->type_size);
            } else {
#if defined(__cpp_aligned_new) && (!defined(_MSC_VER) || _MSC_VER >= 1912)
                if (type->type_align > __STDCPP_DEFAULT_NEW_ALIGNMENT__) {
                    vptr = ::operator new(type->type_size, std::align_val_t(type->type_align));
                } else {
                    vptr = ::operator new(type->type_size);
                }
#else
                vptr = ::operator new(type->type_size);
#endif
            }
        }
        value = vptr;
    }
    bool try_implicit_casts(handle src, bool convert) {
        for (const auto &cast : typeinfo->implicit_casts) {
            type_caster_generic sub_caster(*cast.first);
            if (sub_caster.load(src, convert)) {
                value = cast.second(sub_caster.value);
                return true;
            }
        }
        return false;
    }
    bool try_direct_conversions(handle src) {
        for (auto &converter : *typeinfo->direct_conversions) {
            if (converter(src.ptr(), value)) {
                return true;
            }
        }
        return false;
    }
    bool try_cpp_conduit(handle src) {
        value = try_raw_pointer_ephemeral_from_cpp_conduit(src, cpptype);
        if (value != nullptr) {
            return true;
        }
        return false;
    }
    void check_holder_compat() {}
    bool set_foreign_holder(handle) { return true; }

    PYBIND11_NOINLINE static void *local_load(PyObject *src, const type_info *ti) {
        auto caster = type_caster_generic(ti);
        if (caster.load(src, false)) {
            return caster.value;
        }
        return nullptr;
    }

    /// Try to load with foreign typeinfo, if available. Used when there is no
    /// native typeinfo, or when the native one wasn't able to produce a value.
    PYBIND11_NOINLINE bool try_load_foreign_module_local(handle src) {
        constexpr auto *local_key = PYBIND11_MODULE_LOCAL_ID;
        const auto pytype = type::handle_of(src);
        if (!hasattr(pytype, local_key)) {
            return false;
        }

        type_info *foreign_typeinfo = reinterpret_borrow<capsule>(getattr(pytype, local_key));
        // Only consider this foreign loader if actually foreign and is a loader of the correct cpp
        // type
        if (foreign_typeinfo->module_local_load == &local_load
            || (cpptype && !same_type(*cpptype, *foreign_typeinfo->cpptype))) {
            return false;
        }

        if (auto *result = foreign_typeinfo->module_local_load(src.ptr(), foreign_typeinfo)) {
            value = result;
            return true;
        }
        return false;
    }

    // Implementation of `load`; this takes the type of `this` so that it can dispatch the relevant
    // bits of code between here and copyable_holder_caster where the two classes need different
    // logic (without having to resort to virtual inheritance).
    template <typename ThisT>
    PYBIND11_NOINLINE bool load_impl(handle src, bool convert) {
        auto &this_ = static_cast<ThisT &>(*this);
        if (!src) {
            return false;
        }
        if (!typeinfo) {
            return try_load_foreign_module_local(src) && this_.set_foreign_holder(src);
        }

        this_.check_holder_compat();

        PyTypeObject *srctype = Py_TYPE(src.ptr());

        // Case 1: If src is an exact type match for the target type then we can reinterpret_cast
        // the instance's value pointer to the target type:
        if (srctype == typeinfo->type) {
            this_.load_value(reinterpret_cast<instance *>(src.ptr())->get_value_and_holder());
            return true;
        }
        // Case 2: We have a derived class
        if (PyType_IsSubtype(srctype, typeinfo->type)) {
            const auto &bases = all_type_info(srctype);
            bool no_cpp_mi = typeinfo->simple_type;

            // Case 2a: the python type is a Python-inherited derived class that inherits from just
            // one simple (no MI) pybind11 class, or is an exact match, so the C++ instance is of
            // the right type and we can use reinterpret_cast.
            // (This is essentially the same as case 2b, but because not using multiple inheritance
            // is extremely common, we handle it specially to avoid the loop iterator and type
            // pointer lookup overhead)
            if (bases.size() == 1 && (no_cpp_mi || bases.front()->type == typeinfo->type)) {
                this_.load_value(reinterpret_cast<instance *>(src.ptr())->get_value_and_holder());
                return true;
            }
            // Case 2b: the python type inherits from multiple C++ bases.  Check the bases to see
            // if we can find an exact match (or, for a simple C++ type, an inherited match); if
            // so, we can safely reinterpret_cast to the relevant pointer.
            if (bases.size() > 1) {
                for (auto *base : bases) {
                    if (no_cpp_mi ? PyType_IsSubtype(base->type, typeinfo->type)
                                  : base->type == typeinfo->type) {
                        this_.load_value(
                            reinterpret_cast<instance *>(src.ptr())->get_value_and_holder(base));
                        return true;
                    }
                }
            }

            // Case 2c: C++ multiple inheritance is involved and we couldn't find an exact type
            // match in the registered bases, above, so try implicit casting (needed for proper C++
            // casting when MI is involved).
            if (this_.try_implicit_casts(src, convert)) {
                return true;
            }
        }

        // Perform an implicit conversion
        if (convert) {
            for (const auto &converter : typeinfo->implicit_conversions) {
                auto temp = reinterpret_steal<object>(converter(src.ptr(), typeinfo->type));
                if (load_impl<ThisT>(temp, false)) {
                    loader_life_support::add_patient(temp);
                    return true;
                }
            }
            if (this_.try_direct_conversions(src)) {
                return true;
            }
        }

        // Failed to match local typeinfo. Try again with global.
        if (typeinfo->module_local) {
            if (auto *gtype = get_global_type_info(*typeinfo->cpptype)) {
                typeinfo = gtype;
                return load_impl<ThisT>(src, false);
            }
        }

        // Global typeinfo has precedence over foreign module_local
        if (try_load_foreign_module_local(src)) {
            return this_.set_foreign_holder(src);
        }

        // Custom converters didn't take None, now we convert None to nullptr.
        if (src.is_none()) {
            // Defer accepting None to other overloads (if we aren't in convert mode):
            if (!convert) {
                return false;
            }
            value = nullptr;
            return true;
        }

        if (convert && cpptype && this_.try_cpp_conduit(src)) {
            return this_.set_foreign_holder(src);
        }

        return false;
    }

    const type_info *typeinfo = nullptr;
    const std::type_info *cpptype = nullptr;
    void *value = nullptr;
};

inline object cpp_conduit_method(handle self,
                                 const bytes &pybind11_platform_abi_id,
                                 const capsule &cpp_type_info_capsule,
                                 const bytes &pointer_kind) {
#ifdef PYBIND11_HAS_STRING_VIEW
    using cpp_str = std::string_view;
#else
    using cpp_str = std::string;
#endif
    if (cpp_str(pybind11_platform_abi_id) != PYBIND11_PLATFORM_ABI_ID) {
        return none();
    }
    if (std::strcmp(cpp_type_info_capsule.name(), typeid(std::type_info).name()) != 0) {
        return none();
    }
    if (cpp_str(pointer_kind) != "raw_pointer_ephemeral") {
        throw std::runtime_error("Invalid pointer_kind: \"" + std::string(pointer_kind) + "\"");
    }
    const auto *cpp_type_info = cpp_type_info_capsule.get_pointer<const std::type_info>();
    type_caster_generic caster(*cpp_type_info);
    if (!caster.load(self, false)) {
        return none();
    }
    return capsule(caster.value, cpp_type_info->name());
}

/**
 * Determine suitable casting operator for pointer-or-lvalue-casting type casters.  The type caster
 * needs to provide `operator T*()` and `operator T&()` operators.
 *
 * If the type supports moving the value away via an `operator T&&() &&` method, it should use
 * `movable_cast_op_type` instead.
 */
template <typename T>
using cast_op_type = conditional_t<std::is_pointer<remove_reference_t<T>>::value,
                                   typename std::add_pointer<intrinsic_t<T>>::type,
                                   typename std::add_lvalue_reference<intrinsic_t<T>>::type>;

/**
 * Determine suitable casting operator for a type caster with a movable value.  Such a type caster
 * needs to provide `operator T*()`, `operator T&()`, and `operator T&&() &&`.  The latter will be
 * called in appropriate contexts where the value can be moved rather than copied.
 *
 * These operator are automatically provided when using the PYBIND11_TYPE_CASTER macro.
 */
template <typename T>
using movable_cast_op_type
    = conditional_t<std::is_pointer<typename std::remove_reference<T>::type>::value,
                    typename std::add_pointer<intrinsic_t<T>>::type,
                    conditional_t<std::is_rvalue_reference<T>::value,
                                  typename std::add_rvalue_reference<intrinsic_t<T>>::type,
                                  typename std::add_lvalue_reference<intrinsic_t<T>>::type>>;

// Does the container have a mapped type and is it recursive?
// Implemented by specializations below.
template <typename Container, typename SFINAE = void>
struct container_mapped_type_traits {
    static constexpr bool has_mapped_type = false;
    static constexpr bool has_recursive_mapped_type = false;
};

template <typename Container>
struct container_mapped_type_traits<
    Container,
    typename std::enable_if<
        std::is_same<typename Container::mapped_type, Container>::value>::type> {
    static constexpr bool has_mapped_type = true;
    static constexpr bool has_recursive_mapped_type = true;
};

template <typename Container>
struct container_mapped_type_traits<
    Container,
    typename std::enable_if<
        negation<std::is_same<typename Container::mapped_type, Container>>::value>::type> {
    static constexpr bool has_mapped_type = true;
    static constexpr bool has_recursive_mapped_type = false;
};

// Does the container have a value type and is it recursive?
// Implemented by specializations below.
template <typename Container, typename SFINAE = void>
struct container_value_type_traits : std::false_type {
    static constexpr bool has_value_type = false;
    static constexpr bool has_recursive_value_type = false;
};

template <typename Container>
struct container_value_type_traits<
    Container,
    typename std::enable_if<
        std::is_same<typename Container::value_type, Container>::value>::type> {
    static constexpr bool has_value_type = true;
    static constexpr bool has_recursive_value_type = true;
};

template <typename Container>
struct container_value_type_traits<
    Container,
    typename std::enable_if<
        negation<std::is_same<typename Container::value_type, Container>>::value>::type> {
    static constexpr bool has_value_type = true;
    static constexpr bool has_recursive_value_type = false;
};

/*
 * Tag to be used for representing the bottom of recursively defined types.
 * Define this tag so we don't have to use void.
 */
struct recursive_bottom {};

/*
 * Implementation detail of `recursive_container_traits` below.
 * `T` is the `value_type` of the container, which might need to be modified to
 * avoid recursive types and const types.
 */
template <typename T, bool is_this_a_map>
struct impl_type_to_check_recursively {
    /*
     * If the container is recursive, then no further recursion should be done.
     */
    using if_recursive = recursive_bottom;
    /*
     * Otherwise yield `T` unchanged.
     */
    using if_not_recursive = T;
};

/*
 * For pairs - only as value type of a map -, the first type should remove the `const`.
 * Also, if the map is recursive, then the recursive checking should consider
 * the first type only.
 */
template <typename A, typename B>
struct impl_type_to_check_recursively<std::pair<A, B>, /* is_this_a_map = */ true> {
    using if_recursive = typename std::remove_const<A>::type;
    using if_not_recursive = std::pair<typename std::remove_const<A>::type, B>;
};

/*
 * Implementation of `recursive_container_traits` below.
 */
template <typename Container, typename SFINAE = void>
struct impl_recursive_container_traits {
    using type_to_check_recursively = recursive_bottom;
};

template <typename Container>
struct impl_recursive_container_traits<
    Container,
    typename std::enable_if<container_value_type_traits<Container>::has_value_type>::type> {
    static constexpr bool is_recursive
        = container_mapped_type_traits<Container>::has_recursive_mapped_type
          || container_value_type_traits<Container>::has_recursive_value_type;
    /*
     * This member dictates which type Pybind11 should check recursively in traits
     * such as `is_move_constructible`, `is_copy_constructible`, `is_move_assignable`, ...
     * Direct access to `value_type` should be avoided:
     * 1. `value_type` might recursively contain the type again
     * 2. `value_type` of STL map types is `std::pair<A const, B>`, the `const`
     *    should be removed.
     *
     */
    using type_to_check_recursively = typename std::conditional<
        is_recursive,
        typename impl_type_to_check_recursively<
            typename Container::value_type,
            container_mapped_type_traits<Container>::has_mapped_type>::if_recursive,
        typename impl_type_to_check_recursively<
            typename Container::value_type,
            container_mapped_type_traits<Container>::has_mapped_type>::if_not_recursive>::type;
};

/*
 * This trait defines the `type_to_check_recursively` which is needed to properly
 * handle recursively defined traits such as `is_move_constructible` without going
 * into an infinite recursion.
 * Should be used instead of directly accessing the `value_type`.
 * It cancels the recursion by returning the `recursive_bottom` tag.
 *
 * The default definition of `type_to_check_recursively` is as follows:
 *
 * 1. By default, it is `recursive_bottom`, so that the recursion is canceled.
 * 2. If the type is non-recursive and defines a `value_type`, then the `value_type` is used.
 *    If the `value_type` is a pair and a `mapped_type` is defined,
 *    then the `const` is removed from the first type.
 * 3. If the type is recursive and `value_type` is not a pair, then `recursive_bottom` is returned.
 * 4. If the type is recursive and `value_type` is a pair and a `mapped_type` is defined,
 *    then `const` is removed from the first type and the first type is returned.
 *
 * This behavior can be extended by the user as seen in test_stl_binders.cpp.
 *
 * This struct is exactly the same as impl_recursive_container_traits.
 * The duplication achieves that user-defined specializations don't compete
 * with internal specializations, but take precedence.
 */
template <typename Container, typename SFINAE = void>
struct recursive_container_traits : impl_recursive_container_traits<Container> {};

template <typename T>
struct is_move_constructible
    : all_of<std::is_move_constructible<T>,
             is_move_constructible<
                 typename recursive_container_traits<T>::type_to_check_recursively>> {};

template <>
struct is_move_constructible<recursive_bottom> : std::true_type {};

// Likewise for std::pair
// (after C++17 it is mandatory that the move constructor not exist when the two types aren't
// themselves move constructible, but this can not be relied upon when T1 or T2 are themselves
// containers).
template <typename T1, typename T2>
struct is_move_constructible<std::pair<T1, T2>>
    : all_of<is_move_constructible<T1>, is_move_constructible<T2>> {};

// std::is_copy_constructible isn't quite enough: it lets std::vector<T> (and similar) through when
// T is non-copyable, but code containing such a copy constructor fails to actually compile.
template <typename T>
struct is_copy_constructible
    : all_of<std::is_copy_constructible<T>,
             is_copy_constructible<
                 typename recursive_container_traits<T>::type_to_check_recursively>> {};

template <>
struct is_copy_constructible<recursive_bottom> : std::true_type {};

// Likewise for std::pair
// (after C++17 it is mandatory that the copy constructor not exist when the two types aren't
// themselves copy constructible, but this can not be relied upon when T1 or T2 are themselves
// containers).
template <typename T1, typename T2>
struct is_copy_constructible<std::pair<T1, T2>>
    : all_of<is_copy_constructible<T1>, is_copy_constructible<T2>> {};

// The same problems arise with std::is_copy_assignable, so we use the same workaround.
template <typename T>
struct is_copy_assignable
    : all_of<
          std::is_copy_assignable<T>,
          is_copy_assignable<typename recursive_container_traits<T>::type_to_check_recursively>> {
};

template <>
struct is_copy_assignable<recursive_bottom> : std::true_type {};

template <typename T1, typename T2>
struct is_copy_assignable<std::pair<T1, T2>>
    : all_of<is_copy_assignable<T1>, is_copy_assignable<T2>> {};

PYBIND11_NAMESPACE_END(detail)

// polymorphic_type_hook<itype>::get(src, tinfo) determines whether the object pointed
// to by `src` actually is an instance of some class derived from `itype`.
// If so, it sets `tinfo` to point to the std::type_info representing that derived
// type, and returns a pointer to the start of the most-derived object of that type
// (in which `src` is a subobject; this will be the same address as `src` in most
// single inheritance cases). If not, or if `src` is nullptr, it simply returns `src`
// and leaves `tinfo` at its default value of nullptr.
//
// The default polymorphic_type_hook just returns src. A specialization for polymorphic
// types determines the runtime type of the passed object and adjusts the this-pointer
// appropriately via dynamic_cast<void*>. This is what enables a C++ Animal* to appear
// to Python as a Dog (if Dog inherits from Animal, Animal is polymorphic, Dog is
// registered with pybind11, and this Animal is in fact a Dog).
//
// You may specialize polymorphic_type_hook yourself for types that want to appear
// polymorphic to Python but do not use C++ RTTI. (This is a not uncommon pattern
// in performance-sensitive applications, used most notably in LLVM.)
//
// polymorphic_type_hook_base allows users to specialize polymorphic_type_hook with
// std::enable_if. User provided specializations will always have higher priority than
// the default implementation and specialization provided in polymorphic_type_hook_base.
template <typename itype, typename SFINAE = void>
struct polymorphic_type_hook_base {
    static const void *get(const itype *src, const std::type_info *&) { return src; }
};
template <typename itype>
struct polymorphic_type_hook_base<itype, detail::enable_if_t<std::is_polymorphic<itype>::value>> {
    static const void *get(const itype *src, const std::type_info *&type) {
        type = src ? &typeid(*src) : nullptr;
        return dynamic_cast<const void *>(src);
    }
};
template <typename itype, typename SFINAE = void>
struct polymorphic_type_hook : public polymorphic_type_hook_base<itype> {};

PYBIND11_NAMESPACE_BEGIN(detail)

template <typename itype>
cast_sources::cast_sources(const itype *ptr) : original{ptr, &typeid(itype)} {
    // If this is a base pointer to a derived type, and the derived type is
    // registered with pybind11, we want to make the full derived object
    // available. In the typical case where itype is polymorphic, we get the
    // correct derived pointer (which may be != base pointer) by a dynamic_cast
    // to most derived type. If itype is not polymorphic, a user-provided
    // specialization of polymorphic_type_hook can do the same thing.
    // If there is no downcast to perform, then the default hook will leave
    // derived.type set to nullptr, which causes us to ignore derived.obj.
    downcast.cppobj = polymorphic_type_hook<itype>::get(ptr, downcast.cpptype);
    result = resolve();
}

/// Generic type caster for objects stored on the heap
template <typename type>
class type_caster_base : public type_caster_generic {
    using itype = intrinsic_t<type>;

public:
    static constexpr auto name = const_name<type>();

    type_caster_base() : type_caster_base(typeid(type)) {}
    explicit type_caster_base(const std::type_info &info) : type_caster_generic(info) {}

    // Wrap the generic cast_sources to be only constructible from the type
    // that's correct in this context, so you can't use type_caster_base<A>
    // to convert an unrelated B* to Python.
    struct cast_sources : detail::cast_sources {
        explicit cast_sources(const itype *ptr) : detail::cast_sources(ptr) {}
    };

    static handle cast(const itype &src, return_value_policy policy, handle parent) {
        if (policy == return_value_policy::automatic
            || policy == return_value_policy::automatic_reference) {
            policy = return_value_policy::copy;
        }
        return cast(std::addressof(src), policy, parent);
    }

    static handle cast(itype &&src, return_value_policy, handle parent) {
        return cast(std::addressof(src), return_value_policy::move, parent);
    }

    static handle cast(const itype *src, return_value_policy policy, handle parent) {
        return cast(cast_sources{src}, policy, parent);
    }

    static handle cast(const cast_sources &srcs, return_value_policy policy, handle parent) {
        return type_caster_generic::cast(srcs,
                                         policy,
                                         parent,
                                         make_copy_constructor((const itype *) nullptr),
                                         make_move_constructor((const itype *) nullptr));
    }

    static handle cast_holder(const itype *src, const void *holder) {
        return cast_holder(cast_sources{src}, holder);
    }

    static handle cast_holder(const cast_sources &srcs, const void *holder) {
        auto policy = return_value_policy::take_ownership;
        return type_caster_generic::cast(srcs, policy, {}, nullptr, nullptr, holder);
    }

    template <typename T>
    using cast_op_type = detail::cast_op_type<T>;

    // NOLINTNEXTLINE(google-explicit-constructor)
    operator itype *() { return (type *) value; }
    // NOLINTNEXTLINE(google-explicit-constructor)
    operator itype &() {
        if (!value) {
            throw reference_cast_error();
        }
        return *((itype *) value);
    }

protected:
    using Constructor = void *(*) (const void *);

    /* Only enabled when the types are {copy,move}-constructible *and* when the type
       does not have a private operator new implementation. A comma operator is used in the
       decltype argument to apply SFINAE to the public copy/move constructors.*/
    template <typename T, typename = enable_if_t<is_copy_constructible<T>::value>>
    static auto make_copy_constructor(const T *)
        -> decltype(new T(std::declval<const T>()), Constructor{}) {
        return [](const void *arg) -> void * { return new T(*reinterpret_cast<const T *>(arg)); };
    }

    template <typename T, typename = enable_if_t<is_move_constructible<T>::value>>
    static auto make_move_constructor(const T *)
        -> decltype(new T(std::declval<T &&>()), Constructor{}) {
        return [](const void *arg) -> void * {
            return new T(std::move(*const_cast<T *>(reinterpret_cast<const T *>(arg))));
        };
    }

    static Constructor make_copy_constructor(...) { return nullptr; }
    static Constructor make_move_constructor(...) { return nullptr; }
};

inline std::string quote_cpp_type_name(const std::string &cpp_type_name) {
    return cpp_type_name; // No-op for now. See PR #4888
}

PYBIND11_NOINLINE std::string type_info_description(const std::type_info &ti) {
    if (auto *type_data = get_type_info(ti)) {
        handle th(reinterpret_cast<PyObject *>(type_data->type));
        return th.attr("__module__").cast<std::string>() + '.'
               + th.attr("__qualname__").cast<std::string>();
    }
    return quote_cpp_type_name(clean_type_id(ti.name()));
}

PYBIND11_NAMESPACE_END(detail)
PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)

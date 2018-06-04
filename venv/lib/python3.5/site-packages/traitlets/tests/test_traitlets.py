# encoding: utf-8
"""Tests for traitlets.traitlets."""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
# 
# Adapted from enthought.traits, Copyright (c) Enthought, Inc.,
# also under the terms of the Modified BSD License.

import pickle
import re
import sys
from ._warnings import expected_warnings

from unittest import TestCase
import pytest
from pytest import mark

from traitlets import (
    HasTraits, MetaHasTraits, TraitType, Any, Bool, CBytes, Dict, Enum,
    Int, CInt, Long, CLong, Integer, Float, CFloat, Complex, Bytes, Unicode,
    TraitError, Union, All, Undefined, Type, This, Instance, TCPAddress,
    List, Tuple, ObjectName, DottedObjectName, CRegExp, link, directional_link,
    ForwardDeclaredType, ForwardDeclaredInstance, validate, observe, default,
    observe_compat, BaseDescriptor, HasDescriptors,
)

import six

def change_dict(*ordered_values):
    change_names = ('name', 'old', 'new', 'owner', 'type')
    return dict(zip(change_names, ordered_values))

#-----------------------------------------------------------------------------
# Helper classes for testing
#-----------------------------------------------------------------------------


class HasTraitsStub(HasTraits):

    def notify_change(self, change):
        self._notify_name = change['name']
        self._notify_old = change['old']
        self._notify_new = change['new']
        self._notify_type = change['type']


#-----------------------------------------------------------------------------
# Test classes
#-----------------------------------------------------------------------------


class TestTraitType(TestCase):

    def test_get_undefined(self):
        class A(HasTraits):
            a = TraitType
        a = A()
        with self.assertRaises(TraitError):
            a.a

    def test_set(self):
        class A(HasTraitsStub):
            a = TraitType

        a = A()
        a.a = 10
        self.assertEqual(a.a, 10)
        self.assertEqual(a._notify_name, 'a')
        self.assertEqual(a._notify_old, Undefined)
        self.assertEqual(a._notify_new, 10)

    def test_validate(self):
        class MyTT(TraitType):
            def validate(self, inst, value):
                return -1
        class A(HasTraitsStub):
            tt = MyTT

        a = A()
        a.tt = 10
        self.assertEqual(a.tt, -1)

    def test_default_validate(self):
        class MyIntTT(TraitType):
            def validate(self, obj, value):
                if isinstance(value, int):
                    return value
                self.error(obj, value)
        class A(HasTraits):
            tt = MyIntTT(10)
        a = A()
        self.assertEqual(a.tt, 10)

        # Defaults are validated when the HasTraits is instantiated
        class B(HasTraits):
            tt = MyIntTT('bad default')
        self.assertRaises(TraitError, B)

    def test_info(self):
        class A(HasTraits):
            tt = TraitType
        a = A()
        self.assertEqual(A.tt.info(), 'any value')

    def test_error(self):
        class A(HasTraits):
            tt = TraitType
        a = A()
        self.assertRaises(TraitError, A.tt.error, a, 10)

    def test_deprecated_dynamic_initializer(self):
        class A(HasTraits):
            x = Int(10)
            def _x_default(self):
                return 11
        class B(A):
            x = Int(20)
        class C(A):
            def _x_default(self):
                return 21

        a = A()
        self.assertEqual(a._trait_values, {})
        self.assertEqual(a.x, 11)
        self.assertEqual(a._trait_values, {'x': 11})
        b = B()
        self.assertEqual(b.x, 20)
        self.assertEqual(b._trait_values, {'x': 20})
        c = C()
        self.assertEqual(c._trait_values, {})
        self.assertEqual(c.x, 21)
        self.assertEqual(c._trait_values, {'x': 21})
        # Ensure that the base class remains unmolested when the _default
        # initializer gets overridden in a subclass.
        a = A()
        c = C()
        self.assertEqual(a._trait_values, {})
        self.assertEqual(a.x, 11)
        self.assertEqual(a._trait_values, {'x': 11})

    def test_dynamic_initializer(self):

        class A(HasTraits):
            x = Int(10)

            @default('x')
            def _default_x(self):
                return 11

        class B(A):
            x = Int(20)

        class C(A):

            @default('x')
            def _default_x(self):
                return 21

        a = A()
        self.assertEqual(a._trait_values, {})
        self.assertEqual(a.x, 11)
        self.assertEqual(a._trait_values, {'x': 11})
        b = B()
        self.assertEqual(b.x, 20)
        self.assertEqual(b._trait_values, {'x': 20})
        c = C()
        self.assertEqual(c._trait_values, {})
        self.assertEqual(c.x, 21)
        self.assertEqual(c._trait_values, {'x': 21})
        # Ensure that the base class remains unmolested when the _default
        # initializer gets overridden in a subclass.
        a = A()
        c = C()
        self.assertEqual(a._trait_values, {})
        self.assertEqual(a.x, 11)
        self.assertEqual(a._trait_values, {'x': 11})

    def test_tag_metadata(self):
        class MyIntTT(TraitType):
            metadata = {'a': 1, 'b': 2}
        a = MyIntTT(10).tag(b=3, c=4)
        self.assertEqual(a.metadata, {'a': 1, 'b': 3, 'c': 4})

    def test_metadata_localized_instance(self):
        class MyIntTT(TraitType):
            metadata = {'a': 1, 'b': 2}
        a = MyIntTT(10)
        b = MyIntTT(10)
        a.metadata['c'] = 3
        # make sure that changing a's metadata didn't change b's metadata
        self.assertNotIn('c', b.metadata)

    def test_union_metadata(self):
        class Foo(HasTraits):
            bar = (Int().tag(ta=1) | Dict().tag(ta=2, ti='b')).tag(ti='a')
        foo = Foo()
        # At this point, no value has been set for bar, so value-specific
        # is not set.
        self.assertEqual(foo.trait_metadata('bar', 'ta'), None)
        self.assertEqual(foo.trait_metadata('bar', 'ti'), 'a')
        foo.bar = {}
        self.assertEqual(foo.trait_metadata('bar', 'ta'), 2)
        self.assertEqual(foo.trait_metadata('bar', 'ti'), 'b')
        foo.bar = 1
        self.assertEqual(foo.trait_metadata('bar', 'ta'), 1)
        self.assertEqual(foo.trait_metadata('bar', 'ti'), 'a')

    def test_union_default_value(self):
        class Foo(HasTraits):
            bar = Union([Dict(), Int()], default_value=1)
        foo = Foo()
        self.assertEqual(foo.bar, 1)

    def test_deprecated_metadata_access(self):
        class MyIntTT(TraitType):
            metadata = {'a': 1, 'b': 2}
        a = MyIntTT(10)
        with expected_warnings(["use the instance .metadata dictionary directly"]*2):
            a.set_metadata('key', 'value')
            v = a.get_metadata('key')
        self.assertEqual(v, 'value')
        with expected_warnings(["use the instance .help string directly"]*2):
            a.set_metadata('help', 'some help')
            v = a.get_metadata('help')
        self.assertEqual(v, 'some help')

    def test_trait_types_deprecated(self):
        with expected_warnings(["Traits should be given as instances"]):
            class C(HasTraits):
                t = Int

    def test_trait_types_list_deprecated(self):
        with expected_warnings(["Traits should be given as instances"]):
            class C(HasTraits):
                t = List(Int)

    def test_trait_types_tuple_deprecated(self):
        with expected_warnings(["Traits should be given as instances"]):
            class C(HasTraits):
                t = Tuple(Int)

    def test_trait_types_dict_deprecated(self):
        with expected_warnings(["Traits should be given as instances"]):
            class C(HasTraits):
                t = Dict(Int)

class TestHasDescriptorsMeta(TestCase):

    def test_metaclass(self):
        self.assertEqual(type(HasTraits), MetaHasTraits)

        class A(HasTraits):
            a = Int()

        a = A()
        self.assertEqual(type(a.__class__), MetaHasTraits)
        self.assertEqual(a.a,0)
        a.a = 10
        self.assertEqual(a.a,10)

        class B(HasTraits):
            b = Int()

        b = B()
        self.assertEqual(b.b,0)
        b.b = 10
        self.assertEqual(b.b,10)

        class C(HasTraits):
            c = Int(30)

        c = C()
        self.assertEqual(c.c,30)
        c.c = 10
        self.assertEqual(c.c,10)

    def test_this_class(self):
        class A(HasTraits):
            t = This()
            tt = This()
        class B(A):
            tt = This()
            ttt = This()
        self.assertEqual(A.t.this_class, A)
        self.assertEqual(B.t.this_class, A)
        self.assertEqual(B.tt.this_class, B)
        self.assertEqual(B.ttt.this_class, B)

class TestHasDescriptors(TestCase):

    def test_setup_instance(self):

        class FooDescriptor(BaseDescriptor):

            def instance_init(self, inst):
                foo = inst.foo # instance should have the attr

        class HasFooDescriptors(HasDescriptors):

            fd = FooDescriptor()

            def setup_instance(self, *args, **kwargs):
                self.foo = kwargs.get('foo', None)
                super(HasFooDescriptors, self).setup_instance(*args, **kwargs)

        hfd = HasFooDescriptors(foo='bar')

class TestHasTraitsNotify(TestCase):

    def setUp(self):
        self._notify1 = []
        self._notify2 = []

    def notify1(self, name, old, new):
        self._notify1.append((name, old, new))

    def notify2(self, name, old, new):
        self._notify2.append((name, old, new))

    def test_notify_all(self):

        class A(HasTraits):
            a = Int()
            b = Float()

        a = A()
        a.on_trait_change(self.notify1)
        a.a = 0
        self.assertEqual(len(self._notify1),0)
        a.b = 0.0
        self.assertEqual(len(self._notify1),0)
        a.a = 10
        self.assertTrue(('a',0,10) in self._notify1)
        a.b = 10.0
        self.assertTrue(('b',0.0,10.0) in self._notify1)
        self.assertRaises(TraitError,setattr,a,'a','bad string')
        self.assertRaises(TraitError,setattr,a,'b','bad string')
        self._notify1 = []
        a.on_trait_change(self.notify1,remove=True)
        a.a = 20
        a.b = 20.0
        self.assertEqual(len(self._notify1),0)

    def test_notify_one(self):

        class A(HasTraits):
            a = Int()
            b = Float()

        a = A()
        a.on_trait_change(self.notify1, 'a')
        a.a = 0
        self.assertEqual(len(self._notify1),0)
        a.a = 10
        self.assertTrue(('a',0,10) in self._notify1)
        self.assertRaises(TraitError,setattr,a,'a','bad string')

    def test_subclass(self):

        class A(HasTraits):
            a = Int()

        class B(A):
            b = Float()

        b = B()
        self.assertEqual(b.a,0)
        self.assertEqual(b.b,0.0)
        b.a = 100
        b.b = 100.0
        self.assertEqual(b.a,100)
        self.assertEqual(b.b,100.0)

    def test_notify_subclass(self):

        class A(HasTraits):
            a = Int()

        class B(A):
            b = Float()

        b = B()
        b.on_trait_change(self.notify1, 'a')
        b.on_trait_change(self.notify2, 'b')
        b.a = 0
        b.b = 0.0
        self.assertEqual(len(self._notify1),0)
        self.assertEqual(len(self._notify2),0)
        b.a = 10
        b.b = 10.0
        self.assertTrue(('a',0,10) in self._notify1)
        self.assertTrue(('b',0.0,10.0) in self._notify2)

    def test_static_notify(self):

        class A(HasTraits):
            a = Int()
            _notify1 = []
            def _a_changed(self, name, old, new):
                self._notify1.append((name, old, new))

        a = A()
        a.a = 0
        # This is broken!!!
        self.assertEqual(len(a._notify1),0)
        a.a = 10
        self.assertTrue(('a',0,10) in a._notify1)

        class B(A):
            b = Float()
            _notify2 = []
            def _b_changed(self, name, old, new):
                self._notify2.append((name, old, new))

        b = B()
        b.a = 10
        b.b = 10.0
        self.assertTrue(('a',0,10) in b._notify1)
        self.assertTrue(('b',0.0,10.0) in b._notify2)

    def test_notify_args(self):

        def callback0():
            self.cb = ()
        def callback1(name):
            self.cb = (name,)
        def callback2(name, new):
            self.cb = (name, new)
        def callback3(name, old, new):
            self.cb = (name, old, new)
        def callback4(name, old, new, obj):
            self.cb = (name, old, new, obj)

        class A(HasTraits):
            a = Int()

        a = A()
        a.on_trait_change(callback0, 'a')
        a.a = 10
        self.assertEqual(self.cb,())
        a.on_trait_change(callback0, 'a', remove=True)

        a.on_trait_change(callback1, 'a')
        a.a = 100
        self.assertEqual(self.cb,('a',))
        a.on_trait_change(callback1, 'a', remove=True)

        a.on_trait_change(callback2, 'a')
        a.a = 1000
        self.assertEqual(self.cb,('a',1000))
        a.on_trait_change(callback2, 'a', remove=True)

        a.on_trait_change(callback3, 'a')
        a.a = 10000
        self.assertEqual(self.cb,('a',1000,10000))
        a.on_trait_change(callback3, 'a', remove=True)

        a.on_trait_change(callback4, 'a')
        a.a = 100000
        self.assertEqual(self.cb,('a',10000,100000,a))
        self.assertEqual(len(a._trait_notifiers['a']['change']), 1)
        a.on_trait_change(callback4, 'a', remove=True)

        self.assertEqual(len(a._trait_notifiers['a']['change']), 0)

    def test_notify_only_once(self):

        class A(HasTraits):
            listen_to = ['a']
            
            a = Int(0)
            b = 0
            
            def __init__(self, **kwargs):
                super(A, self).__init__(**kwargs)
                self.on_trait_change(self.listener1, ['a'])
            
            def listener1(self, name, old, new):
                self.b += 1

        class B(A):
                    
            c = 0
            d = 0
            
            def __init__(self, **kwargs):
                super(B, self).__init__(**kwargs)
                self.on_trait_change(self.listener2)
            
            def listener2(self, name, old, new):
                self.c += 1
            
            def _a_changed(self, name, old, new):
                self.d += 1

        b = B()
        b.a += 1
        self.assertEqual(b.b, b.c)
        self.assertEqual(b.b, b.d)
        b.a += 1
        self.assertEqual(b.b, b.c)
        self.assertEqual(b.b, b.d)

class TestObserveDecorator(TestCase):

    def setUp(self):
        self._notify1 = []
        self._notify2 = []

    def notify1(self, change):
        self._notify1.append(change)

    def notify2(self, change):
        self._notify2.append(change)

    def test_notify_all(self):

        class A(HasTraits):
            a = Int()
            b = Float()

        a = A()
        a.observe(self.notify1)
        a.a = 0
        self.assertEqual(len(self._notify1),0)
        a.b = 0.0
        self.assertEqual(len(self._notify1),0)
        a.a = 10
        change = change_dict('a', 0, 10, a, 'change')
        self.assertTrue(change in self._notify1)
        a.b = 10.0
        change = change_dict('b', 0.0, 10.0, a, 'change')
        self.assertTrue(change in self._notify1)
        self.assertRaises(TraitError,setattr,a,'a','bad string')
        self.assertRaises(TraitError,setattr,a,'b','bad string')
        self._notify1 = []
        a.unobserve(self.notify1)
        a.a = 20
        a.b = 20.0
        self.assertEqual(len(self._notify1),0)

    def test_notify_one(self):

        class A(HasTraits):
            a = Int()
            b = Float()

        a = A()
        a.observe(self.notify1, 'a')
        a.a = 0
        self.assertEqual(len(self._notify1),0)
        a.a = 10
        change = change_dict('a', 0, 10, a, 'change')
        self.assertTrue(change in self._notify1)
        self.assertRaises(TraitError,setattr,a,'a','bad string')

    def test_subclass(self):

        class A(HasTraits):
            a = Int()

        class B(A):
            b = Float()

        b = B()
        self.assertEqual(b.a,0)
        self.assertEqual(b.b,0.0)
        b.a = 100
        b.b = 100.0
        self.assertEqual(b.a,100)
        self.assertEqual(b.b,100.0)

    def test_notify_subclass(self):

        class A(HasTraits):
            a = Int()

        class B(A):
            b = Float()

        b = B()
        b.observe(self.notify1, 'a')
        b.observe(self.notify2, 'b')
        b.a = 0
        b.b = 0.0
        self.assertEqual(len(self._notify1),0)
        self.assertEqual(len(self._notify2),0)
        b.a = 10
        b.b = 10.0
        change = change_dict('a', 0, 10, b, 'change')
        self.assertTrue(change in self._notify1)
        change = change_dict('b', 0.0, 10.0, b, 'change')
        self.assertTrue(change in self._notify2)

    def test_static_notify(self):

        class A(HasTraits):
            a = Int()
            b = Int()
            _notify1 = []
            _notify_any = []
            
            @observe('a')
            def _a_changed(self, change):
                self._notify1.append(change)

            @observe(All)
            def _any_changed(self, change):
                self._notify_any.append(change)

        a = A()
        a.a = 0
        self.assertEqual(len(a._notify1),0)
        a.a = 10
        change = change_dict('a', 0, 10, a, 'change')
        self.assertTrue(change in a._notify1)
        a.b = 1
        self.assertEqual(len(a._notify_any), 2)
        change = change_dict('b', 0, 1, a, 'change')
        self.assertTrue(change in a._notify_any)

        class B(A):
            b = Float()
            _notify2 = []
            @observe('b')
            def _b_changed(self, change):
                self._notify2.append(change)

        b = B()
        b.a = 10
        b.b = 10.0
        change = change_dict('a', 0, 10, b, 'change')
        self.assertTrue(change in b._notify1)
        change = change_dict('b', 0.0, 10.0, b, 'change')
        self.assertTrue(change in b._notify2)

    def test_notify_args(self):

        def callback0():
            self.cb = ()
        def callback1(change):
            self.cb = change

        class A(HasTraits):
            a = Int()

        a = A()
        a.on_trait_change(callback0, 'a')
        a.a = 10
        self.assertEqual(self.cb,())
        a.unobserve(callback0, 'a')

        a.observe(callback1, 'a')
        a.a = 100
        change = change_dict('a', 10, 100, a, 'change')
        self.assertEqual(self.cb, change)
        self.assertEqual(len(a._trait_notifiers['a']['change']), 1)
        a.unobserve(callback1, 'a')

        self.assertEqual(len(a._trait_notifiers['a']['change']), 0)

    def test_notify_only_once(self):

        class A(HasTraits):
            listen_to = ['a']
            
            a = Int(0)
            b = 0
            
            def __init__(self, **kwargs):
                super(A, self).__init__(**kwargs)
                self.observe(self.listener1, ['a'])
            
            def listener1(self, change):
                self.b += 1

        class B(A):
                    
            c = 0
            d = 0
            
            def __init__(self, **kwargs):
                super(B, self).__init__(**kwargs)
                self.observe(self.listener2)
            
            def listener2(self, change):
                self.c += 1
            
            @observe('a')
            def _a_changed(self, change):
                self.d += 1

        b = B()
        b.a += 1
        self.assertEqual(b.b, b.c)
        self.assertEqual(b.b, b.d)
        b.a += 1
        self.assertEqual(b.b, b.c)
        self.assertEqual(b.b, b.d)


class TestHasTraits(TestCase):

    def test_trait_names(self):
        class A(HasTraits):
            i = Int()
            f = Float()
        a = A()
        self.assertEqual(sorted(a.trait_names()),['f','i'])
        self.assertEqual(sorted(A.class_trait_names()),['f','i'])
        self.assertTrue(a.has_trait('f'))
        self.assertFalse(a.has_trait('g'))

    def test_trait_metadata_deprecated(self):
        with expected_warnings(['metadata should be set using the \.tag\(\) method']):
            class A(HasTraits):
                i = Int(config_key='MY_VALUE')
        a = A()
        self.assertEqual(a.trait_metadata('i','config_key'), 'MY_VALUE')

    def test_trait_metadata(self):
        class A(HasTraits):
            i = Int().tag(config_key='MY_VALUE')
        a = A()
        self.assertEqual(a.trait_metadata('i','config_key'), 'MY_VALUE')

    def test_trait_metadata_default(self):
        class A(HasTraits):
            i = Int()
        a = A()
        self.assertEqual(a.trait_metadata('i', 'config_key'), None)
        self.assertEqual(a.trait_metadata('i', 'config_key', 'default'), 'default')

    def test_traits(self):
        class A(HasTraits):
            i = Int()
            f = Float()
        a = A()
        self.assertEqual(a.traits(), dict(i=A.i, f=A.f))
        self.assertEqual(A.class_traits(), dict(i=A.i, f=A.f))

    def test_traits_metadata(self):
        class A(HasTraits):
            i = Int().tag(config_key='VALUE1', other_thing='VALUE2')
            f = Float().tag(config_key='VALUE3', other_thing='VALUE2')
            j = Int(0)
        a = A()
        self.assertEqual(a.traits(), dict(i=A.i, f=A.f, j=A.j))
        traits = a.traits(config_key='VALUE1', other_thing='VALUE2')
        self.assertEqual(traits, dict(i=A.i))

        # This passes, but it shouldn't because I am replicating a bug in
        # traits.
        traits = a.traits(config_key=lambda v: True)
        self.assertEqual(traits, dict(i=A.i, f=A.f, j=A.j))

    def test_traits_metadata_deprecated(self):
        with expected_warnings(['metadata should be set using the \.tag\(\) method']*2):
            class A(HasTraits):
                i = Int(config_key='VALUE1', other_thing='VALUE2')
                f = Float(config_key='VALUE3', other_thing='VALUE2')
                j = Int(0)
        a = A()
        self.assertEqual(a.traits(), dict(i=A.i, f=A.f, j=A.j))
        traits = a.traits(config_key='VALUE1', other_thing='VALUE2')
        self.assertEqual(traits, dict(i=A.i))

        # This passes, but it shouldn't because I am replicating a bug in
        # traits.
        traits = a.traits(config_key=lambda v: True)
        self.assertEqual(traits, dict(i=A.i, f=A.f, j=A.j))


    def test_init(self):
        class A(HasTraits):
            i = Int()
            x = Float()
        a = A(i=1, x=10.0)
        self.assertEqual(a.i, 1)
        self.assertEqual(a.x, 10.0)

    def test_positional_args(self):
        class A(HasTraits):
            i = Int(0)
            def __init__(self, i):
                super(A, self).__init__()
                self.i = i
        
        a = A(5)
        self.assertEqual(a.i, 5)
        # should raise TypeError if no positional arg given
        self.assertRaises(TypeError, A)

#-----------------------------------------------------------------------------
# Tests for specific trait types
#-----------------------------------------------------------------------------


class TestType(TestCase):

    def test_default(self):

        class B(object): pass
        class A(HasTraits):
            klass = Type(allow_none=True)

        a = A()
        self.assertEqual(a.klass, object)

        a.klass = B
        self.assertEqual(a.klass, B)
        self.assertRaises(TraitError, setattr, a, 'klass', 10)

    def test_default_options(self):

        class B(object): pass
        class C(B): pass
        class A(HasTraits):
            # Different possible combinations of options for default_value
            # and klass. default_value=None is only valid with allow_none=True.
            k1 = Type()
            k2 = Type(None, allow_none=True)
            k3 = Type(B)
            k4 = Type(klass=B)
            k5 = Type(default_value=None, klass=B, allow_none=True)
            k6 = Type(default_value=C, klass=B)

        self.assertIs(A.k1.default_value, object)
        self.assertIs(A.k1.klass, object)
        self.assertIs(A.k2.default_value, None)
        self.assertIs(A.k2.klass, object)
        self.assertIs(A.k3.default_value, B)
        self.assertIs(A.k3.klass, B)
        self.assertIs(A.k4.default_value, B)
        self.assertIs(A.k4.klass, B)
        self.assertIs(A.k5.default_value, None)
        self.assertIs(A.k5.klass, B)
        self.assertIs(A.k6.default_value, C)
        self.assertIs(A.k6.klass, B)

        a = A()
        self.assertIs(a.k1, object)
        self.assertIs(a.k2, None)
        self.assertIs(a.k3, B)
        self.assertIs(a.k4, B)
        self.assertIs(a.k5, None)
        self.assertIs(a.k6, C)

    def test_value(self):

        class B(object): pass
        class C(object): pass
        class A(HasTraits):
            klass = Type(B)

        a = A()
        self.assertEqual(a.klass, B)
        self.assertRaises(TraitError, setattr, a, 'klass', C)
        self.assertRaises(TraitError, setattr, a, 'klass', object)
        a.klass = B

    def test_allow_none(self):

        class B(object): pass
        class C(B): pass
        class A(HasTraits):
            klass = Type(B)

        a = A()
        self.assertEqual(a.klass, B)
        self.assertRaises(TraitError, setattr, a, 'klass', None)
        a.klass = C
        self.assertEqual(a.klass, C)

    def test_validate_klass(self):

        class A(HasTraits):
            klass = Type('no strings allowed')

        self.assertRaises(ImportError, A)

        class A(HasTraits):
            klass = Type('rub.adub.Duck')

        self.assertRaises(ImportError, A)

    def test_validate_default(self):

        class B(object): pass
        class A(HasTraits):
            klass = Type('bad default', B)

        self.assertRaises(ImportError, A)

        class C(HasTraits):
            klass = Type(None, B)

        self.assertRaises(TraitError, C)

    def test_str_klass(self):

        class A(HasTraits):
            klass = Type('ipython_genutils.ipstruct.Struct')

        from ipython_genutils.ipstruct import Struct
        a = A()
        a.klass = Struct
        self.assertEqual(a.klass, Struct)

        self.assertRaises(TraitError, setattr, a, 'klass', 10)

    def test_set_str_klass(self):

        class A(HasTraits):
            klass = Type()

        a = A(klass='ipython_genutils.ipstruct.Struct')
        from ipython_genutils.ipstruct import Struct
        self.assertEqual(a.klass, Struct)

class TestInstance(TestCase):

    def test_basic(self):
        class Foo(object): pass
        class Bar(Foo): pass
        class Bah(object): pass

        class A(HasTraits):
            inst = Instance(Foo, allow_none=True)

        a = A()
        self.assertTrue(a.inst is None)
        a.inst = Foo()
        self.assertTrue(isinstance(a.inst, Foo))
        a.inst = Bar()
        self.assertTrue(isinstance(a.inst, Foo))
        self.assertRaises(TraitError, setattr, a, 'inst', Foo)
        self.assertRaises(TraitError, setattr, a, 'inst', Bar)
        self.assertRaises(TraitError, setattr, a, 'inst', Bah())

    def test_default_klass(self):
        class Foo(object): pass
        class Bar(Foo): pass
        class Bah(object): pass

        class FooInstance(Instance):
            klass = Foo

        class A(HasTraits):
            inst = FooInstance(allow_none=True)

        a = A()
        self.assertTrue(a.inst is None)
        a.inst = Foo()
        self.assertTrue(isinstance(a.inst, Foo))
        a.inst = Bar()
        self.assertTrue(isinstance(a.inst, Foo))
        self.assertRaises(TraitError, setattr, a, 'inst', Foo)
        self.assertRaises(TraitError, setattr, a, 'inst', Bar)
        self.assertRaises(TraitError, setattr, a, 'inst', Bah())

    def test_unique_default_value(self):
        class Foo(object): pass
        class A(HasTraits):
            inst = Instance(Foo,(),{})

        a = A()
        b = A()
        self.assertTrue(a.inst is not b.inst)

    def test_args_kw(self):
        class Foo(object):
            def __init__(self, c): self.c = c
        class Bar(object): pass
        class Bah(object):
            def __init__(self, c, d):
                self.c = c; self.d = d

        class A(HasTraits):
            inst = Instance(Foo, (10,))
        a = A()
        self.assertEqual(a.inst.c, 10)

        class B(HasTraits):
            inst = Instance(Bah, args=(10,), kw=dict(d=20))
        b = B()
        self.assertEqual(b.inst.c, 10)
        self.assertEqual(b.inst.d, 20)

        class C(HasTraits):
            inst = Instance(Foo, allow_none=True)
        c = C()
        self.assertTrue(c.inst is None)

    def test_bad_default(self):
        class Foo(object): pass

        class A(HasTraits):
            inst = Instance(Foo)

        a = A()
        with self.assertRaises(TraitError):
            a.inst

    def test_instance(self):
        class Foo(object): pass

        def inner():
            class A(HasTraits):
                inst = Instance(Foo())

        self.assertRaises(TraitError, inner)


class TestThis(TestCase):

    def test_this_class(self):
        class Foo(HasTraits):
            this = This()

        f = Foo()
        self.assertEqual(f.this, None)
        g = Foo()
        f.this = g
        self.assertEqual(f.this, g)
        self.assertRaises(TraitError, setattr, f, 'this', 10)

    def test_this_inst(self):
        class Foo(HasTraits):
            this = This()

        f = Foo()
        f.this = Foo()
        self.assertTrue(isinstance(f.this, Foo))

    def test_subclass(self):
        class Foo(HasTraits):
            t = This()
        class Bar(Foo):
            pass
        f = Foo()
        b = Bar()
        f.t = b
        b.t = f
        self.assertEqual(f.t, b)
        self.assertEqual(b.t, f)

    def test_subclass_override(self):
        class Foo(HasTraits):
            t = This()
        class Bar(Foo):
            t = This()
        f = Foo()
        b = Bar()
        f.t = b
        self.assertEqual(f.t, b)
        self.assertRaises(TraitError, setattr, b, 't', f)

    def test_this_in_container(self):

        class Tree(HasTraits):
            value = Unicode()
            leaves = List(This())

        tree = Tree(
            value='foo',
            leaves=[Tree(value='bar'), Tree(value='buzz')]
        )

        with self.assertRaises(TraitError):
            tree.leaves = [1, 2]

class TraitTestBase(TestCase):
    """A best testing class for basic trait types."""

    def assign(self, value):
        self.obj.value = value

    def coerce(self, value):
        return value

    def test_good_values(self):
        if hasattr(self, '_good_values'):
            for value in self._good_values:
                self.assign(value)
                self.assertEqual(self.obj.value, self.coerce(value))

    def test_bad_values(self):
        if hasattr(self, '_bad_values'):
            for value in self._bad_values:
                try:
                    self.assertRaises(TraitError, self.assign, value)
                except AssertionError:
                    assert False, value

    def test_default_value(self):
        if hasattr(self, '_default_value'):
            self.assertEqual(self._default_value, self.obj.value)

    def test_allow_none(self):
        if (hasattr(self, '_bad_values') and hasattr(self, '_good_values') and
        None in self._bad_values):
            trait=self.obj.traits()['value']
            try:
                trait.allow_none = True
                self._bad_values.remove(None)
                #skip coerce. Allow None casts None to None.
                self.assign(None)
                self.assertEqual(self.obj.value,None)
                self.test_good_values()
                self.test_bad_values()
            finally:
                #tear down
                trait.allow_none = False
                self._bad_values.append(None)

    def tearDown(self):
        # restore default value after tests, if set
        if hasattr(self, '_default_value'):
            self.obj.value = self._default_value


class AnyTrait(HasTraits):

    value = Any()

class AnyTraitTest(TraitTestBase):

    obj = AnyTrait()

    _default_value = None
    _good_values   = [10.0, 'ten', u'ten', [10], {'ten': 10},(10,), None, 1j]
    _bad_values    = []

class UnionTrait(HasTraits):

    value = Union([Type(), Bool()])

class UnionTraitTest(TraitTestBase):

    obj = UnionTrait(value='ipython_genutils.ipstruct.Struct')
    _good_values = [int, float, True]
    _bad_values = [[], (0,), 1j]

class OrTrait(HasTraits):

    value = Bool() | Unicode()

class OrTraitTest(TraitTestBase):

    obj = OrTrait()
    _good_values = [True, False, 'ten']
    _bad_values = [[], (0,), 1j]

class IntTrait(HasTraits):

    value = Int(99, min=-100)

class TestInt(TraitTestBase):

    obj = IntTrait()
    _default_value = 99
    _good_values   = [10, -10]
    _bad_values    = ['ten', u'ten', [10], {'ten': 10}, (10,), None, 1j,
                      10.1, -10.1, '10L', '-10L', '10.1', '-10.1', u'10L',
                      u'-10L', u'10.1', u'-10.1',  '10', '-10', u'10', -200]
    if not six.PY3:
        _bad_values.extend([long(10), long(-10), 10*sys.maxint, -10*sys.maxint])


class CIntTrait(HasTraits):
    value = CInt('5')

class TestCInt(TraitTestBase):
    obj = CIntTrait()

    _default_value = 5
    _good_values   = ['10', '-10', u'10', u'-10', 10, 10.0, -10.0, 10.1]
    _bad_values    = ['ten', u'ten', [10], {'ten': 10},(10,),
                      None, 1j, '10.1', u'10.1']

    def coerce(self, n):
        return int(n)


class MinBoundCIntTrait(HasTraits):
    value = CInt('5', min=3)

class TestMinBoundCInt(TestCInt):
    obj = MinBoundCIntTrait()

    _default_value = 5
    _good_values   = [3, 3.0, '3']
    _bad_values    = [2.6, 2, -3, -3.0]


class LongTrait(HasTraits):

    value = Long(99 if six.PY3 else long(99))

class TestLong(TraitTestBase):

    obj = LongTrait()

    _default_value = 99 if six.PY3 else long(99)
    _good_values   = [10, -10]
    _bad_values    = ['ten', u'ten', [10], {'ten': 10},(10,),
                      None, 1j, 10.1, -10.1, '10', '-10', '10L', '-10L', '10.1',
                      '-10.1', u'10', u'-10', u'10L', u'-10L', u'10.1',
                      u'-10.1']
    if not six.PY3:
        # maxint undefined on py3, because int == long
        _good_values.extend([long(10), long(-10), 10*sys.maxint, -10*sys.maxint])
        _bad_values.extend([[long(10)], (long(10),)])

    @mark.skipif(six.PY3, reason="not relevant on py3")
    def test_cast_small(self):
        """Long casts ints to long"""
        self.obj.value = 10
        self.assertEqual(type(self.obj.value), long)


class MinBoundLongTrait(HasTraits):
    value = Long(99 if six.PY3 else long(99), min=5)

class TestMinBoundLong(TraitTestBase):
    obj = MinBoundLongTrait()

    _default_value = 99 if six.PY3 else long(99)
    _good_values   = [5, 10]
    _bad_values    = [4, -10]


class MaxBoundLongTrait(HasTraits):
    value = Long(5 if six.PY3 else long(5), max=10)

class TestMaxBoundLong(TraitTestBase):
    obj = MaxBoundLongTrait()

    _default_value = 5 if six.PY3 else long(5)
    _good_values   = [10, -2]
    _bad_values    = [11, 20]


class CLongTrait(HasTraits):
    value = CLong('5')

class TestCLong(TraitTestBase):
    obj = CLongTrait()

    _default_value = 5 if six.PY3 else long(5)
    _good_values   = ['10', '-10', u'10', u'-10', 10, 10.0, -10.0, 10.1]
    _bad_values    = ['ten', u'ten', [10], {'ten': 10},(10,),
                      None, 1j, '10.1', u'10.1']

    def coerce(self, n):
        return int(n) if six.PY3 else long(n)


class MaxBoundCLongTrait(HasTraits):
    value = CLong('5', max=10)

class TestMaxBoundCLong(TestCLong):
    obj = MaxBoundCLongTrait()

    _default_value = 5 if six.PY3 else long(5)
    _good_values   = [10, '10', 10.3]
    _bad_values    = [11.0, '11']


class IntegerTrait(HasTraits):
    value = Integer(1)

class TestInteger(TestLong):
    obj = IntegerTrait()
    _default_value = 1

    def coerce(self, n):
        return int(n)

    @mark.skipif(six.PY3, reason="not relevant on py3")
    def test_cast_small(self):
        """Integer casts small longs to int"""

        self.obj.value = long(100)
        self.assertEqual(type(self.obj.value), int)


class MinBoundIntegerTrait(HasTraits):
    value = Integer(5, min=3)

class TestMinBoundInteger(TraitTestBase):
    obj = MinBoundIntegerTrait()

    _default_value = 5
    _good_values   = 3, 20
    _bad_values    = [2, -10]


class MaxBoundIntegerTrait(HasTraits):
    value = Integer(1, max=3)

class TestMaxBoundInteger(TraitTestBase):
    obj = MaxBoundIntegerTrait()

    _default_value = 1
    _good_values   = 3, -2
    _bad_values    = [4, 10]


class FloatTrait(HasTraits):

    value = Float(99.0, max=200.0)

class TestFloat(TraitTestBase):

    obj = FloatTrait()

    _default_value = 99.0
    _good_values   = [10, -10, 10.1, -10.1]
    _bad_values    = ['ten', u'ten', [10], {'ten': 10}, (10,), None,
                      1j, '10', '-10', '10L', '-10L', '10.1', '-10.1', u'10',
                      u'-10', u'10L', u'-10L', u'10.1', u'-10.1', 201.0]
    if not six.PY3:
        _bad_values.extend([long(10), long(-10)])


class CFloatTrait(HasTraits):

    value = CFloat('99.0', max=200.0)

class TestCFloat(TraitTestBase):

    obj = CFloatTrait()

    _default_value = 99.0
    _good_values   = [10, 10.0, 10.5, '10.0', '10', '-10', '10.0', u'10']
    _bad_values    = ['ten', u'ten', [10], {'ten': 10}, (10,), None, 1j,
                      200.1, '200.1']

    def coerce(self, v):
        return float(v)


class ComplexTrait(HasTraits):

    value = Complex(99.0-99.0j)

class TestComplex(TraitTestBase):

    obj = ComplexTrait()

    _default_value = 99.0-99.0j
    _good_values   = [10, -10, 10.1, -10.1, 10j, 10+10j, 10-10j,
                      10.1j, 10.1+10.1j, 10.1-10.1j]
    _bad_values    = [u'10L', u'-10L', 'ten', [10], {'ten': 10},(10,), None]
    if not six.PY3:
        _bad_values.extend([long(10), long(-10)])


class BytesTrait(HasTraits):

    value = Bytes(b'string')

class TestBytes(TraitTestBase):

    obj = BytesTrait()

    _default_value = b'string'
    _good_values   = [b'10', b'-10', b'10L',
                      b'-10L', b'10.1', b'-10.1', b'string']
    _bad_values    = [10, -10, 10.1, -10.1, 1j, [10],
                      ['ten'],{'ten': 10},(10,), None,  u'string']
    if not six.PY3:
        _bad_values.extend([long(10), long(-10)])


class UnicodeTrait(HasTraits):

    value = Unicode(u'unicode')

class TestUnicode(TraitTestBase):

    obj = UnicodeTrait()

    _default_value = u'unicode'
    _good_values   = ['10', '-10', '10L', '-10L', '10.1',
                      '-10.1', '', u'', 'string', u'string', u"€"]
    _bad_values    = [10, -10, 10.1, -10.1, 1j,
                      [10], ['ten'], [u'ten'], {'ten': 10},(10,), None]
    if not six.PY3:
        _bad_values.extend([long(10), long(-10)])


class ObjectNameTrait(HasTraits):
    value = ObjectName("abc")

class TestObjectName(TraitTestBase):
    obj = ObjectNameTrait()

    _default_value = "abc"
    _good_values = ["a", "gh", "g9", "g_", "_G", u"a345_"]
    _bad_values = [1, "", u"€", "9g", "!", "#abc", "aj@", "a.b", "a()", "a[0]",
                                                        None, object(), object]
    if sys.version_info[0] < 3:
        _bad_values.append(u"þ")
    else:
        _good_values.append(u"þ")  # þ=1 is valid in Python 3 (PEP 3131).


class DottedObjectNameTrait(HasTraits):
    value = DottedObjectName("a.b")

class TestDottedObjectName(TraitTestBase):
    obj = DottedObjectNameTrait()

    _default_value = "a.b"
    _good_values = ["A", "y.t", "y765.__repr__", "os.path.join", u"os.path.join"]
    _bad_values = [1, u"abc.€", "_.@", ".", ".abc", "abc.", ".abc.", None]
    if sys.version_info[0] < 3:
        _bad_values.append(u"t.þ")
    else:
        _good_values.append(u"t.þ")


class TCPAddressTrait(HasTraits):
    value = TCPAddress()

class TestTCPAddress(TraitTestBase):

    obj = TCPAddressTrait()

    _default_value = ('127.0.0.1',0)
    _good_values = [('localhost',0),('192.168.0.1',1000),('www.google.com',80)]
    _bad_values = [(0,0),('localhost',10.0),('localhost',-1), None]

class ListTrait(HasTraits):

    value = List(Int())

class TestList(TraitTestBase):

    obj = ListTrait()

    _default_value = []
    _good_values = [[], [1], list(range(10)), (1,2)]
    _bad_values = [10, [1,'a'], 'a']
    
    def coerce(self, value):
        if value is not None:
            value = list(value)
        return value

class Foo(object):
    pass

class NoneInstanceListTrait(HasTraits):
    
    value = List(Instance(Foo))

class TestNoneInstanceList(TraitTestBase):

    obj = NoneInstanceListTrait()

    _default_value = []
    _good_values = [[Foo(), Foo()], []]
    _bad_values = [[None], [Foo(), None]]


class InstanceListTrait(HasTraits):

    value = List(Instance(__name__+'.Foo'))

class TestInstanceList(TraitTestBase):

    obj = InstanceListTrait()

    def test_klass(self):
        """Test that the instance klass is properly assigned."""
        self.assertIs(self.obj.traits()['value']._trait.klass, Foo)

    _default_value = []
    _good_values = [[Foo(), Foo()], []]
    _bad_values = [['1', 2,], '1', [Foo], None]

class UnionListTrait(HasTraits):

    value = List(Int() | Bool())

class TestUnionListTrait(HasTraits):

    obj = UnionListTrait()

    _default_value = []
    _good_values = [[True, 1], [False, True]]
    _bad_values = [[1, 'True'], False]


class LenListTrait(HasTraits):

    value = List(Int(), [0], minlen=1, maxlen=2)

class TestLenList(TraitTestBase):

    obj = LenListTrait()

    _default_value = [0]
    _good_values = [[1], [1,2], (1,2)]
    _bad_values = [10, [1,'a'], 'a', [], list(range(3))]

    def coerce(self, value):
        if value is not None:
            value = list(value)
        return value

class TupleTrait(HasTraits):

    value = Tuple(Int(allow_none=True), default_value=(1,))

class TestTupleTrait(TraitTestBase):

    obj = TupleTrait()

    _default_value = (1,) 
    _good_values = [(1,), (0,), [1]]
    _bad_values = [10, (1, 2), ('a'), (), None]

    def coerce(self, value):
        if value is not None:
            value = tuple(value)
        return value

    def test_invalid_args(self):
        self.assertRaises(TypeError, Tuple, 5)
        self.assertRaises(TypeError, Tuple, default_value='hello')
        t = Tuple(Int(), CBytes(), default_value=(1,5))

class LooseTupleTrait(HasTraits):

    value = Tuple((1,2,3))

class TestLooseTupleTrait(TraitTestBase):

    obj = LooseTupleTrait()

    _default_value = (1,2,3)
    _good_values = [(1,), [1], (0,), tuple(range(5)), tuple('hello'), ('a',5), ()]
    _bad_values = [10, 'hello', {}, None]

    def coerce(self, value):
        if value is not None:
            value = tuple(value)
        return value

    def test_invalid_args(self):
        self.assertRaises(TypeError, Tuple, 5)
        self.assertRaises(TypeError, Tuple, default_value='hello')
        t = Tuple(Int(), CBytes(), default_value=(1,5))


class MultiTupleTrait(HasTraits):

    value = Tuple(Int(), Bytes(), default_value=[99,b'bottles'])

class TestMultiTuple(TraitTestBase):

    obj = MultiTupleTrait()

    _default_value = (99,b'bottles')
    _good_values = [(1,b'a'), (2,b'b')]
    _bad_values = ((),10, b'a', (1,b'a',3), (b'a',1), (1, u'a'))

class CRegExpTrait(HasTraits):

    value = CRegExp(r'')

class TestCRegExp(TraitTestBase):

    def coerce(self, value):
        return re.compile(value)

    obj = CRegExpTrait()

    _default_value = re.compile(r'')
    _good_values = [r'\d+', re.compile(r'\d+')]
    _bad_values = ['(', None, ()]

class DictTrait(HasTraits):
    value = Dict()

def test_dict_assignment():
    d = dict()
    c = DictTrait()
    c.value = d
    d['a'] = 5
    assert d == c.value  
    assert c.value is d


class UniformlyValidatedDictTrait(HasTraits):

    value = Dict(trait=Unicode(),
                 default_value={'foo': '1'})


class TestInstanceUniformlyValidatedDict(TraitTestBase):

    obj = UniformlyValidatedDictTrait()

    _default_value = {'foo': '1'}
    _good_values = [{'foo': '0', 'bar': '1'}]
    _bad_values = [{'foo': 0, 'bar': '1'}]


class KeyValidatedDictTrait(HasTraits):

    value = Dict(traits={'foo': Int()},
                 default_value={'foo': 1})


class TestInstanceKeyValidatedDict(TraitTestBase):

    obj = KeyValidatedDictTrait()

    _default_value = {'foo': 1}
    _good_values = [{'foo': 0, 'bar': '1'}, {'foo': 0, 'bar': 1}]
    _bad_values = [{'foo': '0', 'bar': '1'}]


class FullyValidatedDictTrait(HasTraits):

    value = Dict(trait=Unicode(),
                 traits={'foo': Int()},
                 default_value={'foo': 1})


class TestInstanceFullyValidatedDict(TraitTestBase):

    obj = FullyValidatedDictTrait()

    _default_value = {'foo': 1}
    _good_values = [{'foo': 0, 'bar': '1'}, {'foo': 1, 'bar': '2'}]
    _bad_values = [{'foo': 0, 'bar': 1}, {'foo': '0', 'bar': '1'}]


def test_dict_default_value():
    """Check that the `{}` default value of the Dict traitlet constructor is
    actually copied."""

    class Foo(HasTraits):
        d1 = Dict()
        d2 = Dict()

    foo = Foo()
    assert foo.d1 == {}
    assert foo.d2 ==  {}
    assert foo.d1 is not foo.d2


class TestValidationHook(TestCase):

    def test_parity_trait(self):
        """Verify that the early validation hook is effective"""

        class Parity(HasTraits):

            value = Int(0)
            parity = Enum(['odd', 'even'], default_value='even')

            @validate('value')
            def _value_validate(self, proposal):
                value = proposal['value']
                if self.parity == 'even' and value % 2:
                    raise TraitError('Expected an even number')
                if self.parity == 'odd' and (value % 2 == 0):
                    raise TraitError('Expected an odd number')
                return value
        
        u = Parity()
        u.parity = 'odd'
        u.value = 1  # OK
        with self.assertRaises(TraitError):
            u.value = 2  # Trait Error

        u.parity = 'even'
        u.value = 2  # OK

    def test_multiple_validate(self):
        """Verify that we can register the same validator to multiple names"""

        class OddEven(HasTraits):

            odd = Int(1)
            even = Int(0)

            @validate('odd', 'even')
            def check_valid(self, proposal):
                if proposal['trait'].name == 'odd' and not proposal['value'] % 2:
                    raise TraitError('odd should be odd')
                if proposal['trait'].name == 'even' and proposal['value'] % 2:
                    raise TraitError('even should be even')

        u = OddEven()
        u.odd = 3  # OK
        with self.assertRaises(TraitError):
            u.odd = 2  # Trait Error

        u.even = 2  # OK
        with self.assertRaises(TraitError):
            u.even = 3  # Trait Error



class TestLink(TestCase):

    def test_connect_same(self):
        """Verify two traitlets of the same type can be linked together using link."""

        # Create two simple classes with Int traitlets.
        class A(HasTraits):
            value = Int()
        a = A(value=9)
        b = A(value=8)

        # Conenct the two classes.
        c = link((a, 'value'), (b, 'value'))

        # Make sure the values are the same at the point of linking.
        self.assertEqual(a.value, b.value)

        # Change one of the values to make sure they stay in sync.
        a.value = 5
        self.assertEqual(a.value, b.value)
        b.value = 6
        self.assertEqual(a.value, b.value)

    def test_link_different(self):
        """Verify two traitlets of different types can be linked together using link."""

        # Create two simple classes with Int traitlets.
        class A(HasTraits):
            value = Int()
        class B(HasTraits):
            count = Int()
        a = A(value=9)
        b = B(count=8)

        # Conenct the two classes.
        c = link((a, 'value'), (b, 'count'))

        # Make sure the values are the same at the point of linking.
        self.assertEqual(a.value, b.count)

        # Change one of the values to make sure they stay in sync.
        a.value = 5
        self.assertEqual(a.value, b.count)
        b.count = 4
        self.assertEqual(a.value, b.count)

    def test_unlink(self):
        """Verify two linked traitlets can be unlinked."""

        # Create two simple classes with Int traitlets.
        class A(HasTraits):
            value = Int()
        a = A(value=9)
        b = A(value=8)

        # Connect the two classes.
        c = link((a, 'value'), (b, 'value'))
        a.value = 4
        c.unlink()

        # Change one of the values to make sure they don't stay in sync.
        a.value = 5
        self.assertNotEqual(a.value, b.value)

    def test_callbacks(self):
        """Verify two linked traitlets have their callbacks called once."""

        # Create two simple classes with Int traitlets.
        class A(HasTraits):
            value = Int()
        class B(HasTraits):
            count = Int()
        a = A(value=9)
        b = B(count=8)
        
        # Register callbacks that count.
        callback_count = []
        def a_callback(name, old, new):
            callback_count.append('a')
        a.on_trait_change(a_callback, 'value')
        def b_callback(name, old, new):
            callback_count.append('b')
        b.on_trait_change(b_callback, 'count')

        # Connect the two classes.
        c = link((a, 'value'), (b, 'count'))

        # Make sure b's count was set to a's value once.
        self.assertEqual(''.join(callback_count), 'b')
        del callback_count[:]

        # Make sure a's value was set to b's count once.
        b.count = 5
        self.assertEqual(''.join(callback_count), 'ba')
        del callback_count[:]

        # Make sure b's count was set to a's value once.
        a.value = 4
        self.assertEqual(''.join(callback_count), 'ab')
        del callback_count[:]

class TestDirectionalLink(TestCase):
    def test_connect_same(self):
        """Verify two traitlets of the same type can be linked together using directional_link."""

        # Create two simple classes with Int traitlets.
        class A(HasTraits):
            value = Int()
        a = A(value=9)
        b = A(value=8)

        # Conenct the two classes.
        c = directional_link((a, 'value'), (b, 'value'))

        # Make sure the values are the same at the point of linking.
        self.assertEqual(a.value, b.value)

        # Change one the value of the source and check that it synchronizes the target.
        a.value = 5
        self.assertEqual(b.value, 5)
        # Change one the value of the target and check that it has no impact on the source
        b.value = 6
        self.assertEqual(a.value, 5)

    def test_tranform(self):
        """Test transform link."""

        # Create two simple classes with Int traitlets.
        class A(HasTraits):
            value = Int()
        a = A(value=9)
        b = A(value=8)

        # Conenct the two classes.
        c = directional_link((a, 'value'), (b, 'value'), lambda x: 2 * x)

        # Make sure the values are correct at the point of linking.
        self.assertEqual(b.value, 2 * a.value)

        # Change one the value of the source and check that it modifies the target.
        a.value = 5
        self.assertEqual(b.value, 10)
        # Change one the value of the target and check that it has no impact on the source
        b.value = 6
        self.assertEqual(a.value, 5)

    def test_link_different(self):
        """Verify two traitlets of different types can be linked together using link."""

        # Create two simple classes with Int traitlets.
        class A(HasTraits):
            value = Int()
        class B(HasTraits):
            count = Int()
        a = A(value=9)
        b = B(count=8)

        # Conenct the two classes.
        c = directional_link((a, 'value'), (b, 'count'))

        # Make sure the values are the same at the point of linking.
        self.assertEqual(a.value, b.count)

        # Change one the value of the source and check that it synchronizes the target.
        a.value = 5
        self.assertEqual(b.count, 5)
        # Change one the value of the target and check that it has no impact on the source
        b.value = 6
        self.assertEqual(a.value, 5)

    def test_unlink(self):
        """Verify two linked traitlets can be unlinked."""

        # Create two simple classes with Int traitlets.
        class A(HasTraits):
            value = Int()
        a = A(value=9)
        b = A(value=8)

        # Connect the two classes.
        c = directional_link((a, 'value'), (b, 'value'))
        a.value = 4
        c.unlink()

        # Change one of the values to make sure they don't stay in sync.
        a.value = 5
        self.assertNotEqual(a.value, b.value)

class Pickleable(HasTraits):

    i = Int()
    @observe('i')
    def _i_changed(self, change): pass
    @validate('i')
    def _i_validate(self, commit):
        return commit['value']

    j = Int()
    
    def __init__(self):
        with self.hold_trait_notifications():
            self.i = 1
        self.on_trait_change(self._i_changed, 'i')

def test_pickle_hastraits():
    c = Pickleable()
    for protocol in range(pickle.HIGHEST_PROTOCOL + 1):
        p = pickle.dumps(c, protocol)
        c2 = pickle.loads(p)
        assert c2.i == c.i
        assert c2.j == c.j

    c.i = 5
    for protocol in range(pickle.HIGHEST_PROTOCOL + 1):
        p = pickle.dumps(c, protocol)
        c2 = pickle.loads(p)
        assert c2.i == c.i
        assert c2.j == c.j


def test_hold_trait_notifications():
    changes = []

    class Test(HasTraits):
        a = Integer(0)
        b = Integer(0)

        def _a_changed(self, name, old, new):
            changes.append((old, new))

        def _b_validate(self, value, trait):
            if value != 0:
                raise TraitError('Only 0 is a valid value')
            return value

    # Test context manager and nesting
    t = Test()
    with t.hold_trait_notifications():
        with t.hold_trait_notifications():
            t.a = 1
            assert t.a == 1
            assert changes == []
        t.a = 2
        assert t.a == 2
        with t.hold_trait_notifications():
            t.a = 3
            assert t.a == 3
            assert changes == []
            t.a = 4
            assert t.a == 4
            assert changes == []
        t.a = 4
        assert t.a == 4
        assert changes == []

    assert changes == [(0, 4)]
    # Test roll-back
    try:
         with t.hold_trait_notifications():
             t.b = 1  # raises a Trait error
    except:
        pass
    assert t.b == 0


class RollBack(HasTraits):
    bar = Int()
    def _bar_validate(self, value, trait):
        if value:
            raise TraitError('foobar')
        return value


class TestRollback(TestCase):

    def test_roll_back(self):

        def assign_rollback():
            RollBack(bar=1)

        self.assertRaises(TraitError, assign_rollback)


class CacheModification(HasTraits):
    foo = Int()
    bar = Int()

    def _bar_validate(self, value, trait):
        self.foo = value
        return value

    def _foo_validate(self, value, trait):
        self.bar = value
        return value


def test_cache_modification():
    CacheModification(foo=1)
    CacheModification(bar=1)


class OrderTraits(HasTraits):
    notified = Dict()
    
    a = Unicode()
    b = Unicode()
    c = Unicode()
    d = Unicode()
    e = Unicode()
    f = Unicode()
    g = Unicode()
    h = Unicode()
    i = Unicode()
    j = Unicode()
    k = Unicode()
    l = Unicode()
    
    def _notify(self, name, old, new):
        """check the value of all traits when each trait change is triggered
        
        This verifies that the values are not sensitive
        to dict ordering when loaded from kwargs
        """
        # check the value of the other traits
        # when a given trait change notification fires
        self.notified[name] = {
            c: getattr(self, c) for c in 'abcdefghijkl'
        }
    
    def __init__(self, **kwargs):
        self.on_trait_change(self._notify)
        super(OrderTraits, self).__init__(**kwargs)

def test_notification_order():
    d = {c:c for c in 'abcdefghijkl'}
    obj = OrderTraits()
    assert obj.notified == {}
    obj = OrderTraits(**d)
    notifications = {
        c: d for c in 'abcdefghijkl'
    }
    assert obj.notified == notifications



###
# Traits for Forward Declaration Tests
###
class ForwardDeclaredInstanceTrait(HasTraits):

    value = ForwardDeclaredInstance('ForwardDeclaredBar', allow_none=True)

class ForwardDeclaredTypeTrait(HasTraits):

    value = ForwardDeclaredType('ForwardDeclaredBar', allow_none=True)

class ForwardDeclaredInstanceListTrait(HasTraits):

    value = List(ForwardDeclaredInstance('ForwardDeclaredBar'))

class ForwardDeclaredTypeListTrait(HasTraits):

    value = List(ForwardDeclaredType('ForwardDeclaredBar'))
###
# End Traits for Forward Declaration Tests
###

###
# Classes for Forward Declaration Tests
###
class ForwardDeclaredBar(object):
    pass

class ForwardDeclaredBarSub(ForwardDeclaredBar):
    pass
###
# End Classes for Forward Declaration Tests
###

###
# Forward Declaration Tests
###
class TestForwardDeclaredInstanceTrait(TraitTestBase):

    obj = ForwardDeclaredInstanceTrait()
    _default_value = None
    _good_values = [None, ForwardDeclaredBar(), ForwardDeclaredBarSub()]
    _bad_values = ['foo', 3, ForwardDeclaredBar, ForwardDeclaredBarSub]

class TestForwardDeclaredTypeTrait(TraitTestBase):

    obj = ForwardDeclaredTypeTrait()
    _default_value = None
    _good_values = [None, ForwardDeclaredBar, ForwardDeclaredBarSub]
    _bad_values = ['foo', 3, ForwardDeclaredBar(), ForwardDeclaredBarSub()]

class TestForwardDeclaredInstanceList(TraitTestBase):

    obj = ForwardDeclaredInstanceListTrait()

    def test_klass(self):
        """Test that the instance klass is properly assigned."""
        self.assertIs(self.obj.traits()['value']._trait.klass, ForwardDeclaredBar)

    _default_value = []
    _good_values = [
        [ForwardDeclaredBar(), ForwardDeclaredBarSub()],
        [],
    ]
    _bad_values = [
        ForwardDeclaredBar(),
        [ForwardDeclaredBar(), 3, None],
        '1',
        # Note that this is the type, not an instance.
        [ForwardDeclaredBar],
        [None],
        None,
    ]

class TestForwardDeclaredTypeList(TraitTestBase):

    obj = ForwardDeclaredTypeListTrait()

    def test_klass(self):
        """Test that the instance klass is properly assigned."""
        self.assertIs(self.obj.traits()['value']._trait.klass, ForwardDeclaredBar)

    _default_value = []
    _good_values = [
        [ForwardDeclaredBar, ForwardDeclaredBarSub],
        [],
    ]
    _bad_values = [
        ForwardDeclaredBar,
        [ForwardDeclaredBar, 3],
        '1',
        # Note that this is an instance, not the type.
        [ForwardDeclaredBar()],
        [None],
        None,
    ]
###
# End Forward Declaration Tests
###

class TestDynamicTraits(TestCase):

    def setUp(self):
        self._notify1 = []

    def notify1(self, name, old, new):
        self._notify1.append((name, old, new))

    def test_notify_all(self):

        class A(HasTraits):
            pass

        a = A()
        self.assertTrue(not hasattr(a, 'x'))
        self.assertTrue(not hasattr(a, 'y'))

        # Dynamically add trait x.
        a.add_traits(x=Int())
        self.assertTrue(hasattr(a, 'x'))
        self.assertTrue(isinstance(a, (A, )))

        # Dynamically add trait y.
        a.add_traits(y=Float())
        self.assertTrue(hasattr(a, 'y'))
        self.assertTrue(isinstance(a, (A, )))
        self.assertEqual(a.__class__.__name__, A.__name__)

        # Create a new instance and verify that x and y
        # aren't defined.
        b = A()
        self.assertTrue(not hasattr(b, 'x'))
        self.assertTrue(not hasattr(b, 'y'))

        # Verify that notification works like normal.
        a.on_trait_change(self.notify1)
        a.x = 0
        self.assertEqual(len(self._notify1), 0)
        a.y = 0.0
        self.assertEqual(len(self._notify1), 0)
        a.x = 10
        self.assertTrue(('x', 0, 10) in self._notify1)
        a.y = 10.0
        self.assertTrue(('y', 0.0, 10.0) in self._notify1)
        self.assertRaises(TraitError, setattr, a, 'x', 'bad string')
        self.assertRaises(TraitError, setattr, a, 'y', 'bad string')
        self._notify1 = []
        a.on_trait_change(self.notify1, remove=True)
        a.x = 20
        a.y = 20.0
        self.assertEqual(len(self._notify1), 0)


def test_enum_no_default():
    class C(HasTraits):
        t = Enum(['a', 'b'])

    c = C()
    c.t  = 'a'
    assert c.t == 'a'

    c = C()

    with pytest.raises(TraitError):
        t = c.t

    c = C(t='b')
    assert c.t == 'b'


def test_default_value_repr():
    class C(HasTraits):
        t = Type('traitlets.HasTraits')
        t2 = Type(HasTraits)
        n = Integer(0)
        lis = List()
        d = Dict()
    
    assert C.t.default_value_repr() == "'traitlets.HasTraits'"
    assert C.t2.default_value_repr() == "'traitlets.traitlets.HasTraits'"
    assert C.n.default_value_repr() == '0'
    assert C.lis.default_value_repr() == '[]'
    assert C.d.default_value_repr() == '{}'


class TransitionalClass(HasTraits):
    
    d = Any()
    @default('d')
    def _d_default(self):
        return TransitionalClass

    parent_super = False
    calls_super = Integer(0)
    
    @default('calls_super')
    def _calls_super_default(self):
        return -1
    
    @observe('calls_super')
    @observe_compat
    def _calls_super_changed(self, change):
        self.parent_super = change
    
    parent_override = False
    overrides = Integer(0)
    
    @observe('overrides')
    @observe_compat
    def _overrides_changed(self, change):
        self.parent_override = change


class SubClass(TransitionalClass):
    def _d_default(self):
        return SubClass
    
    subclass_super = False
    def _calls_super_changed(self, name, old, new):
        self.subclass_super = True
        super(SubClass, self)._calls_super_changed(name, old, new)
    
    subclass_override = False
    def _overrides_changed(self, name, old, new):
        self.subclass_override = True


def test_subclass_compat():
    obj = SubClass()
    obj.calls_super = 5
    assert obj.parent_super
    assert obj.subclass_super
    obj.overrides = 5
    assert obj.subclass_override
    assert not obj.parent_override
    assert obj.d is SubClass


class DefinesHandler(HasTraits):
    parent_called = False
    
    trait = Integer()
    @observe('trait')
    def handler(self, change):
        self.parent_called = True


class OverridesHandler(DefinesHandler):
    child_called = False
    
    @observe('trait')
    def handler(self, change):
        self.child_called = True


def test_subclass_override_observer():
    obj = OverridesHandler()
    obj.trait = 5
    assert obj.child_called
    assert not obj.parent_called


class DoesntRegisterHandler(DefinesHandler):
    child_called = False
    
    def handler(self, change):
        self.child_called = True


def test_subclass_override_not_registered():
    """Subclass that overrides observer and doesn't re-register unregisters both"""
    obj = DoesntRegisterHandler()
    obj.trait = 5
    assert not obj.child_called
    assert not obj.parent_called


class AddsHandler(DefinesHandler):
    child_called = False
    
    @observe('trait')
    def child_handler(self, change):
        self.child_called = True

def test_subclass_add_observer():
    obj = AddsHandler()
    obj.trait = 5
    assert obj.child_called
    assert obj.parent_called


def test_observe_iterables():

    class C(HasTraits):
        i = Integer()
        s = Unicode()

    c = C()
    recorded = {}
    def record(change):
        recorded['change'] = change

    # observe with names=set
    c.observe(record, names={'i', 's'})
    c.i = 5
    assert recorded['change'].name == 'i'
    assert recorded['change'].new == 5
    c.s = 'hi'
    assert recorded['change'].name == 's'
    assert recorded['change'].new == 'hi'

    # observe with names=custom container with iter, contains
    class MyContainer(object):
        def __init__(self, container):
            self.container = container

        def __iter__(self):
            return iter(self.container)

        def __contains__(self, key):
            return key in self.container

    c.observe(record, names=MyContainer({'i', 's'}))
    c.i = 10
    assert recorded['change'].name == 'i'
    assert recorded['change'].new == 10
    c.s = 'ok'
    assert recorded['change'].name == 's'
    assert recorded['change'].new == 'ok'


def test_super_args():
    class SuperRecorder(object):
        def __init__(self, *args, **kwargs):
            self.super_args = args
            self.super_kwargs = kwargs
    
    class SuperHasTraits(HasTraits, SuperRecorder):
        i = Integer()
    
    obj = SuperHasTraits('a1', 'a2', b=10, i=5, c='x')
    assert obj.i ==  5
    assert not hasattr(obj, 'b')
    assert not hasattr(obj, 'c')
    assert obj.super_args == ('a1' ,  'a2')
    assert obj.super_kwargs == {'b': 10 , 'c': 'x'}

def test_super_bad_args():
    class SuperHasTraits(HasTraits):
        a = Integer()
    
    if sys.version_info < (3,):
        # Legacy Python, object.__init__ warns itself, instead of raising
        w = ['object.__init__']
    else:
        w = ["Passing unrecoginized arguments"]
    with expected_warnings(w):
        obj = SuperHasTraits(a=1, b=2)
    assert obj.a ==  1 
    assert not hasattr(obj, 'b')

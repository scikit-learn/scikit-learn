"""
Filters decide whether something is active or not (they decide about a boolean
state). This is used to enable/disable features, like key bindings, parts of
the layout and other stuff. For instance, we could have a `HasSearch` filter
attached to some part of the layout, in order to show that part of the user
interface only while the user is searching.

Filters are made to avoid having to attach callbacks to all event in order to
propagate state. However, they are lazy, they don't automatically propagate the
state of what they are observing. Only when a filter is called (it's actually a
callable), it will calculate its value. So, its not really reactive
programming, but it's made to fit for this framework.

One class of filters observe a `CommandLineInterface` instance. However, they
are not attached to such an instance. (We have to pass this instance to the
filter when calling it.) The reason for this is to allow declarative
programming: for key bindings, we can attach a filter to a key binding without
knowing yet which `CommandLineInterface` instance it will observe in the end.
Examples are `HasSearch` or `IsExiting`.

Another class of filters doesn't take anything as input. And a third class of
filters are universal, for instance `Always` and `Never`.
It is impossible to mix the first and the second class, because that would mean
mixing filters with a different signature.

Filters can be chained using ``&`` and ``|`` operations, and inverted using the
``~`` operator, for instance::

    filter = HasFocus('default') & ~ HasSelection()
"""
from __future__ import unicode_literals

from .base import *
from .cli import *
from .types import *
from .utils import *

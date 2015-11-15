# -*- coding: utf8 -*-
# Author: Alain Pena <alain.pena.r@gmail.com>
#
# License: BSD 3 clause

import operator
import numpy as np

try:
    basestring
except NameError:
    basestring = (str, bytes)

__all__ = [
    "CheckFalseError",
    "check_all_in",
    "check_is_binary_01",
    "check_is_in",
    "check_is_in_list",
    "check_is_in_range",
    "check_is_iterable",
    "check_is_list_or_tuple",
]


class CheckFalseError(Exception):

    """
    A user-made exception.
    """

    def __init__(self, value=''):
        self.value = value

    def __str__(self):
        return repr(self.value)


def __get_wrapper(launch_exception=False,
                  exception_type=CheckFalseError):
    """
    Wrapper to encapsulate the exception arguments.

    Parameters
    ----------
    launch_exception : bool
        Whether or not an exception should be raised
        by the function returned if the check is `False`.

    exception_type : class
        The exception to raise by the returned function
        if the check is `False`.

    Returns
    -------
    Function(val = False, msg = '')
        A function raising an error (or not) if the value is `False`.
    """
    def __wrapper(val=False, msg=''):
        if (launch_exception and not val):
            raise exception_type(str(msg))
        else:
            return val
    return __wrapper


def check_is_binary_01(obj, launch_exception=False,
                       exception_type=CheckFalseError):
    """
    Checks if an object is an iterable containing only `0` and `1`.

    Parameters
    ----------
    obj : any
        The object to test.

    launch_exception : bool
        Whether or not an exception should be raised
        if the check is `False`.

    exception_type : class
        The exception to raise if the check is `False`.

    Returns
    -------
    bool
        The value of the check.

    Raises
    ------
    exception_type
        In case the check is `False` and the function should launch
        an exception (via the `launch_exception` argument).

        The exception class is the one of the `exception_type`
        argument.
    """
    wrapper = __get_wrapper(launch_exception=launch_exception,
                            exception_type=exception_type)
    check = check_is_iterable(obj,
                              launch_exception=launch_exception,
                              exception_type=CheckFalseError)
    if not check:
        return False
    try:
        uniques = np.unique(obj)
    except TypeError as e:
        return wrapper(False, e)
    contains0 = 0 in uniques
    contains1 = 1 in uniques
    containsboth = contains0 and contains1 and (len(uniques) == 2)
    containsonlyone = ((len(uniques) == 1) and
                       (contains0 or contains1))
    return (wrapper(len(uniques) <= 2,
                    msg=('Error: the array is not binary but has '
                         '{0} different values.').format(
        str(len(uniques)))) and
        wrapper(containsboth or containsonlyone,
                msg=('Error: the array does not contain only '
                     '0 and 1 but {0}.').format(str(uniques))))


def check_is_in_list(obj, reference, launch_exception=False,
                     exception_type=CheckFalseError):
    """
    Checks if an object is in a container and the container is
    not `None`.

    Parameters
    ----------
    obj : any
        The object to test.

    reference : iterable
        The container.

    launch_exception : bool
        Whether or not an exception should be raised
        if the check is `False`.

    exception_type : class
        The exception to raise if the check is `False`.

    msg : any
        The message [`str(msg)`] that should be raised with
        the exception.

    Returns
    -------
    bool
        The value of the check.

    Raises
    ------
    exception_type
        In case the check is `False` and the function should launch
        an exception (via the `launch_exception` argument).

        The exception class is the one of the `exception_type`
        argument.
    """
    wrapper = __get_wrapper(launch_exception=launch_exception,
                            exception_type=exception_type)
    return (wrapper(reference is not None,
                    msg='Error: the reference is None.') and
            wrapper(obj in reference,
                    msg=('Error: the object ({0}) is not in list.'
                         .format(str(obj)))))


def check_is_in(obj, maybe_collection, stringed=False,
                launch_exception=False, exception_type=CheckFalseError):
    """
    Checks if an object is in a container or, if the so-called
    container is not one, if the object is this later.

    Parameters
    ----------
    obj : any
        The object to test.

    maybe_collection : iterable or any
        The "container".

        If this object is iterable (string allowed), the check is
        `True` if `obj` is in this object.

        If this object is not iterable, the check is `True` if
        `obj` == this object.

    stringed : bool
        Specifies what to do when the container is a string.

        When `stringed` is `False`, tests `obj in maybe_collection`
        (fails if `obj` is not a string), while if `stringed` is `True`,
        tests `str(obj) in maybe_collection`.

    launch_exception : bool
        Whether or not an exception should be raised
        if the check is `False`.

    exception_type : class
        The exception to raise if the check is `False`.

    Returns
    -------
    bool
        The value of the check.

    Raises
    ------
    exception_type
        In case the check is `False` and the function should launch
        an exception (via the `launch_exception` argument).

        The exception class is the one of the `exception_type`
        argument.
    """
    wrapper = __get_wrapper(launch_exception=launch_exception,
                            exception_type=exception_type)
    if (maybe_collection is None):
        return wrapper(obj is None,
                       msg=('Error: the collection is None'
                            ' but not the object.'))
    elif (check_is_iterable(maybe_collection,
                            string_allowed=False)):
        return wrapper(obj in maybe_collection,
                       msg=('Error: the object ({0}) is not in the'
                            ' collection.'.format(str(obj))))
    elif (check_is_iterable(maybe_collection,
                            string_allowed=True)):
        if (stringed):
            return wrapper(str(obj) in maybe_collection,
                           msg=('Error: the object ({0}) is not in the'
                                ' collection.'.format(str(obj))))
        try:
            i = obj in maybe_collection
        except TypeError:
            return wrapper(False,
                           msg=('Error: the object ({0}) is not in the'
                                ' collection.'.format(str(obj))))
        else:
            return wrapper(i,
                           msg=('Error: the object ({0}) is not in the'
                                ' collection.'.format(str(obj))))
    else:
        return wrapper(obj == maybe_collection,
                       msg=('Error: the object ({0}) is not equal '
                            'to the collection.'.format(str(obj))))


def check_all_in(first, second, both=False,
                 first_in_second=False, second_in_first=False,
                 launch_exception=False,
                 exception_type=CheckFalseError):
    """
    Checks that all objects in a container are in the other one.
    Returns True in case
    `both` == `first_in_second` == `second_in_first` == `False`.

    Parameters
    ----------
    first : any
        The first object to test.

    second : any
        The second object to test.

    both : bool
        Whether or not the check should be done as
            `first` :math:`\\subseteq` `second` AND
            `second` :math:`\\subseteq` `first`

        This is equivalent to
            `first_in_second` = `second_in_first` = `True`.

    first_in_second : bool
        Whether or not the check should be done as
            `first` :math:`\\subseteq` `second`

    second_in_first : bool
        Whether or not the check should be done as
            `second` :math:`\\subseteq` `first`

    launch_exception : bool
        Whether or not an exception should be raised
        if the check is `False`.

    exception_type : class
        The exception to raise if the check is `False`.

    Returns
    -------
    bool
        The value of the check.

    Raises
    ------
    exception_type
        In case the check is `False` and the function should launch
        an exception (via the `launch_exception` argument).

        The exception class is the one of the `exception_type`
        argument.
    """
    first_ = first
    second_ = second
    if (both or first_in_second):
        if not check_is_iterable(first_):
            first_ = (first,)
        for elem in first_:
            if not check_is_in(elem, second,
                               launch_exception=launch_exception,
                               exception_type=exception_type):
                return False
    if (both or second_in_first):
        if not check_is_iterable(second_):
            second_ = (second,)
        for elem in second_:
            if not check_is_in(elem, first,
                               launch_exception=launch_exception,
                               exception_type=exception_type):
                return False
    return True


def check_is_list_or_tuple(elem=None, launch_exception=False,
                           length_req=None,
                           exception_type=CheckFalseError):
    """
    Checks that the object is a `list` or a `tuple`
    of correct length.

    Parameters
    ----------
    elem : any
        The object to check

    launch_exception : bool
        Whether or not an exception should be raised
        if the check is `False`.

    length_req : Function(x) or `None`
        If `None`, no check on the length will be performed.

        If `not None`, the function should return `True`
        when the length of the `list` or `tuple` (which is passed
        to this function) meets the requirement.

    exception_type : class
        The exception to raise if the check is `False`.

    Returns
    -------
    bool
        The value of the check.

    Raises
    ------
    exception_type
        In case the check is `False` and the function should launch
        an exception (via the `launch_exception` argument).

        The exception class is the one of the `exception_type`
        argument.
    """
    wrapper = __get_wrapper(launch_exception=launch_exception,
                            exception_type=exception_type)
    return (wrapper(elem is not None,
                    msg='Error: the element argument is None.') and
            wrapper((type(elem) is tuple) or (type(elem) is list),
                    msg=('Error: the element argument is not a '
                         'list or a tuple but a {0}.'
                         ).format(str(type(elem)))) and
            wrapper((length_req is None) or length_req(len(elem)),
                    msg=('Error: the element argument is a list '
                         'or a tuple of incorrect size: {0}.'
                         ).format(str(len(elem)))))


def check_is_iterable(elem=None, string_allowed=False,
                      launch_exception=False,
                      exception_type=CheckFalseError):
    """
    Checks that the object is iterable.

    Parameters
    ----------
    elem : any
        The object to check.

    string_allowed : bool
        If set to `True`, just checks that the object is iterable.

        If set to `False`, checks if the object is iterable but not
        a `basestring` instance.

    launch_exception : bool
        Whether or not an exception should be raised
        if the check is `False`.

    exception_type : class
        The exception to raise if the check is `False`.

    Returns
    -------
    bool
        The value of the check.

    Raises
    ------
    exception_type
        In case the check is `False` and the function should launch
        an exception (via the `launch_exception` argument).

        The exception class is the one of the `exception_type`
        argument.
    """
    wrapper = __get_wrapper(launch_exception=launch_exception,
                            exception_type=exception_type)

    def is_iterable(elem):
        try:
            iter(elem)
        except TypeError:
            return False
        else:
            return True
    return (wrapper(elem is not None,
                    msg='Error: the element argument is None.') and
            wrapper(is_iterable(elem),
                    msg=('Error: the element argument '
                         'is not iterable.')) and
            wrapper(string_allowed or
                    not isinstance(elem, basestring),
                    msg='Error: the element argument is a string.'))


def check_is_in_range(elem=0.0, lower_bound=None,
                      higher_bound=None, low_exclusive=False,
                      high_exclusive=True, launch_exception=False,
                      msg="Error: the object is not in range.",
                      exception_type=CheckFalseError):
    """
    Checks if the object is in range.

    Parameters
    ----------
    elem : any
        The object to check.

        Should support `operator.lt` and `operator.le`.

    lower_bound : any
        If set to `None`, no test is done for the lower bound.

        If `not None`, checks that
        `operator.l*`(`lower_bound`, `elem`) evaluates to `True`

        Should support `operator.lt` and `operator.le` or be `None`.

    higher_bound : any
        If set to `None`, no test is done for the higher bound.

        If `not None`, checks that
        `operator.l*`(`elem`, `higher_bound`) evaluates to `True`

        Should support `operator.lt` and `operator.le` or be `None`.

    low_exclusive : bool
        For the `lower_bound`.
        If `True`, the function `operator.lt` will be used,
        `operator.le` otherwise.

    high_exclusive : bool
        For the `higher_bound`.
        If `True`, the function `operator.lt` will be used,
        `operator.le` otherwise.

    launch_exception : bool
        Whether or not an exception should be raised
        if the check is `False`.

    msg : any
        The message [`str(msg)`] that should be raised
        with the exception.

    exception_type : class
        The exception to raise if the check is `False`.

    Returns
    -------
    bool
        The value of the check.

    Raises
    ------
    exception_type
        In case the check is `False` and the function should launch
        an exception (via the `launch_exception` argument).

        The exception class is the one of the `exception_type`
        argument.
    """
    low = operator.lt if low_exclusive else operator.le
    high = operator.lt if high_exclusive else operator.le
    wrapper = __get_wrapper(launch_exception=launch_exception,
                            exception_type=exception_type)
    try:
        return wrapper(((lower_bound is None) or
                        low(lower_bound, elem)) and
                       ((higher_bound is None) or
                        high(elem, higher_bound)), msg=msg)
    except TypeError as e:
        return wrapper(False, msg=e)

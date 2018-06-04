from __future__ import division, print_function, absolute_import

from distutils.version import LooseVersion
from functools import wraps

import numpy as np
from numpy.compat import basestring

from .core import (atop, asarray)
from . import chunk

einsum_symbols = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
einsum_symbols_set = set(einsum_symbols)


# This function duplicates numpy's _parse_einsum_input() function
# See https://github.com/numpy/numpy/blob/master/LICENSE.txt
# or NUMPY_LICENSE.txt within this directory
def parse_einsum_input(operands):
    """
    A reproduction of numpy's _parse_einsum_input()
    which in itself is a reproduction of
    c side einsum parsing in python.

    Returns
    -------
    input_strings : str
        Parsed input strings
    output_string : str
        Parsed output string
    operands : list of array_like
        The operands to use in the numpy contraction
    Examples
    --------
    The operand list is simplified to reduce printing:
    >> a = np.random.rand(4, 4)
    >> b = np.random.rand(4, 4, 4)
    >> __parse_einsum_input(('...a,...a->...', a, b))
    ('za,xza', 'xz', [a, b])
    >> __parse_einsum_input((a, [Ellipsis, 0], b, [Ellipsis, 0]))
    ('za,xza', 'xz', [a, b])
    """

    if len(operands) == 0:
        raise ValueError("No input operands")

    if isinstance(operands[0], basestring):
        subscripts = operands[0].replace(" ", "")
        operands = [asarray(o) for o in operands[1:]]

        # Ensure all characters are valid
        for s in subscripts:
            if s in '.,->':
                continue
            if s not in einsum_symbols_set:
                raise ValueError("Character %s is not a valid symbol." % s)

    else:
        tmp_operands = list(operands)
        operand_list = []
        subscript_list = []
        for p in range(len(operands) // 2):
            operand_list.append(tmp_operands.pop(0))
            subscript_list.append(tmp_operands.pop(0))

        output_list = tmp_operands[-1] if len(tmp_operands) else None
        operands = [asarray(v) for v in operand_list]
        subscripts = ""
        last = len(subscript_list) - 1
        for num, sub in enumerate(subscript_list):
            for s in sub:
                if s is Ellipsis:
                    subscripts += "..."
                elif isinstance(s, int):
                    subscripts += einsum_symbols[s]
                else:
                    raise TypeError("For this input type lists must contain "
                                    "either int or Ellipsis")
            if num != last:
                subscripts += ","

        if output_list is not None:
            subscripts += "->"
            for s in output_list:
                if s is Ellipsis:
                    subscripts += "..."
                elif isinstance(s, int):
                    subscripts += einsum_symbols[s]
                else:
                    raise TypeError("For this input type lists must contain "
                                    "either int or Ellipsis")
    # Check for proper "->"
    if ("-" in subscripts) or (">" in subscripts):
        invalid = (subscripts.count("-") > 1) or (subscripts.count(">") > 1)
        if invalid or (subscripts.count("->") != 1):
            raise ValueError("Subscripts can only contain one '->'.")

    # Parse ellipses
    if "." in subscripts:
        used = subscripts.replace(".", "").replace(",", "").replace("->", "")
        unused = list(einsum_symbols_set - set(used))
        ellipse_inds = "".join(unused)
        longest = 0

        if "->" in subscripts:
            input_tmp, output_sub = subscripts.split("->")
            split_subscripts = input_tmp.split(",")
            out_sub = True
        else:
            split_subscripts = subscripts.split(',')
            out_sub = False

        for num, sub in enumerate(split_subscripts):
            if "." in sub:
                if (sub.count(".") != 3) or (sub.count("...") != 1):
                    raise ValueError("Invalid Ellipses.")

                # Take into account numerical values
                if operands[num].shape == ():
                    ellipse_count = 0
                else:
                    ellipse_count = max(operands[num].ndim, 1)
                    ellipse_count -= (len(sub) - 3)

                if ellipse_count > longest:
                    longest = ellipse_count

                if ellipse_count < 0:
                    raise ValueError("Ellipses lengths do not match.")
                elif ellipse_count == 0:
                    split_subscripts[num] = sub.replace('...', '')
                else:
                    rep_inds = ellipse_inds[-ellipse_count:]
                    split_subscripts[num] = sub.replace('...', rep_inds)

        subscripts = ",".join(split_subscripts)
        if longest == 0:
            out_ellipse = ""
        else:
            out_ellipse = ellipse_inds[-longest:]

        if out_sub:
            subscripts += "->" + output_sub.replace("...", out_ellipse)
        else:
            # Special care for outputless ellipses
            output_subscript = ""
            tmp_subscripts = subscripts.replace(",", "")
            for s in sorted(set(tmp_subscripts)):
                if s not in einsum_symbols_set:
                    raise ValueError("Character %s is not a valid symbol." % s)
                if tmp_subscripts.count(s) == 1:
                    output_subscript += s
            normal_inds = ''.join(sorted(set(output_subscript) -
                                         set(out_ellipse)))

            subscripts += "->" + out_ellipse + normal_inds

    # Build output string if does not exist
    if "->" in subscripts:
        input_subscripts, output_subscript = subscripts.split("->")
    else:
        input_subscripts = subscripts
        # Build output subscripts
        tmp_subscripts = subscripts.replace(",", "")
        output_subscript = ""
        for s in sorted(set(tmp_subscripts)):
            if s not in einsum_symbols_set:
                raise ValueError("Character %s is not a valid symbol." % s)
            if tmp_subscripts.count(s) == 1:
                output_subscript += s

    # Make sure output subscripts are in the input
    for char in output_subscript:
        if char not in input_subscripts:
            raise ValueError("Output character %s did not appear in the input"
                             % char)

    # Make sure number operands is equivalent to the number of terms
    if len(input_subscripts.split(',')) != len(operands):
        raise ValueError("Number of einsum subscripts must be equal to the "
                         "number of operands.")

    return (input_subscripts, output_subscript, operands)


einsum_can_optimize = LooseVersion(np.__version__) >= LooseVersion("1.12.0")


@wraps(np.einsum)
def einsum(*operands, **kwargs):
    casting = kwargs.pop('casting', 'safe')
    dtype = kwargs.pop('dtype', None)
    optimize = kwargs.pop('optimize', False)
    order = kwargs.pop('order', 'K')
    split_every = kwargs.pop('split_every', None)
    if kwargs:
        raise TypeError("einsum() got unexpected keyword "
                        "argument(s) %s" % ",".join(kwargs))

    einsum_dtype = dtype

    inputs, outputs, ops = parse_einsum_input(operands)
    subscripts = '->'.join((inputs, outputs))

    # Infer the output dtype from operands
    if dtype is None:
        dtype = np.result_type(*[o.dtype for o in ops])

    if einsum_can_optimize:
        if optimize is not False:
            # Avoid computation of dask arrays within np.einsum_path
            # by passing in small numpy arrays broadcasted
            # up to the right shape
            fake_ops = [np.broadcast_to(o.dtype.type(0), shape=o.shape)
                        for o in ops]
            optimize, _ = np.einsum_path(subscripts, *fake_ops,
                                         optimize=optimize)
        kwargs = {'optimize': optimize}
    else:
        kwargs = {}

    inputs = [tuple(i) for i in inputs.split(",")]

    # Set of all indices
    all_inds = set(a for i in inputs for a in i)

    # Which indices are contracted?
    contract_inds = all_inds - set(outputs)
    ncontract_inds = len(contract_inds)

    # Introduce the contracted indices into the atop product
    # so that we get numpy arrays, not lists
    result = atop(chunk.einsum, tuple(outputs) + tuple(contract_inds),
                  *(a for ap in zip(ops, inputs) for a in ap),
                  # atop parameters
                  adjust_chunks={ind: 1 for ind in contract_inds}, dtype=dtype,
                  # np.einsum parameters
                  subscripts=subscripts, kernel_dtype=einsum_dtype,
                  ncontract_inds=ncontract_inds, order=order,
                  casting=casting, **kwargs)

    # Now reduce over any extra contraction dimensions
    if ncontract_inds > 0:
        size = len(outputs)
        return result.sum(axis=list(range(size, size + ncontract_inds)),
                          split_every=split_every)

    return result

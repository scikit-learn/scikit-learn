##########
Some plots
##########

Plot 1 does not use context:

.. plot::

    plt.plot(range(10))
    a = 10

Plot 2 doesn't use context either; has length 6:

.. plot::

    plt.plot(range(6))

Plot 3 has length 4:

.. plot::

    plt.plot(range(4))

Plot 4 shows that a new block with context does not see the variable defined
in the no-context block:

.. plot::
    :context:

    assert 'a' not in globals()

Plot 5 defines ``a`` in a context block:

.. plot::
    :context:

    plt.plot(range(6))
    a = 10

Plot 6 shows that a block with context sees the new variable.  It also uses
``:nofigs:``:

.. plot::
    :context:
    :nofigs:

    assert a == 10
    b = 4

Plot 7 uses a variable previously defined in previous ``nofigs`` context. It
also closes any previous figures to create a fresh figure:

.. plot::
    :context: close-figs

    assert b == 4
    plt.plot(range(b))

Plot 8 shows that a non-context block still doesn't have ``a``:

.. plot::
    :nofigs:

    assert 'a' not in globals()

Plot 9 has a context block, and does have ``a``:

.. plot::
    :context:
    :nofigs:

    assert a == 10

Plot 10 resets context, and ``a`` has gone again:

.. plot::
    :context: reset
    :nofigs:

    assert 'a' not in globals()
    c = 10

Plot 11 continues the context, we have the new value, but not the old:

.. plot::
    :context:

    assert c == 10
    assert 'a' not in globals()
    plt.plot(range(c))

Plot 12 opens a new figure.  By default the directive will plot both the first
and the second figure:

.. plot::
    :context:

    plt.figure()
    plt.plot(range(6))

Plot 13 shows ``close-figs`` in action.  ``close-figs`` closes all figures
previous to this plot directive, so we get always plot the figure we create in
the directive:

.. plot::
    :context: close-figs

    plt.figure()
    plt.plot(range(4))

Plot 14 uses ``include-source``:

.. plot::
    :include-source:

    # Only a comment

Plot 15 uses an external file with the plot commands and a caption:

.. plot:: range4.py

   This is the caption for plot 15.


Plot 16 uses a specific function in a file with plot commands:

.. plot:: range6.py range6



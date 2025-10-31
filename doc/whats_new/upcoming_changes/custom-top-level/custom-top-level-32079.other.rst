Free-threaded CPython 3.14 support
----------------------------------

scikit-learn has support for free-threaded CPython, in particular
free-threaded wheels are available for all of our supported platforms on Python
3.14.

Free-threaded (also known as nogil) CPython is a version of CPython that aims at
enabling efficient multi-threaded use cases by removing the Global Interpreter
Lock (GIL).

If you want to try out free-threaded Python, the recommendation is to use
Python 3.14, that has fixed a number of issues compared to Python 3.13. Feel
free to try free-threaded on your use case and report any issues!

For more details about free-threaded CPython see `py-free-threading doc <https://py-free-threading.github.io>`_,
in particular `how to install a free-threaded CPython <https://py-free-threading.github.io/installing_cpython/>`_
and `Ecosystem compatibility tracking <https://py-free-threading.github.io/tracking/>`_.

By :user:`Loïc Estève <lesteve>` and :user:`Olivier Grisel <ogrisel>` and many
other people in the wider Scientific Python and CPython ecosystem, for example
:user:`Nathan Goldbaum <ngoldbaum>`, :user:`Ralf Gommers <rgommers>`,
:user:`Edgar Andrés Margffoy Tuay <andfoy>`.

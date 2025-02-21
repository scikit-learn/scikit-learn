# cython: auto_pickle=False

r"""
Implements a buffer with insertion points. When you know you need to
"get back" to a place and write more later, simply call insertion_point()
at that spot and get a new StringIOTree object that is "left behind".

EXAMPLE:

>>> a = StringIOTree()
>>> _= a.write('first\n')
>>> b = a.insertion_point()
>>> _= a.write('third\n')
>>> _= b.write('second\n')
>>> a.getvalue().split()
['first', 'second', 'third']

>>> c = b.insertion_point()
>>> d = c.insertion_point()
>>> _= d.write('alpha\n')
>>> _= b.write('gamma\n')
>>> _= c.write('beta\n')
>>> b.getvalue().split()
['second', 'alpha', 'beta', 'gamma']

>>> try: from cStringIO import StringIO
... except ImportError: from io import StringIO

>>> i = StringIOTree()
>>> d.insert(i)
>>> _= i.write('inserted\n')
>>> out = StringIO()
>>> a.copyto(out)
>>> out.getvalue().split()
['first', 'second', 'alpha', 'inserted', 'beta', 'gamma', 'third']
"""

from __future__ import absolute_import  #, unicode_literals

try:
    # Prefer cStringIO since io.StringIO() does not support writing 'str' in Py2.
    from cStringIO import StringIO
except ImportError:
    from io import StringIO


class StringIOTree(object):
    """
    See module docs.
    """

    def __init__(self, stream=None):
        self.prepended_children = []
        if stream is None:
            stream = StringIO()
        self.stream = stream
        self.write = stream.write
        self.markers = []

    def empty(self):
        if self.stream.tell():
            return False
        return all([child.empty() for child in self.prepended_children]) if self.prepended_children else True

    def getvalue(self):
        content = []
        self._collect_in(content)
        return "".join(content)

    def _collect_in(self, target_list):
        for x in self.prepended_children:
            x._collect_in(target_list)
        stream_content = self.stream.getvalue()
        if stream_content:
            target_list.append(stream_content)

    def copyto(self, target):
        """Potentially cheaper than getvalue as no string concatenation
        needs to happen."""
        for child in self.prepended_children:
            child.copyto(target)
        stream_content = self.stream.getvalue()
        if stream_content:
            target.write(stream_content)

    def commit(self):
        # Save what we have written until now so that the buffer
        # itself is empty -- this makes it ready for insertion
        if self.stream.tell():
            self.prepended_children.append(StringIOTree(self.stream))
            self.prepended_children[-1].markers = self.markers
            self.markers = []
            self.stream = StringIO()
            self.write = self.stream.write

    def reset(self):
        self.prepended_children = []
        self.markers = []
        self.stream = StringIO()
        self.write = self.stream.write

    def insert(self, iotree):
        """
        Insert a StringIOTree (and all of its contents) at this location.
        Further writing to self appears after what is inserted.
        """
        self.commit()
        self.prepended_children.append(iotree)

    def insertion_point(self):
        """
        Returns a new StringIOTree, which is left behind at the current position
        (it what is written to the result will appear right before whatever is
        next written to self).

        Calling getvalue() or copyto() on the result will only return the
        contents written to it.
        """
        # Save what we have written until now
        # This is so that getvalue on the result doesn't include it.
        self.commit()
        # Construct the new forked object to return
        other = StringIOTree()
        self.prepended_children.append(other)
        return other

    def allmarkers(self):
        children = self.prepended_children
        return [m for c in children for m in c.allmarkers()] + self.markers

    """
    # Print the result of allmarkers in a nice human-readable form. Use it only for debugging.
    # Prints e.g.
    # /path/to/source.pyx:
    #     cython line 2 maps to 3299-3343
    #     cython line 4 maps to 2236-2245  2306  3188-3201
    # /path/to/othersource.pyx:
    #     cython line 3 maps to 1234-1270
    # ...
    # Note: In the example above, 3343 maps to line 2, 3344 does not.
    def print_hr_allmarkers(self):
        from collections import defaultdict
        markers = self.allmarkers()
        totmap = defaultdict(lambda: defaultdict(list))
        for c_lineno, (cython_desc, cython_lineno) in enumerate(markers):
            if cython_lineno > 0 and cython_desc.filename is not None:
                totmap[cython_desc.filename][cython_lineno].append(c_lineno + 1)
        reprstr = ""
        if totmap == 0:
            reprstr += "allmarkers is empty\n"
        try:
            sorted(totmap.items())
        except:
            print(totmap)
            print(totmap.items())
        for cython_path, filemap in sorted(totmap.items()):
            reprstr += cython_path + ":\n"
            for cython_lineno, c_linenos in sorted(filemap.items()):
                reprstr += "\tcython line " + str(cython_lineno) + " maps to "
                i = 0
                while i < len(c_linenos):
                    reprstr += str(c_linenos[i])
                    flag = False
                    while i+1 < len(c_linenos) and c_linenos[i+1] == c_linenos[i]+1:
                        i += 1
                        flag = True
                    if flag:
                        reprstr += "-" + str(c_linenos[i]) + " "
                    i += 1
                reprstr += "\n"

        import sys
        sys.stdout.write(reprstr)
    """

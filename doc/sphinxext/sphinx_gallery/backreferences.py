# -*- coding: utf-8 -*-
# Author: Óscar Nájera
# License: 3-clause BSD
"""
========================
Backreferences Generator
========================

Reviews generated example files in order to keep track of used modules
"""

from __future__ import print_function
import ast
import os


# Try Python 2 first, otherwise load from Python 3
try:
    import cPickle as pickle
except ImportError:
    import pickle


class NameFinder(ast.NodeVisitor):
    """Finds the longest form of variable names and their imports in code

    Only retains names from imported modules.
    """

    def __init__(self):
        super(NameFinder, self).__init__()
        self.imported_names = {}
        self.accessed_names = set()

    def visit_Import(self, node, prefix=''):
        for alias in node.names:
            local_name = alias.asname or alias.name
            self.imported_names[local_name] = prefix + alias.name

    def visit_ImportFrom(self, node):
        self.visit_Import(node, node.module + '.')

    def visit_Name(self, node):
        self.accessed_names.add(node.id)

    def visit_Attribute(self, node):
        attrs = []
        while isinstance(node, ast.Attribute):
            attrs.append(node.attr)
            node = node.value

        if isinstance(node, ast.Name):
            # This is a.b, not e.g. a().b
            attrs.append(node.id)
            self.accessed_names.add('.'.join(reversed(attrs)))
        else:
            # need to get a in a().b
            self.visit(node)

    def get_mapping(self):
        for name in self.accessed_names:
            local_name = name.split('.', 1)[0]
            remainder = name[len(local_name):]
            if local_name in self.imported_names:
                # Join import path to relative path
                full_name = self.imported_names[local_name] + remainder
                yield name, full_name


def get_short_module_name(module_name, obj_name):
    """ Get the shortest possible module name """
    parts = module_name.split('.')
    short_name = module_name
    for i in range(len(parts) - 1, 0, -1):
        short_name = '.'.join(parts[:i])
        try:
            exec('from %s import %s' % (short_name, obj_name))
        except ImportError:
            # get the last working module name
            short_name = '.'.join(parts[:(i + 1)])
            break
    return short_name


def identify_names(code):
    """Builds a codeobj summary by identifying and resolving used names

    >>> code = '''
    ... from a.b import c
    ... import d as e
    ... print(c)
    ... e.HelloWorld().f.g
    ... '''
    >>> for name, o in sorted(identify_names(code).items()):
    ...     print(name, o['name'], o['module'], o['module_short'])
    c c a.b a.b
    e.HelloWorld HelloWorld d d
    """
    finder = NameFinder()
    finder.visit(ast.parse(code))

    example_code_obj = {}
    for name, full_name in finder.get_mapping():
        # name is as written in file (e.g. np.asarray)
        # full_name includes resolved import path (e.g. numpy.asarray)
        module, attribute = full_name.rsplit('.', 1)
        # get shortened module name
        module_short = get_short_module_name(module, attribute)
        cobj = {'name': attribute, 'module': module,
                'module_short': module_short}
        example_code_obj[name] = cobj
    return example_code_obj


def scan_used_functions(example_file, gallery_conf):
    """save variables so we can later add links to the documentation"""
    example_code_obj = identify_names(open(example_file).read())
    if example_code_obj:
        codeobj_fname = example_file[:-3] + '_codeobj.pickle'
        with open(codeobj_fname, 'wb') as fid:
            pickle.dump(example_code_obj, fid, pickle.HIGHEST_PROTOCOL)

    backrefs = set('{module_short}.{name}'.format(**entry)
                   for entry in example_code_obj.values()
                   if entry['module'].startswith(gallery_conf['doc_module']))

    return backrefs


THUMBNAIL_TEMPLATE = """
.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="{snippet}">

.. only:: html

    .. figure:: /{thumbnail}

        :ref:`sphx_glr_{ref_name}`

.. raw:: html

    </div>
"""

BACKREF_THUMBNAIL_TEMPLATE = THUMBNAIL_TEMPLATE + """
.. only:: not html

    * :ref:`sphx_glr_{ref_name}`
"""


def _thumbnail_div(full_dir, fname, snippet, is_backref=False):
    """Generates RST to place a thumbnail in a gallery"""
    thumb = os.path.join(full_dir, 'images', 'thumb',
                         'sphx_glr_%s_thumb.png' % fname[:-3])
    ref_name = os.path.join(full_dir, fname).replace(os.path.sep, '_')

    template = BACKREF_THUMBNAIL_TEMPLATE if is_backref else THUMBNAIL_TEMPLATE
    return template.format(snippet=snippet, thumbnail=thumb, ref_name=ref_name)


def write_backreferences(seen_backrefs, gallery_conf,
                         target_dir, fname, snippet):
    """Writes down back reference files, which include a thumbnail list
    of examples using a certain module"""
    example_file = os.path.join(target_dir, fname)
    backrefs = scan_used_functions(example_file, gallery_conf)
    for backref in backrefs:
        include_path = os.path.join(gallery_conf['mod_example_dir'],
                                    '%s.examples' % backref)
        seen = backref in seen_backrefs
        with open(include_path, 'a' if seen else 'w') as ex_file:
            if not seen:
                heading = '\n\nExamples using ``%s``' % backref
                ex_file.write(heading + '\n')
                ex_file.write('^' * len(heading) + '\n')
            ex_file.write(_thumbnail_div(target_dir, fname, snippet,
                                         is_backref=True))
            seen_backrefs.add(backref)

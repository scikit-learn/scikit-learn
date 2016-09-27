"""
==============
phantom_import
==============

Sphinx extension to make directives from ``sphinx.ext.autodoc`` and similar
extensions to use docstrings loaded from an XML file.

This extension loads an XML file in the Pydocweb format [1] and
creates a dummy module that contains the specified docstrings. This
can be used to get the current docstrings from a Pydocweb instance
without needing to rebuild the documented module.

.. [1] http://code.google.com/p/pydocweb

"""
from __future__ import division, absolute_import, print_function

import imp, sys, compiler, types, os, inspect, re

def setup(app):
    app.connect('builder-inited', initialize)
    app.add_config_value('phantom_import_file', None, True)

def initialize(app):
    fn = app.config.phantom_import_file
    if (fn and os.path.isfile(fn)):
        print("[numpydoc] Phantom importing modules from", fn, "...")
        import_phantom_module(fn)

#------------------------------------------------------------------------------
# Creating 'phantom' modules from an XML description
#------------------------------------------------------------------------------
def import_phantom_module(xml_file):
    """
    Insert a fake Python module to sys.modules, based on a XML file.

    The XML file is expected to conform to Pydocweb DTD. The fake
    module will contain dummy objects, which guarantee the following:

    - Docstrings are correct.
    - Class inheritance relationships are correct (if present in XML).
    - Function argspec is *NOT* correct (even if present in XML).
      Instead, the function signature is prepended to the function docstring.
    - Class attributes are *NOT* correct; instead, they are dummy objects.

    Parameters
    ----------
    xml_file : str
        Name of an XML file to read
    
    """
    import lxml.etree as etree

    object_cache = {}

    tree = etree.parse(xml_file)
    root = tree.getroot()

    # Sort items so that
    # - Base classes come before classes inherited from them
    # - Modules come before their contents
    all_nodes = dict([(n.attrib['id'], n) for n in root])
    
    def _get_bases(node, recurse=False):
        bases = [x.attrib['ref'] for x in node.findall('base')]
        if recurse:
            j = 0
            while True:
                try:
                    b = bases[j]
                except IndexError: break
                if b in all_nodes:
                    bases.extend(_get_bases(all_nodes[b]))
                j += 1
        return bases

    type_index = ['module', 'class', 'callable', 'object']
    
    def base_cmp(a, b):
        x = cmp(type_index.index(a.tag), type_index.index(b.tag))
        if x != 0: return x

        if a.tag == 'class' and b.tag == 'class':
            a_bases = _get_bases(a, recurse=True)
            b_bases = _get_bases(b, recurse=True)
            x = cmp(len(a_bases), len(b_bases))
            if x != 0: return x
            if a.attrib['id'] in b_bases: return -1
            if b.attrib['id'] in a_bases: return 1
        
        return cmp(a.attrib['id'].count('.'), b.attrib['id'].count('.'))

    nodes = root.getchildren()
    nodes.sort(base_cmp)

    # Create phantom items
    for node in nodes:
        name = node.attrib['id']
        doc = (node.text or '').decode('string-escape') + "\n"
        if doc == "\n": doc = ""

        # create parent, if missing
        parent = name
        while True:
            parent = '.'.join(parent.split('.')[:-1])
            if not parent: break
            if parent in object_cache: break
            obj = imp.new_module(parent)
            object_cache[parent] = obj
            sys.modules[parent] = obj

        # create object
        if node.tag == 'module':
            obj = imp.new_module(name)
            obj.__doc__ = doc
            sys.modules[name] = obj
        elif node.tag == 'class':
            bases = [object_cache[b] for b in _get_bases(node)
                     if b in object_cache]
            bases.append(object)
            init = lambda self: None
            init.__doc__ = doc
            obj = type(name, tuple(bases), {'__doc__': doc, '__init__': init})
            obj.__name__ = name.split('.')[-1]
        elif node.tag == 'callable':
            funcname = node.attrib['id'].split('.')[-1]
            argspec = node.attrib.get('argspec')
            if argspec:
                argspec = re.sub('^[^(]*', '', argspec)
                doc = "%s%s\n\n%s" % (funcname, argspec, doc)
            obj = lambda: 0
            obj.__argspec_is_invalid_ = True
            if sys.version_info[0] >= 3:
                obj.__name__ = funcname
            else:
                obj.func_name = funcname
            obj.__name__ = name
            obj.__doc__ = doc
            if inspect.isclass(object_cache[parent]):
                obj.__objclass__ = object_cache[parent]
        else:
            class Dummy(object): pass
            obj = Dummy()
            obj.__name__ = name
            obj.__doc__ = doc
            if inspect.isclass(object_cache[parent]):
                obj.__get__ = lambda: None
        object_cache[name] = obj

        if parent:
            if inspect.ismodule(object_cache[parent]):
                obj.__module__ = parent
                setattr(object_cache[parent], name.split('.')[-1], obj)

    # Populate items
    for node in root:
        obj = object_cache.get(node.attrib['id'])
        if obj is None: continue
        for ref in node.findall('ref'):
            if node.tag == 'class':
                if ref.attrib['ref'].startswith(node.attrib['id'] + '.'):
                    setattr(obj, ref.attrib['name'],
                            object_cache.get(ref.attrib['ref']))
            else:
                setattr(obj, ref.attrib['name'],
                        object_cache.get(ref.attrib['ref']))

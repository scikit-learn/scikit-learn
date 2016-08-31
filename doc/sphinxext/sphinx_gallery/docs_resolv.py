# -*- coding: utf-8 -*-
# Author: Óscar Nájera
# License: 3-clause BSD
###############################################################################
# Documentation link resolver objects
from __future__ import print_function
import gzip
import os
import posixpath
import re
import shelve
import sys

# Try Python 2 first, otherwise load from Python 3
try:
    from StringIO import StringIO
    import cPickle as pickle
    import urllib2 as urllib
    from urllib2 import HTTPError, URLError
except ImportError:
    from io import StringIO
    import pickle
    import urllib.request
    import urllib.error
    import urllib.parse
    from urllib.error import HTTPError, URLError


def _get_data(url):
    """Helper function to get data over http or from a local file"""
    if url.startswith('http://'):
        # Try Python 2, use Python 3 on exception
        try:
            resp = urllib.urlopen(url)
            encoding = resp.headers.dict.get('content-encoding', 'plain')
        except AttributeError:
            resp = urllib.request.urlopen(url)
            encoding = resp.headers.get('content-encoding', 'plain')
        data = resp.read()
        if encoding == 'plain':
            pass
        elif encoding == 'gzip':
            data = StringIO(data)
            data = gzip.GzipFile(fileobj=data).read()
        else:
            raise RuntimeError('unknown encoding')
    else:
        with open(url, 'r') as fid:
            data = fid.read()

    return data


def get_data(url, gallery_dir):
    """Persistent dictionary usage to retrieve the search indexes"""

    # shelve keys need to be str in python 2
    if sys.version_info[0] == 2 and isinstance(url, unicode):
        url = url.encode('utf-8')

    cached_file = os.path.join(gallery_dir, 'searchindex')
    search_index = shelve.open(cached_file)
    if url in search_index:
        data = search_index[url]
    else:
        data = _get_data(url)
        search_index[url] = data
    search_index.close()

    return data


def _select_block(str_in, start_tag, end_tag):
    """Select first block delimited by start_tag and end_tag"""
    start_pos = str_in.find(start_tag)
    if start_pos < 0:
        raise ValueError('start_tag not found')
    depth = 0
    for pos in range(start_pos, len(str_in)):
        if str_in[pos] == start_tag:
            depth += 1
        elif str_in[pos] == end_tag:
            depth -= 1

        if depth == 0:
            break
    sel = str_in[start_pos + 1:pos]
    return sel


def _parse_dict_recursive(dict_str):
    """Parse a dictionary from the search index"""
    dict_out = dict()
    pos_last = 0
    pos = dict_str.find(':')
    while pos >= 0:
        key = dict_str[pos_last:pos]
        if dict_str[pos + 1] == '[':
            # value is a list
            pos_tmp = dict_str.find(']', pos + 1)
            if pos_tmp < 0:
                raise RuntimeError('error when parsing dict')
            value = dict_str[pos + 2: pos_tmp].split(',')
            # try to convert elements to int
            for i in range(len(value)):
                try:
                    value[i] = int(value[i])
                except ValueError:
                    pass
        elif dict_str[pos + 1] == '{':
            # value is another dictionary
            subdict_str = _select_block(dict_str[pos:], '{', '}')
            value = _parse_dict_recursive(subdict_str)
            pos_tmp = pos + len(subdict_str)
        else:
            raise ValueError('error when parsing dict: unknown elem')

        key = key.strip('"')
        if len(key) > 0:
            dict_out[key] = value

        pos_last = dict_str.find(',', pos_tmp)
        if pos_last < 0:
            break
        pos_last += 1
        pos = dict_str.find(':', pos_last)

    return dict_out


def parse_sphinx_searchindex(searchindex):
    """Parse a Sphinx search index

    Parameters
    ----------
    searchindex : str
        The Sphinx search index (contents of searchindex.js)

    Returns
    -------
    filenames : list of str
        The file names parsed from the search index.
    objects : dict
        The objects parsed from the search index.
    """
    # Make sure searchindex uses UTF-8 encoding
    if hasattr(searchindex, 'decode'):
        searchindex = searchindex.decode('UTF-8')

    # parse objects
    query = 'objects:'
    pos = searchindex.find(query)
    if pos < 0:
        raise ValueError('"objects:" not found in search index')

    sel = _select_block(searchindex[pos:], '{', '}')
    objects = _parse_dict_recursive(sel)

    # parse filenames
    query = 'filenames:'
    pos = searchindex.find(query)
    if pos < 0:
        raise ValueError('"filenames:" not found in search index')
    filenames = searchindex[pos + len(query) + 1:]
    filenames = filenames[:filenames.find(']')]
    filenames = [f.strip('"') for f in filenames.split(',')]

    return filenames, objects


class SphinxDocLinkResolver(object):
    """ Resolve documentation links using searchindex.js generated by Sphinx

    Parameters
    ----------
    doc_url : str
        The base URL of the project website.
    searchindex : str
        Filename of searchindex, relative to doc_url.
    extra_modules_test : list of str
        List of extra module names to test.
    relative : bool
        Return relative links (only useful for links to documentation of this
        package).
    """

    def __init__(self, doc_url, gallery_dir, searchindex='searchindex.js',
                 extra_modules_test=None, relative=False):
        self.doc_url = doc_url
        self.gallery_dir = gallery_dir
        self.relative = relative
        self._link_cache = {}

        self.extra_modules_test = extra_modules_test
        self._page_cache = {}
        if doc_url.startswith('http://'):
            if relative:
                raise ValueError('Relative links are only supported for local '
                                 'URLs (doc_url cannot start with "http://)"')
            searchindex_url = doc_url + '/' + searchindex
        else:
            searchindex_url = os.path.join(doc_url, searchindex)

        # detect if we are using relative links on a Windows system
        if os.name.lower() == 'nt' and not doc_url.startswith('http://'):
            if not relative:
                raise ValueError('You have to use relative=True for the local'
                                 ' package on a Windows system.')
            self._is_windows = True
        else:
            self._is_windows = False

        # download and initialize the search index
        sindex = get_data(searchindex_url, gallery_dir)
        filenames, objects = parse_sphinx_searchindex(sindex)

        self._searchindex = dict(filenames=filenames, objects=objects)

    def _get_link(self, cobj):
        """Get a valid link, False if not found"""

        fname_idx = None
        full_name = cobj['module_short'] + '.' + cobj['name']
        if full_name in self._searchindex['objects']:
            value = self._searchindex['objects'][full_name]
            if isinstance(value, dict):
                value = value[next(iter(value.keys()))]
            fname_idx = value[0]
        elif cobj['module_short'] in self._searchindex['objects']:
            value = self._searchindex['objects'][cobj['module_short']]
            if cobj['name'] in value.keys():
                fname_idx = value[cobj['name']][0]

        if fname_idx is not None:
            fname = self._searchindex['filenames'][fname_idx] + '.html'

            if self._is_windows:
                fname = fname.replace('/', '\\')
                link = os.path.join(self.doc_url, fname)
            else:
                link = posixpath.join(self.doc_url, fname)

            if hasattr(link, 'decode'):
                link = link.decode('utf-8', 'replace')

            if link in self._page_cache:
                html = self._page_cache[link]
            else:
                html = get_data(link, self.gallery_dir)
                self._page_cache[link] = html

            # test if cobj appears in page
            comb_names = [cobj['module_short'] + '.' + cobj['name']]
            if self.extra_modules_test is not None:
                for mod in self.extra_modules_test:
                    comb_names.append(mod + '.' + cobj['name'])
            url = False
            if hasattr(html, 'decode'):
                # Decode bytes under Python 3
                html = html.decode('utf-8', 'replace')

            for comb_name in comb_names:
                if hasattr(comb_name, 'decode'):
                    # Decode bytes under Python 3
                    comb_name = comb_name.decode('utf-8', 'replace')
                if comb_name in html:
                    url = link + u'#' + comb_name
            link = url
        else:
            link = False

        return link

    def resolve(self, cobj, this_url):
        """Resolve the link to the documentation, returns None if not found

        Parameters
        ----------
        cobj : dict
            Dict with information about the "code object" for which we are
            resolving a link.
            cobi['name'] : function or class name (str)
            cobj['module_short'] : shortened module name (str)
            cobj['module'] : module name (str)
        this_url: str
            URL of the current page. Needed to construct relative URLs
            (only used if relative=True in constructor).

        Returns
        -------
        link : str | None
            The link (URL) to the documentation.
        """
        full_name = cobj['module_short'] + '.' + cobj['name']
        link = self._link_cache.get(full_name, None)
        if link is None:
            # we don't have it cached
            link = self._get_link(cobj)
            # cache it for the future
            self._link_cache[full_name] = link

        if link is False or link is None:
            # failed to resolve
            return None

        if self.relative:
            link = os.path.relpath(link, start=this_url)
            if self._is_windows:
                # replace '\' with '/' so it on the web
                link = link.replace('\\', '/')

            # for some reason, the relative link goes one directory too high up
            link = link[3:]

        return link


def _embed_code_links(app, gallery_conf, gallery_dir):
    # Add resolvers for the packages for which we want to show links
    doc_resolvers = {}

    for this_module, url in gallery_conf['reference_url'].items():
        try:
            if url is None:
                doc_resolvers[this_module] = SphinxDocLinkResolver(
                    app.builder.outdir,
                    gallery_dir,
                    relative=True)
            else:
                doc_resolvers[this_module] = SphinxDocLinkResolver(url,
                                                                   gallery_dir)

        except HTTPError as e:
            print("The following HTTP Error has occurred:\n")
            print(e.code)
        except URLError as e:
            print("\n...\n"
                  "Warning: Embedding the documentation hyperlinks requires "
                  "Internet access.\nPlease check your network connection.\n"
                  "Unable to continue embedding `{0}` links due to a URL "
                  "Error:\n".format(this_module))
            print(e.args)

    html_gallery_dir = os.path.abspath(os.path.join(app.builder.outdir,
                                                    gallery_dir))

    # patterns for replacement
    link_pattern = '<a href="%s">%s</a>'
    orig_pattern = '<span class="n">%s</span>'
    period = '<span class="o">.</span>'

    for dirpath, _, filenames in os.walk(html_gallery_dir):
        for fname in filenames:
            print('\tprocessing: %s' % fname)
            full_fname = os.path.join(html_gallery_dir, dirpath, fname)
            subpath = dirpath[len(html_gallery_dir) + 1:]
            pickle_fname = os.path.join(gallery_dir, subpath,
                                        fname[:-5] + '_codeobj.pickle')

            if os.path.exists(pickle_fname):
                # we have a pickle file with the objects to embed links for
                with open(pickle_fname, 'rb') as fid:
                    example_code_obj = pickle.load(fid)
                fid.close()
                str_repl = {}
                # generate replacement strings with the links
                for name, cobj in example_code_obj.items():
                    this_module = cobj['module'].split('.')[0]

                    if this_module not in doc_resolvers:
                        continue

                    try:
                        link = doc_resolvers[this_module].resolve(cobj,
                                                                  full_fname)
                    except (HTTPError, URLError) as e:
                        print("The following error has occurred:\n")
                        print(repr(e))
                        continue

                    if link is not None:
                        parts = name.split('.')
                        name_html = period.join(orig_pattern % part
                                                for part in parts)
                        str_repl[name_html] = link_pattern % (link, name_html)
                # do the replacement in the html file

                # ensure greediness
                names = sorted(str_repl, key=len, reverse=True)
                expr = re.compile(r'(?<!\.)\b' +  # don't follow . or word
                                  '|'.join(re.escape(name)
                                           for name in names))

                def substitute_link(match):
                    return str_repl[match.group()]

                if len(str_repl) > 0:
                    with open(full_fname, 'rb') as fid:
                        lines_in = fid.readlines()
                    with open(full_fname, 'wb') as fid:
                        for line in lines_in:
                            line = line.decode('utf-8')
                            line = expr.sub(substitute_link, line)
                            fid.write(line.encode('utf-8'))
    print('[done]')


def embed_code_links(app, exception):
    """Embed hyperlinks to documentation into example code"""
    if exception is not None:
        return

    # No need to waste time embedding hyperlinks when not running the examples
    # XXX: also at the time of writing this fixes make html-noplot
    # for some reason I don't fully understand
    if not app.builder.config.plot_gallery:
        return

    # XXX: Whitelist of builders for which it makes sense to embed
    # hyperlinks inside the example html. Note that the link embedding
    # require searchindex.js to exist for the links to the local doc
    # and there does not seem to be a good way of knowing which
    # builders creates a searchindex.js.
    if app.builder.name not in ['html', 'readthedocs']:
        return

    print('Embedding documentation hyperlinks in examples..')

    gallery_conf = app.config.sphinx_gallery_conf

    gallery_dirs = gallery_conf['gallery_dirs']
    if not isinstance(gallery_dirs, list):
        gallery_dirs = [gallery_dirs]

    for gallery_dir in gallery_dirs:
        _embed_code_links(app, gallery_conf, gallery_dir)

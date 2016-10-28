
# simple python script to collect text paragraphs from various languages on the
# same topic namely the Wikipedia encyclopedia itself

import os
try:
    # Python 2 compat
    from urllib2 import Request, build_opener
except ImportError:
    # Python 3
    from urllib.request import Request, build_opener

import lxml.html
from lxml.etree import ElementTree
import numpy as np

import codecs

pages = {
    u'ar': u'http://ar.wikipedia.org/wiki/%D9%88%D9%8A%D9%83%D9%8A%D8%A8%D9%8A%D8%AF%D9%8A%D8%A7',
    u'de': u'http://de.wikipedia.org/wiki/Wikipedia',
    u'en': u'https://en.wikipedia.org/wiki/Wikipedia',
    u'es': u'http://es.wikipedia.org/wiki/Wikipedia',
    u'fr': u'http://fr.wikipedia.org/wiki/Wikip%C3%A9dia',
    u'it': u'http://it.wikipedia.org/wiki/Wikipedia',
    u'ja': u'http://ja.wikipedia.org/wiki/Wikipedia',
    u'nl': u'http://nl.wikipedia.org/wiki/Wikipedia',
    u'pl': u'http://pl.wikipedia.org/wiki/Wikipedia',
    u'pt': u'http://pt.wikipedia.org/wiki/Wikip%C3%A9dia',
    u'ru': u'http://ru.wikipedia.org/wiki/%D0%92%D0%B8%D0%BA%D0%B8%D0%BF%D0%B5%D0%B4%D0%B8%D1%8F',
#    u'zh': u'http://zh.wikipedia.org/wiki/Wikipedia',
}

html_folder = u'html'
text_folder = u'paragraphs'
short_text_folder = u'short_paragraphs'
n_words_per_short_text = 5


if not os.path.exists(html_folder):
    os.makedirs(html_folder)

for lang, page in pages.items():

    text_lang_folder = os.path.join(text_folder, lang)
    if not os.path.exists(text_lang_folder):
        os.makedirs(text_lang_folder)

    short_text_lang_folder = os.path.join(short_text_folder, lang)
    if not os.path.exists(short_text_lang_folder):
        os.makedirs(short_text_lang_folder)

    opener = build_opener()
    html_filename = os.path.join(html_folder, lang + '.html')
    if not os.path.exists(html_filename):
        print("Downloading %s" % page)
        request = Request(page)
        # change the User Agent to avoid being blocked by Wikipedia
        # downloading a couple of articles should not be considered abusive
        request.add_header('User-Agent', 'OpenAnything/1.0')
        html_content = opener.open(request).read()
        open(html_filename, 'wb').write(html_content)

    # decode the payload explicitly as UTF-8 since lxml is confused for some
    # reason
    with codecs.open(html_filename,'r','utf-8') as html_file:
        html_content = html_file.read()
    tree = ElementTree(lxml.html.document_fromstring(html_content))
    i = 0
    j = 0
    for p in tree.findall('//p'):
        content = p.text_content()
        if len(content) < 100:
            # skip paragraphs that are too short - probably too noisy and not
            # representative of the actual language
            continue

        text_filename = os.path.join(text_lang_folder,
                                     '%s_%04d.txt' % (lang, i))
        print("Writing %s" % text_filename)
        open(text_filename, 'wb').write(content.encode('utf-8', 'ignore'))
        i += 1

        # split the paragraph into fake smaller paragraphs to make the
        # problem harder e.g. more similar to tweets
        if lang in ('zh', 'ja'):
        # FIXME: whitespace tokenizing does not work on chinese and japanese
            continue
        words = content.split()
        n_groups = len(words) / n_words_per_short_text
        if n_groups < 1:
            continue
        groups = np.array_split(words, n_groups)

        for group in groups:
            small_content = u" ".join(group)

            short_text_filename = os.path.join(short_text_lang_folder,
                                               '%s_%04d.txt' % (lang, j))
            print("Writing %s" % short_text_filename)
            open(short_text_filename, 'wb').write(
                small_content.encode('utf-8', 'ignore'))
            j += 1
            if j >= 1000:
                break


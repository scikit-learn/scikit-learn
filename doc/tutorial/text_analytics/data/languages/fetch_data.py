
# simple python script to collect text paragraphs from various languages on the
# same topic namely the Wikipedia encyclopedia itself

import os
from urllib.request import Request, build_opener

import lxml.html
from lxml.etree import ElementTree
import numpy as np

import codecs

pages = {
    'ar': 'http://ar.wikipedia.org/wiki/%D9%88%D9%8A%D9%83%D9%8A%D8%A8%D9%8A%D8%AF%D9%8A%D8%A7',   # noqa: E501
    'de': 'http://de.wikipedia.org/wiki/Wikipedia',
    'en': 'https://en.wikipedia.org/wiki/Wikipedia',
    'es': 'http://es.wikipedia.org/wiki/Wikipedia',
    'fr': 'http://fr.wikipedia.org/wiki/Wikip%C3%A9dia',
    'it': 'http://it.wikipedia.org/wiki/Wikipedia',
    'ja': 'http://ja.wikipedia.org/wiki/Wikipedia',
    'nl': 'http://nl.wikipedia.org/wiki/Wikipedia',
    'pl': 'http://pl.wikipedia.org/wiki/Wikipedia',
    'pt': 'http://pt.wikipedia.org/wiki/Wikip%C3%A9dia',
    'ru': 'http://ru.wikipedia.org/wiki/%D0%92%D0%B8%D0%BA%D0%B8%D0%BF%D0%B5%D0%B4%D0%B8%D1%8F',  # noqa: E501
#    u'zh': u'http://zh.wikipedia.org/wiki/Wikipedia',
}

html_folder = 'html'
text_folder = 'paragraphs'
short_text_folder = 'short_paragraphs'
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
        with open(html_filename, 'wb') as f:
            f.write(html_content)

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
        with open(text_filename, 'wb') as f:
            f.write(content.encode('utf-8', 'ignore'))
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
            small_content = " ".join(group)

            short_text_filename = os.path.join(short_text_lang_folder,
                                               '%s_%04d.txt' % (lang, j))
            print("Writing %s" % short_text_filename)
            with open(short_text_filename, 'wb') as f:
                f.write(small_content.encode('utf-8', 'ignore'))
            j += 1
            if j >= 1000:
                break


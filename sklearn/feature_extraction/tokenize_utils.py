"""
This module contains utility functions for breaking down
a given block of text into it's component syntactic parts.

Author: Max Mautner

https://github.com/mmautner/readability/blob/master/utils.py
"""

import nltk

from nltk.tokenize import RegexpTokenizer
from . import syllables_en

TOKENIZER = RegexpTokenizer('(?u)\W+|\$[\d\.]+|\S+')
SPECIAL_CHARS = ['.', ',', '!', '?']

def get_char_count(words):
    characters = 0
    for word in words:
        characters += len(word.decode("utf-8"))
    return characters

def get_words(text=''):
    words = []
    words = TOKENIZER.tokenize(text)
    filtered_words = []
    for word in words:
        if word in SPECIAL_CHARS or word == " ":
            pass
        else:
            new_word = word.replace(",","").replace(".","")
            new_word = new_word.replace("!","").replace("?","")
            filtered_words.append(new_word)
    return filtered_words

def get_sentences(text=''):
    tokenizer = nltk.data.load('tokenizers/punkt/english.pickle')
    sentences = tokenizer.tokenize(text)
    return sentences

def count_syllables(words):
    syllableCount = 0
    for word in words:
        syllableCount += syllables_en.count(word)
    return syllableCount

#This method must be enhanced. At the moment it only
#considers the number of syllables in a word.
#This often results in that too many complex words are detected.
def count_complex_words(text=''):
    words = get_words(text)
    sentences = get_sentences(text)
    complex_words = 0
    found = False
    cur_word = []

    for word in words:
        cur_word.append(word)
        if count_syllables(cur_word)>= 3:

            #Checking proper nouns. If a word starts with a capital letter
            #and is NOT at the beginning of a sentence we don't add it
            #as a complex word.
            if not(word[0].isupper()):
                complex_words += 1
            else:
                for sentence in sentences:
                    if str(sentence).startswith(word):
                        found = True
                        break
                if found:
                    complex_words += 1
                    found = False

        cur_word.remove(word)
    return complex_words


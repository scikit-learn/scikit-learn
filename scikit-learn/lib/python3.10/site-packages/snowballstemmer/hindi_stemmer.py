#-*- coding: utf-8 -*-
# Generated from hindi.sbl by Snowball 3.0.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class HindiStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from hindi.sbl by Snowball 3.0.1 - https://snowballstem.org/
    '''

    g_consonant = {u"\u0915", u"\u0916", u"\u0917", u"\u0918", u"\u0919", u"\u091A", u"\u091B", u"\u091C", u"\u091D", u"\u091E", u"\u091F", u"\u0920", u"\u0921", u"\u0922", u"\u0923", u"\u0924", u"\u0925", u"\u0926", u"\u0927", u"\u0928", u"\u092A", u"\u092B", u"\u092C", u"\u092D", u"\u092E", u"\u092F", u"\u0930", u"\u0932", u"\u0935", u"\u0936", u"\u0937", u"\u0938", u"\u0939", u"\u0933", u"\u093C", u"\u0929", u"\u0931", u"\u0934", u"\u0958", u"\u0959", u"\u095A", u"\u095B", u"\u095C", u"\u095D", u"\u095E", u"\u095F"}


    def __r_CONSONANT(self):
        if not self.in_grouping_b(HindiStemmer.g_consonant):
            return False
        return True

    def _stem(self):
        if self.cursor >= self.limit:
            return False
        self.cursor += 1
        self.limit_backward = self.cursor
        self.cursor = self.limit
        self.ket = self.cursor
        if self.find_among_b(HindiStemmer.a_0) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among(u"\u0906\u0901", -1, -1),
        Among(u"\u093E\u0901", -1, -1),
        Among(u"\u0907\u092F\u093E\u0901", 1, -1),
        Among(u"\u0906\u0907\u092F\u093E\u0901", 2, -1),
        Among(u"\u093E\u0907\u092F\u093E\u0901", 2, -1),
        Among(u"\u093F\u092F\u093E\u0901", 1, -1),
        Among(u"\u0906\u0902", -1, -1),
        Among(u"\u0909\u0906\u0902", 6, -1),
        Among(u"\u0941\u0906\u0902", 6, -1),
        Among(u"\u0908\u0902", -1, -1),
        Among(u"\u0906\u0908\u0902", 9, -1),
        Among(u"\u093E\u0908\u0902", 9, -1),
        Among(u"\u090F\u0902", -1, -1),
        Among(u"\u0906\u090F\u0902", 12, -1),
        Among(u"\u0909\u090F\u0902", 12, -1),
        Among(u"\u093E\u090F\u0902", 12, -1),
        Among(u"\u0924\u093E\u090F\u0902", 15, -1, __r_CONSONANT),
        Among(u"\u0905\u0924\u093E\u090F\u0902", 16, -1),
        Among(u"\u0928\u093E\u090F\u0902", 15, -1, __r_CONSONANT),
        Among(u"\u0905\u0928\u093E\u090F\u0902", 18, -1),
        Among(u"\u0941\u090F\u0902", 12, -1),
        Among(u"\u0913\u0902", -1, -1),
        Among(u"\u0906\u0913\u0902", 21, -1),
        Among(u"\u0909\u0913\u0902", 21, -1),
        Among(u"\u093E\u0913\u0902", 21, -1),
        Among(u"\u0924\u093E\u0913\u0902", 24, -1, __r_CONSONANT),
        Among(u"\u0905\u0924\u093E\u0913\u0902", 25, -1),
        Among(u"\u0928\u093E\u0913\u0902", 24, -1, __r_CONSONANT),
        Among(u"\u0905\u0928\u093E\u0913\u0902", 27, -1),
        Among(u"\u0941\u0913\u0902", 21, -1),
        Among(u"\u093E\u0902", -1, -1),
        Among(u"\u0907\u092F\u093E\u0902", 30, -1),
        Among(u"\u0906\u0907\u092F\u093E\u0902", 31, -1),
        Among(u"\u093E\u0907\u092F\u093E\u0902", 31, -1),
        Among(u"\u093F\u092F\u093E\u0902", 30, -1),
        Among(u"\u0940\u0902", -1, -1),
        Among(u"\u0924\u0940\u0902", 35, -1, __r_CONSONANT),
        Among(u"\u0905\u0924\u0940\u0902", 36, -1),
        Among(u"\u0906\u0924\u0940\u0902", 36, -1),
        Among(u"\u093E\u0924\u0940\u0902", 36, -1),
        Among(u"\u0947\u0902", -1, -1),
        Among(u"\u094B\u0902", -1, -1),
        Among(u"\u0907\u092F\u094B\u0902", 41, -1),
        Among(u"\u0906\u0907\u092F\u094B\u0902", 42, -1),
        Among(u"\u093E\u0907\u092F\u094B\u0902", 42, -1),
        Among(u"\u093F\u092F\u094B\u0902", 41, -1),
        Among(u"\u0905", -1, -1),
        Among(u"\u0906", -1, -1),
        Among(u"\u0907", -1, -1),
        Among(u"\u0908", -1, -1),
        Among(u"\u0906\u0908", 49, -1),
        Among(u"\u093E\u0908", 49, -1),
        Among(u"\u0909", -1, -1),
        Among(u"\u090A", -1, -1),
        Among(u"\u090F", -1, -1),
        Among(u"\u0906\u090F", 54, -1),
        Among(u"\u0907\u090F", 54, -1),
        Among(u"\u0906\u0907\u090F", 56, -1),
        Among(u"\u093E\u0907\u090F", 56, -1),
        Among(u"\u093E\u090F", 54, -1),
        Among(u"\u093F\u090F", 54, -1),
        Among(u"\u0913", -1, -1),
        Among(u"\u0906\u0913", 61, -1),
        Among(u"\u093E\u0913", 61, -1),
        Among(u"\u0915\u0930", -1, -1, __r_CONSONANT),
        Among(u"\u0905\u0915\u0930", 64, -1),
        Among(u"\u0906\u0915\u0930", 64, -1),
        Among(u"\u093E\u0915\u0930", 64, -1),
        Among(u"\u093E", -1, -1),
        Among(u"\u090A\u0902\u0917\u093E", 68, -1),
        Among(u"\u0906\u090A\u0902\u0917\u093E", 69, -1),
        Among(u"\u093E\u090A\u0902\u0917\u093E", 69, -1),
        Among(u"\u0942\u0902\u0917\u093E", 68, -1),
        Among(u"\u090F\u0917\u093E", 68, -1),
        Among(u"\u0906\u090F\u0917\u093E", 73, -1),
        Among(u"\u093E\u090F\u0917\u093E", 73, -1),
        Among(u"\u0947\u0917\u093E", 68, -1),
        Among(u"\u0924\u093E", 68, -1, __r_CONSONANT),
        Among(u"\u0905\u0924\u093E", 77, -1),
        Among(u"\u0906\u0924\u093E", 77, -1),
        Among(u"\u093E\u0924\u093E", 77, -1),
        Among(u"\u0928\u093E", 68, -1, __r_CONSONANT),
        Among(u"\u0905\u0928\u093E", 81, -1),
        Among(u"\u0906\u0928\u093E", 81, -1),
        Among(u"\u093E\u0928\u093E", 81, -1),
        Among(u"\u0906\u092F\u093E", 68, -1),
        Among(u"\u093E\u092F\u093E", 68, -1),
        Among(u"\u093F", -1, -1),
        Among(u"\u0940", -1, -1),
        Among(u"\u090A\u0902\u0917\u0940", 88, -1),
        Among(u"\u0906\u090A\u0902\u0917\u0940", 89, -1),
        Among(u"\u093E\u090A\u0902\u0917\u0940", 89, -1),
        Among(u"\u090F\u0902\u0917\u0940", 88, -1),
        Among(u"\u0906\u090F\u0902\u0917\u0940", 92, -1),
        Among(u"\u093E\u090F\u0902\u0917\u0940", 92, -1),
        Among(u"\u0942\u0902\u0917\u0940", 88, -1),
        Among(u"\u0947\u0902\u0917\u0940", 88, -1),
        Among(u"\u090F\u0917\u0940", 88, -1),
        Among(u"\u0906\u090F\u0917\u0940", 97, -1),
        Among(u"\u093E\u090F\u0917\u0940", 97, -1),
        Among(u"\u0913\u0917\u0940", 88, -1),
        Among(u"\u0906\u0913\u0917\u0940", 100, -1),
        Among(u"\u093E\u0913\u0917\u0940", 100, -1),
        Among(u"\u0947\u0917\u0940", 88, -1),
        Among(u"\u094B\u0917\u0940", 88, -1),
        Among(u"\u0924\u0940", 88, -1, __r_CONSONANT),
        Among(u"\u0905\u0924\u0940", 105, -1),
        Among(u"\u0906\u0924\u0940", 105, -1),
        Among(u"\u093E\u0924\u0940", 105, -1),
        Among(u"\u0928\u0940", 88, -1, __r_CONSONANT),
        Among(u"\u0905\u0928\u0940", 109, -1),
        Among(u"\u0941", -1, -1),
        Among(u"\u0942", -1, -1),
        Among(u"\u0947", -1, -1),
        Among(u"\u090F\u0902\u0917\u0947", 113, -1),
        Among(u"\u0906\u090F\u0902\u0917\u0947", 114, -1),
        Among(u"\u093E\u090F\u0902\u0917\u0947", 114, -1),
        Among(u"\u0947\u0902\u0917\u0947", 113, -1),
        Among(u"\u0913\u0917\u0947", 113, -1),
        Among(u"\u0906\u0913\u0917\u0947", 118, -1),
        Among(u"\u093E\u0913\u0917\u0947", 118, -1),
        Among(u"\u094B\u0917\u0947", 113, -1),
        Among(u"\u0924\u0947", 113, -1, __r_CONSONANT),
        Among(u"\u0905\u0924\u0947", 122, -1),
        Among(u"\u0906\u0924\u0947", 122, -1),
        Among(u"\u093E\u0924\u0947", 122, -1),
        Among(u"\u0928\u0947", 113, -1, __r_CONSONANT),
        Among(u"\u0905\u0928\u0947", 126, -1),
        Among(u"\u0906\u0928\u0947", 126, -1),
        Among(u"\u093E\u0928\u0947", 126, -1),
        Among(u"\u094B", -1, -1),
        Among(u"\u094D", -1, -1)
    ]


class lab0(BaseException): pass

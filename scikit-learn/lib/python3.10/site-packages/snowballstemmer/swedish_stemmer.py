#-*- coding: utf-8 -*-
# Generated from swedish.sbl by Snowball 3.0.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class SwedishStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from swedish.sbl by Snowball 3.0.1 - https://snowballstem.org/
    '''

    g_v = {u"a", u"e", u"i", u"o", u"u", u"y", u"ä", u"å", u"ö"}

    g_s_ending = {u"b", u"c", u"d", u"f", u"g", u"h", u"j", u"k", u"l", u"m", u"n", u"o", u"p", u"r", u"t", u"v", u"y"}

    g_ost_ending = {u"i", u"k", u"l", u"n", u"p", u"r", u"t", u"u", u"v"}

    I_x = 0
    I_p1 = 0

    def __r_mark_regions(self):
        self.I_p1 = self.limit
        v_1 = self.cursor
        c = self.cursor + 3
        if c > self.limit:
            return False
        self.cursor = c
        self.I_x = self.cursor
        self.cursor = v_1
        if not self.go_out_grouping(SwedishStemmer.g_v):
            return False
        self.cursor += 1
        if not self.go_in_grouping(SwedishStemmer.g_v):
            return False
        self.cursor += 1
        self.I_p1 = self.cursor
        try:
            if self.I_p1 >= self.I_x:
                raise lab0()
            self.I_p1 = self.I_x
        except lab0: pass
        return True

    def __r_et_condition(self):
        v_1 = self.limit - self.cursor
        if not self.out_grouping_b(SwedishStemmer.g_v):
            return False
        if not self.in_grouping_b(SwedishStemmer.g_v):
            return False
        try:
            if self.cursor > self.limit_backward:
                raise lab0()
            return False
        except lab0: pass
        self.cursor = self.limit - v_1
        v_3 = self.limit - self.cursor
        try:
            if self.find_among_b(SwedishStemmer.a_0) == 0:
                raise lab1()
            return False
        except lab1: pass
        self.cursor = self.limit - v_3
        return True

    def __r_main_suffix(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        among_var = self.find_among_b(SwedishStemmer.a_1)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if among_var == 1:
            if not self.slice_del():
                return False

        elif among_var == 2:
            try:
                v_3 = self.limit - self.cursor
                try:
                    if not self.eq_s_b(u"et"):
                        raise lab1()
                    if not self.__r_et_condition():
                        raise lab1()
                    self.bra = self.cursor
                    raise lab0()
                except lab1: pass
                self.cursor = self.limit - v_3
                if not self.in_grouping_b(SwedishStemmer.g_s_ending):
                    return False
            except lab0: pass
            if not self.slice_del():
                return False

        else:
            if not self.__r_et_condition():
                return False
            if not self.slice_del():
                return False

        return True

    def __r_consonant_pair(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        v_3 = self.limit - self.cursor
        if self.find_among_b(SwedishStemmer.a_2) == 0:
            self.limit_backward = v_2
            return False
        self.cursor = self.limit - v_3
        self.ket = self.cursor
        if self.cursor <= self.limit_backward:
            self.limit_backward = v_2
            return False
        self.cursor -= 1
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.limit_backward = v_2
        return True

    def __r_other_suffix(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        among_var = self.find_among_b(SwedishStemmer.a_3)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if among_var == 1:
            if not self.slice_del():
                return False

        elif among_var == 2:
            if not self.in_grouping_b(SwedishStemmer.g_ost_ending):
                return False
            if not self.slice_from(u"ös"):
                return False
        else:
            if not self.slice_from(u"full"):
                return False
        return True

    def _stem(self):
        v_1 = self.cursor
        self.__r_mark_regions()
        self.cursor = v_1
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_2 = self.limit - self.cursor
        self.__r_main_suffix()
        self.cursor = self.limit - v_2
        v_3 = self.limit - self.cursor
        self.__r_consonant_pair()
        self.cursor = self.limit - v_3
        v_4 = self.limit - self.cursor
        self.__r_other_suffix()
        self.cursor = self.limit - v_4
        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among(u"fab", -1, -1),
        Among(u"h", -1, -1),
        Among(u"pak", -1, -1),
        Among(u"rak", -1, -1),
        Among(u"stak", -1, -1),
        Among(u"kom", -1, -1),
        Among(u"iet", -1, -1),
        Among(u"cit", -1, -1),
        Among(u"dit", -1, -1),
        Among(u"alit", -1, -1),
        Among(u"ilit", -1, -1),
        Among(u"mit", -1, -1),
        Among(u"nit", -1, -1),
        Among(u"pit", -1, -1),
        Among(u"rit", -1, -1),
        Among(u"sit", -1, -1),
        Among(u"tit", -1, -1),
        Among(u"uit", -1, -1),
        Among(u"ivit", -1, -1),
        Among(u"kvit", -1, -1),
        Among(u"xit", -1, -1)
    ]

    a_1 = [
        Among(u"a", -1, 1),
        Among(u"arna", 0, 1),
        Among(u"erna", 0, 1),
        Among(u"heterna", 2, 1),
        Among(u"orna", 0, 1),
        Among(u"ad", -1, 1),
        Among(u"e", -1, 1),
        Among(u"ade", 6, 1),
        Among(u"ande", 6, 1),
        Among(u"arne", 6, 1),
        Among(u"are", 6, 1),
        Among(u"aste", 6, 1),
        Among(u"en", -1, 1),
        Among(u"anden", 12, 1),
        Among(u"aren", 12, 1),
        Among(u"heten", 12, 1),
        Among(u"ern", -1, 1),
        Among(u"ar", -1, 1),
        Among(u"er", -1, 1),
        Among(u"heter", 18, 1),
        Among(u"or", -1, 1),
        Among(u"s", -1, 2),
        Among(u"as", 21, 1),
        Among(u"arnas", 22, 1),
        Among(u"ernas", 22, 1),
        Among(u"ornas", 22, 1),
        Among(u"es", 21, 1),
        Among(u"ades", 26, 1),
        Among(u"andes", 26, 1),
        Among(u"ens", 21, 1),
        Among(u"arens", 29, 1),
        Among(u"hetens", 29, 1),
        Among(u"erns", 21, 1),
        Among(u"at", -1, 1),
        Among(u"et", -1, 3),
        Among(u"andet", 34, 1),
        Among(u"het", 34, 1),
        Among(u"ast", -1, 1)
    ]

    a_2 = [
        Among(u"dd", -1, -1),
        Among(u"gd", -1, -1),
        Among(u"nn", -1, -1),
        Among(u"dt", -1, -1),
        Among(u"gt", -1, -1),
        Among(u"kt", -1, -1),
        Among(u"tt", -1, -1)
    ]

    a_3 = [
        Among(u"ig", -1, 1),
        Among(u"lig", 0, 1),
        Among(u"els", -1, 1),
        Among(u"fullt", -1, 3),
        Among(u"öst", -1, 2)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass

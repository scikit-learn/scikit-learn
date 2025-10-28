#-*- coding: utf-8 -*-
# Generated from norwegian.sbl by Snowball 3.0.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class NorwegianStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from norwegian.sbl by Snowball 3.0.1 - https://snowballstem.org/
    '''

    g_v = {u"a", u"e", u"ê", u"i", u"o", u"ò", u"ó", u"ô", u"u", u"y", u"æ", u"å", u"ø"}

    g_s_ending = {u"b", u"c", u"d", u"f", u"g", u"h", u"j", u"l", u"m", u"n", u"o", u"p", u"t", u"v", u"y", u"z"}

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
        if not self.go_out_grouping(NorwegianStemmer.g_v):
            return False
        self.cursor += 1
        if not self.go_in_grouping(NorwegianStemmer.g_v):
            return False
        self.cursor += 1
        self.I_p1 = self.cursor
        try:
            if self.I_p1 >= self.I_x:
                raise lab0()
            self.I_p1 = self.I_x
        except lab0: pass
        return True

    def __r_main_suffix(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        among_var = self.find_among_b(NorwegianStemmer.a_1)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if among_var == 1:
            if not self.slice_del():
                return False

        elif among_var == 2:
            among_var = self.find_among_b(NorwegianStemmer.a_0)
            if among_var == 1:
                if not self.slice_del():
                    return False

        elif among_var == 3:
            try:
                v_3 = self.limit - self.cursor
                try:
                    if not self.in_grouping_b(NorwegianStemmer.g_s_ending):
                        raise lab1()
                    raise lab0()
                except lab1: pass
                self.cursor = self.limit - v_3
                try:
                    if not self.eq_s_b(u"r"):
                        raise lab2()
                    v_4 = self.limit - self.cursor
                    try:
                        if not self.eq_s_b(u"e"):
                            raise lab3()
                        raise lab2()
                    except lab3: pass
                    self.cursor = self.limit - v_4
                    raise lab0()
                except lab2: pass
                self.cursor = self.limit - v_3
                if not self.eq_s_b(u"k"):
                    return False
                if not self.out_grouping_b(NorwegianStemmer.g_v):
                    return False
            except lab0: pass
            if not self.slice_del():
                return False

        else:
            if not self.slice_from(u"er"):
                return False
        return True

    def __r_consonant_pair(self):
        v_1 = self.limit - self.cursor
        if self.cursor < self.I_p1:
            return False
        v_3 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        if self.find_among_b(NorwegianStemmer.a_2) == 0:
            self.limit_backward = v_3
            return False
        self.bra = self.cursor
        self.limit_backward = v_3
        self.cursor = self.limit - v_1
        if self.cursor <= self.limit_backward:
            return False
        self.cursor -= 1
        self.bra = self.cursor
        if not self.slice_del():
            return False

        return True

    def __r_other_suffix(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        if self.find_among_b(NorwegianStemmer.a_3) == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if not self.slice_del():
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
        Among(u"", -1, 1),
        Among(u"ind", 0, -1),
        Among(u"kk", 0, -1),
        Among(u"nk", 0, -1),
        Among(u"amm", 0, -1),
        Among(u"omm", 0, -1),
        Among(u"kap", 0, -1),
        Among(u"skap", 6, 1),
        Among(u"pp", 0, -1),
        Among(u"lt", 0, -1),
        Among(u"ast", 0, -1),
        Among(u"øst", 0, -1),
        Among(u"v", 0, -1),
        Among(u"hav", 12, 1),
        Among(u"giv", 12, 1)
    ]

    a_1 = [
        Among(u"a", -1, 1),
        Among(u"e", -1, 1),
        Among(u"ede", 1, 1),
        Among(u"ande", 1, 1),
        Among(u"ende", 1, 1),
        Among(u"ane", 1, 1),
        Among(u"ene", 1, 1),
        Among(u"hetene", 6, 1),
        Among(u"erte", 1, 4),
        Among(u"en", -1, 1),
        Among(u"heten", 9, 1),
        Among(u"ar", -1, 1),
        Among(u"er", -1, 1),
        Among(u"heter", 12, 1),
        Among(u"s", -1, 3),
        Among(u"as", 14, 1),
        Among(u"es", 14, 1),
        Among(u"edes", 16, 1),
        Among(u"endes", 16, 1),
        Among(u"enes", 16, 1),
        Among(u"hetenes", 19, 1),
        Among(u"ens", 14, 1),
        Among(u"hetens", 21, 1),
        Among(u"ers", 14, 2),
        Among(u"ets", 14, 1),
        Among(u"et", -1, 1),
        Among(u"het", 25, 1),
        Among(u"ert", -1, 4),
        Among(u"ast", -1, 1)
    ]

    a_2 = [
        Among(u"dt", -1, -1),
        Among(u"vt", -1, -1)
    ]

    a_3 = [
        Among(u"leg", -1, 1),
        Among(u"eleg", 0, 1),
        Among(u"ig", -1, 1),
        Among(u"eig", 2, 1),
        Among(u"lig", 2, 1),
        Among(u"elig", 4, 1),
        Among(u"els", -1, 1),
        Among(u"lov", -1, 1),
        Among(u"elov", 7, 1),
        Among(u"slov", 7, 1),
        Among(u"hetslov", 9, 1)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass


class lab2(BaseException): pass


class lab3(BaseException): pass

#-*- coding: utf-8 -*-
# Generated from finnish.sbl by Snowball 3.0.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class FinnishStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from finnish.sbl by Snowball 3.0.1 - https://snowballstem.org/
    '''

    g_AEI = {u"a", u"ä", u"e", u"i"}

    g_C = {u"b", u"c", u"d", u"f", u"g", u"h", u"j", u"k", u"l", u"m", u"n", u"p", u"q", u"r", u"s", u"t", u"v", u"w", u"x", u"z"}

    g_V1 = {u"a", u"e", u"i", u"o", u"u", u"y", u"ä", u"ö"}

    g_V2 = {u"a", u"e", u"i", u"o", u"u", u"ä", u"ö"}

    g_particle_end = {u"a", u"e", u"i", u"o", u"u", u"y", u"ä", u"ö", u"n", u"t"}

    B_ending_removed = False
    S_x = ""
    I_p2 = 0
    I_p1 = 0

    def __r_mark_regions(self):
        self.I_p1 = self.limit
        self.I_p2 = self.limit
        if not self.go_out_grouping(FinnishStemmer.g_V1):
            return False
        self.cursor += 1
        if not self.go_in_grouping(FinnishStemmer.g_V1):
            return False
        self.cursor += 1
        self.I_p1 = self.cursor
        if not self.go_out_grouping(FinnishStemmer.g_V1):
            return False
        self.cursor += 1
        if not self.go_in_grouping(FinnishStemmer.g_V1):
            return False
        self.cursor += 1
        self.I_p2 = self.cursor
        return True

    def __r_R2(self):
        return self.I_p2 <= self.cursor

    def __r_particle_etc(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        among_var = self.find_among_b(FinnishStemmer.a_0)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if among_var == 1:
            if not self.in_grouping_b(FinnishStemmer.g_particle_end):
                return False
        else:
            if not self.__r_R2():
                return False
        if not self.slice_del():
            return False

        return True

    def __r_possessive(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        among_var = self.find_among_b(FinnishStemmer.a_4)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if among_var == 1:
            v_3 = self.limit - self.cursor
            try:
                if not self.eq_s_b(u"k"):
                    raise lab0()
                return False
            except lab0: pass
            self.cursor = self.limit - v_3
            if not self.slice_del():
                return False

        elif among_var == 2:
            if not self.slice_del():
                return False

            self.ket = self.cursor
            if not self.eq_s_b(u"kse"):
                return False
            self.bra = self.cursor
            if not self.slice_from(u"ksi"):
                return False
        elif among_var == 3:
            if not self.slice_del():
                return False

        elif among_var == 4:
            if self.find_among_b(FinnishStemmer.a_1) == 0:
                return False
            if not self.slice_del():
                return False

        elif among_var == 5:
            if self.find_among_b(FinnishStemmer.a_2) == 0:
                return False
            if not self.slice_del():
                return False

        else:
            if self.find_among_b(FinnishStemmer.a_3) == 0:
                return False
            if not self.slice_del():
                return False

        return True

    def __r_LONG(self):
        if self.find_among_b(FinnishStemmer.a_5) == 0:
            return False
        return True

    def __r_VI(self):
        if not self.eq_s_b(u"i"):
            return False
        if not self.in_grouping_b(FinnishStemmer.g_V2):
            return False
        return True

    def __r_case_ending(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        among_var = self.find_among_b(FinnishStemmer.a_6)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if among_var == 1:
            if not self.eq_s_b(u"a"):
                return False
        elif among_var == 2:
            if not self.eq_s_b(u"e"):
                return False
        elif among_var == 3:
            if not self.eq_s_b(u"i"):
                return False
        elif among_var == 4:
            if not self.eq_s_b(u"o"):
                return False
        elif among_var == 5:
            if not self.eq_s_b(u"ä"):
                return False
        elif among_var == 6:
            if not self.eq_s_b(u"ö"):
                return False
        elif among_var == 7:
            v_3 = self.limit - self.cursor
            try:
                v_4 = self.limit - self.cursor
                try:
                    v_5 = self.limit - self.cursor
                    try:
                        if not self.__r_LONG():
                            raise lab2()
                        raise lab1()
                    except lab2: pass
                    self.cursor = self.limit - v_5
                    if not self.eq_s_b(u"ie"):
                        self.cursor = self.limit - v_3
                        raise lab0()
                except lab1: pass
                self.cursor = self.limit - v_4
                if self.cursor <= self.limit_backward:
                    self.cursor = self.limit - v_3
                    raise lab0()
                self.cursor -= 1
                self.bra = self.cursor
            except lab0: pass
        elif among_var == 8:
            if not self.in_grouping_b(FinnishStemmer.g_V1):
                return False
            if not self.in_grouping_b(FinnishStemmer.g_C):
                return False
        if not self.slice_del():
            return False

        self.B_ending_removed = True
        return True

    def __r_other_endings(self):
        if self.cursor < self.I_p2:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p2
        self.ket = self.cursor
        among_var = self.find_among_b(FinnishStemmer.a_7)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if among_var == 1:
            v_3 = self.limit - self.cursor
            try:
                if not self.eq_s_b(u"po"):
                    raise lab0()
                return False
            except lab0: pass
            self.cursor = self.limit - v_3
        if not self.slice_del():
            return False

        return True

    def __r_i_plural(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        if self.find_among_b(FinnishStemmer.a_8) == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if not self.slice_del():
            return False

        return True

    def __r_t_plural(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        if not self.eq_s_b(u"t"):
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        v_3 = self.limit - self.cursor
        if not self.in_grouping_b(FinnishStemmer.g_V1):
            self.limit_backward = v_2
            return False
        self.cursor = self.limit - v_3
        if not self.slice_del():
            return False

        self.limit_backward = v_2
        if self.cursor < self.I_p2:
            return False
        v_5 = self.limit_backward
        self.limit_backward = self.I_p2
        self.ket = self.cursor
        among_var = self.find_among_b(FinnishStemmer.a_9)
        if among_var == 0:
            self.limit_backward = v_5
            return False
        self.bra = self.cursor
        self.limit_backward = v_5
        if among_var == 1:
            v_6 = self.limit - self.cursor
            try:
                if not self.eq_s_b(u"po"):
                    raise lab0()
                return False
            except lab0: pass
            self.cursor = self.limit - v_6
        if not self.slice_del():
            return False

        return True

    def __r_tidy(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        v_3 = self.limit - self.cursor
        try:
            v_4 = self.limit - self.cursor
            if not self.__r_LONG():
                raise lab0()
            self.cursor = self.limit - v_4
            self.ket = self.cursor
            if self.cursor <= self.limit_backward:
                raise lab0()
            self.cursor -= 1
            self.bra = self.cursor
            if not self.slice_del():
                return False

        except lab0: pass
        self.cursor = self.limit - v_3
        v_5 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if not self.in_grouping_b(FinnishStemmer.g_AEI):
                raise lab1()
            self.bra = self.cursor
            if not self.in_grouping_b(FinnishStemmer.g_C):
                raise lab1()
            if not self.slice_del():
                return False

        except lab1: pass
        self.cursor = self.limit - v_5
        v_6 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if not self.eq_s_b(u"j"):
                raise lab2()
            self.bra = self.cursor
            try:
                v_7 = self.limit - self.cursor
                try:
                    if not self.eq_s_b(u"o"):
                        raise lab4()
                    raise lab3()
                except lab4: pass
                self.cursor = self.limit - v_7
                if not self.eq_s_b(u"u"):
                    raise lab2()
            except lab3: pass
            if not self.slice_del():
                return False

        except lab2: pass
        self.cursor = self.limit - v_6
        v_8 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if not self.eq_s_b(u"o"):
                raise lab5()
            self.bra = self.cursor
            if not self.eq_s_b(u"j"):
                raise lab5()
            if not self.slice_del():
                return False

        except lab5: pass
        self.cursor = self.limit - v_8
        self.limit_backward = v_2
        if not self.go_in_grouping_b(FinnishStemmer.g_V1):
            return False
        self.ket = self.cursor
        if not self.in_grouping_b(FinnishStemmer.g_C):
            return False
        self.bra = self.cursor
        self.S_x = self.slice_to()
        if self.S_x == '':
            return False
        if not self.eq_s_b(self.S_x):
            return False
        if not self.slice_del():
            return False

        return True

    def _stem(self):
        v_1 = self.cursor
        self.__r_mark_regions()
        self.cursor = v_1
        self.B_ending_removed = False
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_2 = self.limit - self.cursor
        self.__r_particle_etc()
        self.cursor = self.limit - v_2
        v_3 = self.limit - self.cursor
        self.__r_possessive()
        self.cursor = self.limit - v_3
        v_4 = self.limit - self.cursor
        self.__r_case_ending()
        self.cursor = self.limit - v_4
        v_5 = self.limit - self.cursor
        self.__r_other_endings()
        self.cursor = self.limit - v_5
        try:
            try:
                if not self.B_ending_removed:
                    raise lab1()
                v_7 = self.limit - self.cursor
                self.__r_i_plural()
                self.cursor = self.limit - v_7
                raise lab0()
            except lab1: pass
            v_8 = self.limit - self.cursor
            self.__r_t_plural()
            self.cursor = self.limit - v_8
        except lab0: pass
        v_9 = self.limit - self.cursor
        self.__r_tidy()
        self.cursor = self.limit - v_9
        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among(u"pa", -1, 1),
        Among(u"sti", -1, 2),
        Among(u"kaan", -1, 1),
        Among(u"han", -1, 1),
        Among(u"kin", -1, 1),
        Among(u"hän", -1, 1),
        Among(u"kään", -1, 1),
        Among(u"ko", -1, 1),
        Among(u"pä", -1, 1),
        Among(u"kö", -1, 1)
    ]

    a_1 = [
        Among(u"lla", -1, -1),
        Among(u"na", -1, -1),
        Among(u"ssa", -1, -1),
        Among(u"ta", -1, -1),
        Among(u"lta", 3, -1),
        Among(u"sta", 3, -1)
    ]

    a_2 = [
        Among(u"llä", -1, -1),
        Among(u"nä", -1, -1),
        Among(u"ssä", -1, -1),
        Among(u"tä", -1, -1),
        Among(u"ltä", 3, -1),
        Among(u"stä", 3, -1)
    ]

    a_3 = [
        Among(u"lle", -1, -1),
        Among(u"ine", -1, -1)
    ]

    a_4 = [
        Among(u"nsa", -1, 3),
        Among(u"mme", -1, 3),
        Among(u"nne", -1, 3),
        Among(u"ni", -1, 2),
        Among(u"si", -1, 1),
        Among(u"an", -1, 4),
        Among(u"en", -1, 6),
        Among(u"än", -1, 5),
        Among(u"nsä", -1, 3)
    ]

    a_5 = [
        Among(u"aa", -1, -1),
        Among(u"ee", -1, -1),
        Among(u"ii", -1, -1),
        Among(u"oo", -1, -1),
        Among(u"uu", -1, -1),
        Among(u"ää", -1, -1),
        Among(u"öö", -1, -1)
    ]

    a_6 = [
        Among(u"a", -1, 8),
        Among(u"lla", 0, -1),
        Among(u"na", 0, -1),
        Among(u"ssa", 0, -1),
        Among(u"ta", 0, -1),
        Among(u"lta", 4, -1),
        Among(u"sta", 4, -1),
        Among(u"tta", 4, 2),
        Among(u"lle", -1, -1),
        Among(u"ine", -1, -1),
        Among(u"ksi", -1, -1),
        Among(u"n", -1, 7),
        Among(u"han", 11, 1),
        Among(u"den", 11, -1, __r_VI),
        Among(u"seen", 11, -1, __r_LONG),
        Among(u"hen", 11, 2),
        Among(u"tten", 11, -1, __r_VI),
        Among(u"hin", 11, 3),
        Among(u"siin", 11, -1, __r_VI),
        Among(u"hon", 11, 4),
        Among(u"hän", 11, 5),
        Among(u"hön", 11, 6),
        Among(u"ä", -1, 8),
        Among(u"llä", 22, -1),
        Among(u"nä", 22, -1),
        Among(u"ssä", 22, -1),
        Among(u"tä", 22, -1),
        Among(u"ltä", 26, -1),
        Among(u"stä", 26, -1),
        Among(u"ttä", 26, 2)
    ]

    a_7 = [
        Among(u"eja", -1, -1),
        Among(u"mma", -1, 1),
        Among(u"imma", 1, -1),
        Among(u"mpa", -1, 1),
        Among(u"impa", 3, -1),
        Among(u"mmi", -1, 1),
        Among(u"immi", 5, -1),
        Among(u"mpi", -1, 1),
        Among(u"impi", 7, -1),
        Among(u"ejä", -1, -1),
        Among(u"mmä", -1, 1),
        Among(u"immä", 10, -1),
        Among(u"mpä", -1, 1),
        Among(u"impä", 12, -1)
    ]

    a_8 = [
        Among(u"i", -1, -1),
        Among(u"j", -1, -1)
    ]

    a_9 = [
        Among(u"mma", -1, 1),
        Among(u"imma", 0, -1)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass


class lab2(BaseException): pass


class lab3(BaseException): pass


class lab4(BaseException): pass


class lab5(BaseException): pass

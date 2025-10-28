#-*- coding: utf-8 -*-
# Generated from esperanto.sbl by Snowball 3.0.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class EsperantoStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from esperanto.sbl by Snowball 3.0.1 - https://snowballstem.org/
    '''

    g_vowel = {u"a", u"e", u"i", u"o", u"u"}

    g_aou = {u"a", u"o", u"u"}

    g_digit = {u"0", u"1", u"2", u"3", u"4", u"5", u"6", u"7", u"8", u"9"}

    B_foreign = False

    def __r_canonical_form(self):
        self.B_foreign = False
        while True:
            v_1 = self.cursor
            try:
                self.bra = self.cursor
                among_var = self.find_among(EsperantoStemmer.a_0)
                self.ket = self.cursor
                if among_var == 1:
                    if not self.slice_from(u"ĉ"):
                        return False
                elif among_var == 2:
                    if not self.slice_from(u"ĝ"):
                        return False
                elif among_var == 3:
                    if not self.slice_from(u"ĥ"):
                        return False
                elif among_var == 4:
                    if not self.slice_from(u"ĵ"):
                        return False
                elif among_var == 5:
                    if not self.slice_from(u"ŝ"):
                        return False
                elif among_var == 6:
                    if not self.slice_from(u"ŭ"):
                        return False
                elif among_var == 7:
                    if not self.slice_from(u"a"):
                        return False
                    self.B_foreign = True
                elif among_var == 8:
                    if not self.slice_from(u"e"):
                        return False
                    self.B_foreign = True
                elif among_var == 9:
                    if not self.slice_from(u"i"):
                        return False
                    self.B_foreign = True
                elif among_var == 10:
                    if not self.slice_from(u"o"):
                        return False
                    self.B_foreign = True
                elif among_var == 11:
                    if not self.slice_from(u"u"):
                        return False
                    self.B_foreign = True
                elif among_var == 12:
                    self.B_foreign = True
                elif among_var == 13:
                    self.B_foreign = False
                else:
                    if self.cursor >= self.limit:
                        raise lab0()
                    self.cursor += 1
                continue
            except lab0: pass
            self.cursor = v_1
            break
        try:
            if not self.B_foreign:
                raise lab1()
            return False
        except lab1: pass
        return True

    def __r_initial_apostrophe(self):
        self.bra = self.cursor
        if not self.eq_s(u"'"):
            return False
        self.ket = self.cursor
        if not self.eq_s(u"st"):
            return False
        if self.find_among(EsperantoStemmer.a_1) == 0:
            return False
        if self.cursor < self.limit:
            return False
        if not self.slice_from(u"e"):
            return False
        return True

    def __r_pronoun(self):
        self.ket = self.cursor
        v_1 = self.limit - self.cursor
        try:
            if not self.eq_s_b(u"n"):
                self.cursor = self.limit - v_1
                raise lab0()
        except lab0: pass
        self.bra = self.cursor
        if self.find_among_b(EsperantoStemmer.a_2) == 0:
            return False
        try:
            v_2 = self.limit - self.cursor
            try:
                if self.cursor > self.limit_backward:
                    raise lab2()
                raise lab1()
            except lab2: pass
            self.cursor = self.limit - v_2
            if not self.eq_s_b(u"-"):
                return False
        except lab1: pass
        if not self.slice_del():
            return False

        return True

    def __r_final_apostrophe(self):
        self.ket = self.cursor
        if not self.eq_s_b(u"'"):
            return False
        self.bra = self.cursor
        try:
            v_1 = self.limit - self.cursor
            try:
                if not self.eq_s_b(u"l"):
                    raise lab1()
                if self.cursor > self.limit_backward:
                    raise lab1()
                if not self.slice_from(u"a"):
                    return False
                raise lab0()
            except lab1: pass
            self.cursor = self.limit - v_1
            try:
                if not self.eq_s_b(u"un"):
                    raise lab2()
                if self.cursor > self.limit_backward:
                    raise lab2()
                if not self.slice_from(u"u"):
                    return False
                raise lab0()
            except lab2: pass
            self.cursor = self.limit - v_1
            try:
                if self.find_among_b(EsperantoStemmer.a_3) == 0:
                    raise lab3()
                try:
                    v_2 = self.limit - self.cursor
                    try:
                        if self.cursor > self.limit_backward:
                            raise lab5()
                        raise lab4()
                    except lab5: pass
                    self.cursor = self.limit - v_2
                    if not self.eq_s_b(u"-"):
                        raise lab3()
                except lab4: pass
                if not self.slice_from(u"aŭ"):
                    return False
                raise lab0()
            except lab3: pass
            self.cursor = self.limit - v_1
            if not self.slice_from(u"o"):
                return False
        except lab0: pass
        return True

    def __r_ujn_suffix(self):
        self.ket = self.cursor
        v_1 = self.limit - self.cursor
        try:
            if not self.eq_s_b(u"n"):
                self.cursor = self.limit - v_1
                raise lab0()
        except lab0: pass
        v_2 = self.limit - self.cursor
        try:
            if not self.eq_s_b(u"j"):
                self.cursor = self.limit - v_2
                raise lab1()
        except lab1: pass
        self.bra = self.cursor
        if self.find_among_b(EsperantoStemmer.a_4) == 0:
            return False
        try:
            v_3 = self.limit - self.cursor
            try:
                if self.cursor > self.limit_backward:
                    raise lab3()
                raise lab2()
            except lab3: pass
            self.cursor = self.limit - v_3
            if not self.eq_s_b(u"-"):
                return False
        except lab2: pass
        if not self.slice_del():
            return False

        return True

    def __r_uninflected(self):
        if self.find_among_b(EsperantoStemmer.a_5) == 0:
            return False
        try:
            v_1 = self.limit - self.cursor
            try:
                if self.cursor > self.limit_backward:
                    raise lab1()
                raise lab0()
            except lab1: pass
            self.cursor = self.limit - v_1
            if not self.eq_s_b(u"-"):
                return False
        except lab0: pass
        return True

    def __r_merged_numeral(self):
        if self.find_among_b(EsperantoStemmer.a_6) == 0:
            return False
        if self.find_among_b(EsperantoStemmer.a_7) == 0:
            return False
        return True

    def __r_correlative(self):
        self.ket = self.cursor
        self.bra = self.cursor
        v_1 = self.limit - self.cursor
        try:
            v_2 = self.limit - self.cursor
            try:
                v_3 = self.limit - self.cursor
                try:
                    if not self.eq_s_b(u"n"):
                        self.cursor = self.limit - v_3
                        raise lab2()
                except lab2: pass
                self.bra = self.cursor
                if not self.eq_s_b(u"e"):
                    raise lab1()
                raise lab0()
            except lab1: pass
            self.cursor = self.limit - v_2
            v_4 = self.limit - self.cursor
            try:
                if not self.eq_s_b(u"n"):
                    self.cursor = self.limit - v_4
                    raise lab3()
            except lab3: pass
            v_5 = self.limit - self.cursor
            try:
                if not self.eq_s_b(u"j"):
                    self.cursor = self.limit - v_5
                    raise lab4()
            except lab4: pass
            self.bra = self.cursor
            if not self.in_grouping_b(EsperantoStemmer.g_aou):
                return False
        except lab0: pass
        if not self.eq_s_b(u"i"):
            return False
        v_6 = self.limit - self.cursor
        try:
            if self.find_among_b(EsperantoStemmer.a_8) == 0:
                self.cursor = self.limit - v_6
                raise lab5()
        except lab5: pass
        try:
            v_7 = self.limit - self.cursor
            try:
                if self.cursor > self.limit_backward:
                    raise lab7()
                raise lab6()
            except lab7: pass
            self.cursor = self.limit - v_7
            if not self.eq_s_b(u"-"):
                return False
        except lab6: pass
        self.cursor = self.limit - v_1
        if not self.slice_del():
            return False

        return True

    def __r_long_word(self):
        try:
            v_1 = self.limit - self.cursor
            try:
                for v_2 in 0, 0:

                    if not self.go_out_grouping_b(EsperantoStemmer.g_vowel):
                        raise lab1()
                    self.cursor -= 1
                raise lab0()
            except lab1: pass
            self.cursor = self.limit - v_1
            try:
                try:
                    while True:
                        try:
                            if not self.eq_s_b(u"-"):
                                raise lab4()
                            raise lab3()
                        except lab4: pass
                        if self.cursor <= self.limit_backward:
                            raise lab2()
                        self.cursor -= 1
                except lab3: pass
                if self.cursor <= self.limit_backward:
                    raise lab2()
                self.cursor -= 1
                raise lab0()
            except lab2: pass
            self.cursor = self.limit - v_1
            if not self.go_out_grouping_b(EsperantoStemmer.g_digit):
                return False
            self.cursor -= 1
        except lab0: pass
        return True

    def __r_not_after_letter(self):
        try:
            v_1 = self.limit - self.cursor
            try:
                if not self.eq_s_b(u"-"):
                    raise lab1()
                raise lab0()
            except lab1: pass
            self.cursor = self.limit - v_1
            if not self.in_grouping_b(EsperantoStemmer.g_digit):
                return False
        except lab0: pass
        return True

    def __r_standard_suffix(self):
        self.ket = self.cursor
        if self.find_among_b(EsperantoStemmer.a_9) == 0:
            return False
        v_1 = self.limit - self.cursor
        try:
            if not self.eq_s_b(u"-"):
                self.cursor = self.limit - v_1
                raise lab0()
        except lab0: pass
        self.bra = self.cursor
        if not self.slice_del():
            return False

        return True

    def _stem(self):
        v_1 = self.cursor
        if not self.__r_canonical_form():
            return False
        self.cursor = v_1
        v_2 = self.cursor
        self.__r_initial_apostrophe()
        self.cursor = v_2
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_3 = self.limit - self.cursor
        try:
            if not self.__r_pronoun():
                raise lab0()
            return False
        except lab0: pass
        self.cursor = self.limit - v_3
        v_4 = self.limit - self.cursor
        self.__r_final_apostrophe()
        self.cursor = self.limit - v_4
        v_5 = self.limit - self.cursor
        try:
            if not self.__r_correlative():
                raise lab1()
            return False
        except lab1: pass
        self.cursor = self.limit - v_5
        v_6 = self.limit - self.cursor
        try:
            if not self.__r_uninflected():
                raise lab2()
            return False
        except lab2: pass
        self.cursor = self.limit - v_6
        v_7 = self.limit - self.cursor
        try:
            if not self.__r_merged_numeral():
                raise lab3()
            return False
        except lab3: pass
        self.cursor = self.limit - v_7
        v_8 = self.limit - self.cursor
        try:
            if not self.__r_ujn_suffix():
                raise lab4()
            return False
        except lab4: pass
        self.cursor = self.limit - v_8
        v_9 = self.limit - self.cursor
        if not self.__r_long_word():
            return False
        self.cursor = self.limit - v_9
        if not self.__r_standard_suffix():
            return False
        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among(u"", -1, 14),
        Among(u"-", 0, 13),
        Among(u"cx", 0, 1),
        Among(u"gx", 0, 2),
        Among(u"hx", 0, 3),
        Among(u"jx", 0, 4),
        Among(u"q", 0, 12),
        Among(u"sx", 0, 5),
        Among(u"ux", 0, 6),
        Among(u"w", 0, 12),
        Among(u"x", 0, 12),
        Among(u"y", 0, 12),
        Among(u"á", 0, 7),
        Among(u"é", 0, 8),
        Among(u"í", 0, 9),
        Among(u"ó", 0, 10),
        Among(u"ú", 0, 11)
    ]

    a_1 = [
        Among(u"as", -1, -1),
        Among(u"i", -1, -1),
        Among(u"is", 1, -1),
        Among(u"os", -1, -1),
        Among(u"u", -1, -1),
        Among(u"us", 4, -1)
    ]

    a_2 = [
        Among(u"ci", -1, -1),
        Among(u"gi", -1, -1),
        Among(u"hi", -1, -1),
        Among(u"li", -1, -1),
        Among(u"ili", 3, -1),
        Among(u"ŝli", 3, -1),
        Among(u"mi", -1, -1),
        Among(u"ni", -1, -1),
        Among(u"oni", 7, -1),
        Among(u"ri", -1, -1),
        Among(u"si", -1, -1),
        Among(u"vi", -1, -1),
        Among(u"ivi", 11, -1),
        Among(u"ĝi", -1, -1),
        Among(u"ŝi", -1, -1),
        Among(u"iŝi", 14, -1),
        Among(u"malŝi", 14, -1)
    ]

    a_3 = [
        Among(u"amb", -1, -1),
        Among(u"bald", -1, -1),
        Among(u"malbald", 1, -1),
        Among(u"morg", -1, -1),
        Among(u"postmorg", 3, -1),
        Among(u"adi", -1, -1),
        Among(u"hodi", -1, -1),
        Among(u"ank", -1, -1),
        Among(u"ĉirk", -1, -1),
        Among(u"tutĉirk", 8, -1),
        Among(u"presk", -1, -1),
        Among(u"almen", -1, -1),
        Among(u"apen", -1, -1),
        Among(u"hier", -1, -1),
        Among(u"antaŭhier", 13, -1),
        Among(u"malgr", -1, -1),
        Among(u"ankor", -1, -1),
        Among(u"kontr", -1, -1),
        Among(u"anstat", -1, -1),
        Among(u"kvaz", -1, -1)
    ]

    a_4 = [
        Among(u"aliu", -1, -1),
        Among(u"unu", -1, -1)
    ]

    a_5 = [
        Among(u"aha", -1, -1),
        Among(u"haha", 0, -1),
        Among(u"haleluja", -1, -1),
        Among(u"hola", -1, -1),
        Among(u"hosana", -1, -1),
        Among(u"maltra", -1, -1),
        Among(u"hura", -1, -1),
        Among(u"ĥaĥa", -1, -1),
        Among(u"ekde", -1, -1),
        Among(u"elde", -1, -1),
        Among(u"disde", -1, -1),
        Among(u"ehe", -1, -1),
        Among(u"maltre", -1, -1),
        Among(u"dirlididi", -1, -1),
        Among(u"malpli", -1, -1),
        Among(u"malĉi", -1, -1),
        Among(u"malkaj", -1, -1),
        Among(u"amen", -1, -1),
        Among(u"tamen", 17, -1),
        Among(u"oho", -1, -1),
        Among(u"maltro", -1, -1),
        Among(u"minus", -1, -1),
        Among(u"uhu", -1, -1),
        Among(u"muu", -1, -1)
    ]

    a_6 = [
        Among(u"tri", -1, -1),
        Among(u"du", -1, -1),
        Among(u"unu", -1, -1)
    ]

    a_7 = [
        Among(u"dek", -1, -1),
        Among(u"cent", -1, -1)
    ]

    a_8 = [
        Among(u"k", -1, -1),
        Among(u"kelk", 0, -1),
        Among(u"nen", -1, -1),
        Among(u"t", -1, -1),
        Among(u"mult", 3, -1),
        Among(u"samt", 3, -1),
        Among(u"ĉ", -1, -1)
    ]

    a_9 = [
        Among(u"a", -1, -1),
        Among(u"e", -1, -1),
        Among(u"i", -1, -1),
        Among(u"j", -1, -1, __r_not_after_letter),
        Among(u"aj", 3, -1),
        Among(u"oj", 3, -1),
        Among(u"n", -1, -1, __r_not_after_letter),
        Among(u"an", 6, -1),
        Among(u"en", 6, -1),
        Among(u"jn", 6, -1, __r_not_after_letter),
        Among(u"ajn", 9, -1),
        Among(u"ojn", 9, -1),
        Among(u"on", 6, -1),
        Among(u"o", -1, -1),
        Among(u"as", -1, -1),
        Among(u"is", -1, -1),
        Among(u"os", -1, -1),
        Among(u"us", -1, -1),
        Among(u"u", -1, -1)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass


class lab2(BaseException): pass


class lab3(BaseException): pass


class lab4(BaseException): pass


class lab5(BaseException): pass


class lab6(BaseException): pass


class lab7(BaseException): pass

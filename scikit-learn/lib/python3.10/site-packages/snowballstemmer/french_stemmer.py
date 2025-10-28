#-*- coding: utf-8 -*-
# Generated from french.sbl by Snowball 3.0.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class FrenchStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from french.sbl by Snowball 3.0.1 - https://snowballstem.org/
    '''

    g_v = {u"a", u"e", u"i", u"o", u"u", u"y", u"â", u"à", u"ë", u"é", u"ê", u"è", u"ï", u"î", u"ô", u"û", u"ù"}

    g_oux_ending = {u"b", u"h", u"j", u"l", u"n", u"p"}

    g_elision_char = {u"c", u"d", u"j", u"l", u"m", u"n", u"s", u"t"}

    g_keep_with_s = {u"a", u"i", u"o", u"u", u"è", u"s"}

    I_p2 = 0
    I_p1 = 0
    I_pV = 0

    def __r_elisions(self):
        self.bra = self.cursor
        try:
            v_1 = self.cursor
            try:
                if not self.in_grouping(FrenchStemmer.g_elision_char):
                    raise lab1()
                raise lab0()
            except lab1: pass
            self.cursor = v_1
            if not self.eq_s(u"qu"):
                return False
        except lab0: pass
        if not self.eq_s(u"'"):
            return False
        self.ket = self.cursor
        try:
            if self.cursor < self.limit:
                raise lab2()
            return False
        except lab2: pass
        if not self.slice_del():
            return False

        return True

    def __r_prelude(self):
        while True:
            v_1 = self.cursor
            try:
                try:
                    while True:
                        v_2 = self.cursor
                        try:
                            try:
                                v_3 = self.cursor
                                try:
                                    if not self.in_grouping(FrenchStemmer.g_v):
                                        raise lab4()
                                    self.bra = self.cursor
                                    try:
                                        v_4 = self.cursor
                                        try:
                                            if not self.eq_s(u"u"):
                                                raise lab6()
                                            self.ket = self.cursor
                                            if not self.in_grouping(FrenchStemmer.g_v):
                                                raise lab6()
                                            if not self.slice_from(u"U"):
                                                return False
                                            raise lab5()
                                        except lab6: pass
                                        self.cursor = v_4
                                        try:
                                            if not self.eq_s(u"i"):
                                                raise lab7()
                                            self.ket = self.cursor
                                            if not self.in_grouping(FrenchStemmer.g_v):
                                                raise lab7()
                                            if not self.slice_from(u"I"):
                                                return False
                                            raise lab5()
                                        except lab7: pass
                                        self.cursor = v_4
                                        if not self.eq_s(u"y"):
                                            raise lab4()
                                        self.ket = self.cursor
                                        if not self.slice_from(u"Y"):
                                            return False
                                    except lab5: pass
                                    raise lab3()
                                except lab4: pass
                                self.cursor = v_3
                                try:
                                    self.bra = self.cursor
                                    if not self.eq_s(u"ë"):
                                        raise lab8()
                                    self.ket = self.cursor
                                    if not self.slice_from(u"He"):
                                        return False
                                    raise lab3()
                                except lab8: pass
                                self.cursor = v_3
                                try:
                                    self.bra = self.cursor
                                    if not self.eq_s(u"ï"):
                                        raise lab9()
                                    self.ket = self.cursor
                                    if not self.slice_from(u"Hi"):
                                        return False
                                    raise lab3()
                                except lab9: pass
                                self.cursor = v_3
                                try:
                                    self.bra = self.cursor
                                    if not self.eq_s(u"y"):
                                        raise lab10()
                                    self.ket = self.cursor
                                    if not self.in_grouping(FrenchStemmer.g_v):
                                        raise lab10()
                                    if not self.slice_from(u"Y"):
                                        return False
                                    raise lab3()
                                except lab10: pass
                                self.cursor = v_3
                                if not self.eq_s(u"q"):
                                    raise lab2()
                                self.bra = self.cursor
                                if not self.eq_s(u"u"):
                                    raise lab2()
                                self.ket = self.cursor
                                if not self.slice_from(u"U"):
                                    return False
                            except lab3: pass
                            self.cursor = v_2
                            raise lab1()
                        except lab2: pass
                        self.cursor = v_2
                        if self.cursor >= self.limit:
                            raise lab0()
                        self.cursor += 1
                except lab1: pass
                continue
            except lab0: pass
            self.cursor = v_1
            break
        return True

    def __r_mark_regions(self):
        self.I_pV = self.limit
        self.I_p1 = self.limit
        self.I_p2 = self.limit
        v_1 = self.cursor
        try:
            try:
                v_2 = self.cursor
                try:
                    if not self.in_grouping(FrenchStemmer.g_v):
                        raise lab2()
                    if not self.in_grouping(FrenchStemmer.g_v):
                        raise lab2()
                    if self.cursor >= self.limit:
                        raise lab2()
                    self.cursor += 1
                    raise lab1()
                except lab2: pass
                self.cursor = v_2
                try:
                    among_var = self.find_among(FrenchStemmer.a_0)
                    if among_var == 0:
                        raise lab3()
                    if among_var == 1:
                        if not self.in_grouping(FrenchStemmer.g_v):
                            raise lab3()
                    raise lab1()
                except lab3: pass
                self.cursor = v_2
                if self.cursor >= self.limit:
                    raise lab0()
                self.cursor += 1
                if not self.go_out_grouping(FrenchStemmer.g_v):
                    raise lab0()
                self.cursor += 1
            except lab1: pass
            self.I_pV = self.cursor
        except lab0: pass
        self.cursor = v_1
        v_3 = self.cursor
        try:
            if not self.go_out_grouping(FrenchStemmer.g_v):
                raise lab4()
            self.cursor += 1
            if not self.go_in_grouping(FrenchStemmer.g_v):
                raise lab4()
            self.cursor += 1
            self.I_p1 = self.cursor
            if not self.go_out_grouping(FrenchStemmer.g_v):
                raise lab4()
            self.cursor += 1
            if not self.go_in_grouping(FrenchStemmer.g_v):
                raise lab4()
            self.cursor += 1
            self.I_p2 = self.cursor
        except lab4: pass
        self.cursor = v_3
        return True

    def __r_postlude(self):
        while True:
            v_1 = self.cursor
            try:
                self.bra = self.cursor
                among_var = self.find_among(FrenchStemmer.a_1)
                self.ket = self.cursor
                if among_var == 1:
                    if not self.slice_from(u"i"):
                        return False
                elif among_var == 2:
                    if not self.slice_from(u"u"):
                        return False
                elif among_var == 3:
                    if not self.slice_from(u"y"):
                        return False
                elif among_var == 4:
                    if not self.slice_from(u"ë"):
                        return False
                elif among_var == 5:
                    if not self.slice_from(u"ï"):
                        return False
                elif among_var == 6:
                    if not self.slice_del():
                        return False

                else:
                    if self.cursor >= self.limit:
                        raise lab0()
                    self.cursor += 1
                continue
            except lab0: pass
            self.cursor = v_1
            break
        return True

    def __r_RV(self):
        return self.I_pV <= self.cursor

    def __r_R1(self):
        return self.I_p1 <= self.cursor

    def __r_R2(self):
        return self.I_p2 <= self.cursor

    def __r_standard_suffix(self):
        self.ket = self.cursor
        among_var = self.find_among_b(FrenchStemmer.a_4)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if not self.__r_R2():
                return False
            if not self.slice_del():
                return False

        elif among_var == 2:
            if not self.__r_R2():
                return False
            if not self.slice_del():
                return False

            v_1 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                if not self.eq_s_b(u"ic"):
                    self.cursor = self.limit - v_1
                    raise lab0()
                self.bra = self.cursor
                try:
                    v_2 = self.limit - self.cursor
                    try:
                        if not self.__r_R2():
                            raise lab2()
                        if not self.slice_del():
                            return False

                        raise lab1()
                    except lab2: pass
                    self.cursor = self.limit - v_2
                    if not self.slice_from(u"iqU"):
                        return False
                except lab1: pass
            except lab0: pass
        elif among_var == 3:
            if not self.__r_R2():
                return False
            if not self.slice_from(u"log"):
                return False
        elif among_var == 4:
            if not self.__r_R2():
                return False
            if not self.slice_from(u"u"):
                return False
        elif among_var == 5:
            if not self.__r_R2():
                return False
            if not self.slice_from(u"ent"):
                return False
        elif among_var == 6:
            if not self.__r_RV():
                return False
            if not self.slice_del():
                return False

            v_3 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                among_var = self.find_among_b(FrenchStemmer.a_2)
                if among_var == 0:
                    self.cursor = self.limit - v_3
                    raise lab3()
                self.bra = self.cursor
                if among_var == 1:
                    if not self.__r_R2():
                        self.cursor = self.limit - v_3
                        raise lab3()
                    if not self.slice_del():
                        return False

                    self.ket = self.cursor
                    if not self.eq_s_b(u"at"):
                        self.cursor = self.limit - v_3
                        raise lab3()
                    self.bra = self.cursor
                    if not self.__r_R2():
                        self.cursor = self.limit - v_3
                        raise lab3()
                    if not self.slice_del():
                        return False

                elif among_var == 2:
                    try:
                        v_4 = self.limit - self.cursor
                        try:
                            if not self.__r_R2():
                                raise lab5()
                            if not self.slice_del():
                                return False

                            raise lab4()
                        except lab5: pass
                        self.cursor = self.limit - v_4
                        if not self.__r_R1():
                            self.cursor = self.limit - v_3
                            raise lab3()
                        if not self.slice_from(u"eux"):
                            return False
                    except lab4: pass
                elif among_var == 3:
                    if not self.__r_R2():
                        self.cursor = self.limit - v_3
                        raise lab3()
                    if not self.slice_del():
                        return False

                else:
                    if not self.__r_RV():
                        self.cursor = self.limit - v_3
                        raise lab3()
                    if not self.slice_from(u"i"):
                        return False
            except lab3: pass
        elif among_var == 7:
            if not self.__r_R2():
                return False
            if not self.slice_del():
                return False

            v_5 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                among_var = self.find_among_b(FrenchStemmer.a_3)
                if among_var == 0:
                    self.cursor = self.limit - v_5
                    raise lab6()
                self.bra = self.cursor
                if among_var == 1:
                    try:
                        v_6 = self.limit - self.cursor
                        try:
                            if not self.__r_R2():
                                raise lab8()
                            if not self.slice_del():
                                return False

                            raise lab7()
                        except lab8: pass
                        self.cursor = self.limit - v_6
                        if not self.slice_from(u"abl"):
                            return False
                    except lab7: pass
                elif among_var == 2:
                    try:
                        v_7 = self.limit - self.cursor
                        try:
                            if not self.__r_R2():
                                raise lab10()
                            if not self.slice_del():
                                return False

                            raise lab9()
                        except lab10: pass
                        self.cursor = self.limit - v_7
                        if not self.slice_from(u"iqU"):
                            return False
                    except lab9: pass
                else:
                    if not self.__r_R2():
                        self.cursor = self.limit - v_5
                        raise lab6()
                    if not self.slice_del():
                        return False

            except lab6: pass
        elif among_var == 8:
            if not self.__r_R2():
                return False
            if not self.slice_del():
                return False

            v_8 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                if not self.eq_s_b(u"at"):
                    self.cursor = self.limit - v_8
                    raise lab11()
                self.bra = self.cursor
                if not self.__r_R2():
                    self.cursor = self.limit - v_8
                    raise lab11()
                if not self.slice_del():
                    return False

                self.ket = self.cursor
                if not self.eq_s_b(u"ic"):
                    self.cursor = self.limit - v_8
                    raise lab11()
                self.bra = self.cursor
                try:
                    v_9 = self.limit - self.cursor
                    try:
                        if not self.__r_R2():
                            raise lab13()
                        if not self.slice_del():
                            return False

                        raise lab12()
                    except lab13: pass
                    self.cursor = self.limit - v_9
                    if not self.slice_from(u"iqU"):
                        return False
                except lab12: pass
            except lab11: pass
        elif among_var == 9:
            if not self.slice_from(u"eau"):
                return False
        elif among_var == 10:
            if not self.__r_R1():
                return False
            if not self.slice_from(u"al"):
                return False
        elif among_var == 11:
            if not self.in_grouping_b(FrenchStemmer.g_oux_ending):
                return False
            if not self.slice_from(u"ou"):
                return False
        elif among_var == 12:
            try:
                v_10 = self.limit - self.cursor
                try:
                    if not self.__r_R2():
                        raise lab15()
                    if not self.slice_del():
                        return False

                    raise lab14()
                except lab15: pass
                self.cursor = self.limit - v_10
                if not self.__r_R1():
                    return False
                if not self.slice_from(u"eux"):
                    return False
            except lab14: pass
        elif among_var == 13:
            if not self.__r_R1():
                return False
            if not self.out_grouping_b(FrenchStemmer.g_v):
                return False
            if not self.slice_del():
                return False

        elif among_var == 14:
            if not self.__r_RV():
                return False
            if not self.slice_from(u"ant"):
                return False
            return False
        elif among_var == 15:
            if not self.__r_RV():
                return False
            if not self.slice_from(u"ent"):
                return False
            return False
        else:
            v_11 = self.limit - self.cursor
            if not self.in_grouping_b(FrenchStemmer.g_v):
                return False
            if not self.__r_RV():
                return False
            self.cursor = self.limit - v_11
            if not self.slice_del():
                return False

            return False
        return True

    def __r_i_verb_suffix(self):
        if self.cursor < self.I_pV:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_pV
        self.ket = self.cursor
        if self.find_among_b(FrenchStemmer.a_5) == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        v_3 = self.limit - self.cursor
        try:
            if not self.eq_s_b(u"H"):
                raise lab0()
            self.limit_backward = v_2
            return False
        except lab0: pass
        self.cursor = self.limit - v_3
        if not self.out_grouping_b(FrenchStemmer.g_v):
            self.limit_backward = v_2
            return False
        if not self.slice_del():
            return False

        self.limit_backward = v_2
        return True

    def __r_verb_suffix(self):
        if self.cursor < self.I_pV:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_pV
        self.ket = self.cursor
        among_var = self.find_among_b(FrenchStemmer.a_7)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if among_var == 1:
            if not self.__r_R2():
                return False
            if not self.slice_del():
                return False

        elif among_var == 2:
            if not self.slice_del():
                return False

        elif among_var == 3:
            v_3 = self.limit - self.cursor
            try:
                if not self.eq_s_b(u"e"):
                    self.cursor = self.limit - v_3
                    raise lab0()
                if not self.__r_RV():
                    self.cursor = self.limit - v_3
                    raise lab0()
                self.bra = self.cursor
            except lab0: pass
            if not self.slice_del():
                return False

        else:
            v_4 = self.limit - self.cursor
            try:
                among_var = self.find_among_b(FrenchStemmer.a_6)
                if among_var == 0:
                    raise lab1()
                if among_var == 1:
                    if self.cursor <= self.limit_backward:
                        raise lab1()
                    self.cursor -= 1
                    if self.cursor > self.limit_backward:
                        raise lab1()
                return False
            except lab1: pass
            self.cursor = self.limit - v_4
            if not self.slice_del():
                return False

        return True

    def __r_residual_suffix(self):
        v_1 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if not self.eq_s_b(u"s"):
                self.cursor = self.limit - v_1
                raise lab0()
            self.bra = self.cursor
            v_2 = self.limit - self.cursor
            try:
                v_3 = self.limit - self.cursor
                try:
                    if not self.eq_s_b(u"Hi"):
                        raise lab2()
                    raise lab1()
                except lab2: pass
                self.cursor = self.limit - v_3
                if not self.out_grouping_b(FrenchStemmer.g_keep_with_s):
                    self.cursor = self.limit - v_1
                    raise lab0()
            except lab1: pass
            self.cursor = self.limit - v_2
            if not self.slice_del():
                return False

        except lab0: pass
        if self.cursor < self.I_pV:
            return False
        v_5 = self.limit_backward
        self.limit_backward = self.I_pV
        self.ket = self.cursor
        among_var = self.find_among_b(FrenchStemmer.a_8)
        if among_var == 0:
            self.limit_backward = v_5
            return False
        self.bra = self.cursor
        if among_var == 1:
            if not self.__r_R2():
                self.limit_backward = v_5
                return False
            try:
                v_6 = self.limit - self.cursor
                try:
                    if not self.eq_s_b(u"s"):
                        raise lab4()
                    raise lab3()
                except lab4: pass
                self.cursor = self.limit - v_6
                if not self.eq_s_b(u"t"):
                    self.limit_backward = v_5
                    return False
            except lab3: pass
            if not self.slice_del():
                return False

        elif among_var == 2:
            if not self.slice_from(u"i"):
                return False
        else:
            if not self.slice_del():
                return False

        self.limit_backward = v_5
        return True

    def __r_un_double(self):
        v_1 = self.limit - self.cursor
        if self.find_among_b(FrenchStemmer.a_9) == 0:
            return False
        self.cursor = self.limit - v_1
        self.ket = self.cursor
        if self.cursor <= self.limit_backward:
            return False
        self.cursor -= 1
        self.bra = self.cursor
        if not self.slice_del():
            return False

        return True

    def __r_un_accent(self):
        v_1 = 1
        while True:
            try:
                if not self.out_grouping_b(FrenchStemmer.g_v):
                    raise lab0()
                v_1 -= 1
                continue
            except lab0: pass
            break
        if v_1 > 0:
            return False
        self.ket = self.cursor
        try:
            v_3 = self.limit - self.cursor
            try:
                if not self.eq_s_b(u"é"):
                    raise lab2()
                raise lab1()
            except lab2: pass
            self.cursor = self.limit - v_3
            if not self.eq_s_b(u"è"):
                return False
        except lab1: pass
        self.bra = self.cursor
        if not self.slice_from(u"e"):
            return False
        return True

    def _stem(self):
        v_1 = self.cursor
        self.__r_elisions()
        self.cursor = v_1
        v_2 = self.cursor
        self.__r_prelude()
        self.cursor = v_2
        self.__r_mark_regions()
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_4 = self.limit - self.cursor
        try:
            try:
                v_5 = self.limit - self.cursor
                try:
                    v_6 = self.limit - self.cursor
                    try:
                        v_7 = self.limit - self.cursor
                        try:
                            if not self.__r_standard_suffix():
                                raise lab4()
                            raise lab3()
                        except lab4: pass
                        self.cursor = self.limit - v_7
                        try:
                            if not self.__r_i_verb_suffix():
                                raise lab5()
                            raise lab3()
                        except lab5: pass
                        self.cursor = self.limit - v_7
                        if not self.__r_verb_suffix():
                            raise lab2()
                    except lab3: pass
                    self.cursor = self.limit - v_6
                    v_8 = self.limit - self.cursor
                    try:
                        self.ket = self.cursor
                        try:
                            v_9 = self.limit - self.cursor
                            try:
                                if not self.eq_s_b(u"Y"):
                                    raise lab8()
                                self.bra = self.cursor
                                if not self.slice_from(u"i"):
                                    return False
                                raise lab7()
                            except lab8: pass
                            self.cursor = self.limit - v_9
                            if not self.eq_s_b(u"ç"):
                                self.cursor = self.limit - v_8
                                raise lab6()
                            self.bra = self.cursor
                            if not self.slice_from(u"c"):
                                return False
                        except lab7: pass
                    except lab6: pass
                    raise lab1()
                except lab2: pass
                self.cursor = self.limit - v_5
                if not self.__r_residual_suffix():
                    raise lab0()
            except lab1: pass
        except lab0: pass
        self.cursor = self.limit - v_4
        v_10 = self.limit - self.cursor
        self.__r_un_double()
        self.cursor = self.limit - v_10
        v_11 = self.limit - self.cursor
        self.__r_un_accent()
        self.cursor = self.limit - v_11
        self.cursor = self.limit_backward
        v_12 = self.cursor
        self.__r_postlude()
        self.cursor = v_12
        return True

    a_0 = [
        Among(u"col", -1, -1),
        Among(u"ni", -1, 1),
        Among(u"par", -1, -1),
        Among(u"tap", -1, -1)
    ]

    a_1 = [
        Among(u"", -1, 7),
        Among(u"H", 0, 6),
        Among(u"He", 1, 4),
        Among(u"Hi", 1, 5),
        Among(u"I", 0, 1),
        Among(u"U", 0, 2),
        Among(u"Y", 0, 3)
    ]

    a_2 = [
        Among(u"iqU", -1, 3),
        Among(u"abl", -1, 3),
        Among(u"Ièr", -1, 4),
        Among(u"ièr", -1, 4),
        Among(u"eus", -1, 2),
        Among(u"iv", -1, 1)
    ]

    a_3 = [
        Among(u"ic", -1, 2),
        Among(u"abil", -1, 1),
        Among(u"iv", -1, 3)
    ]

    a_4 = [
        Among(u"iqUe", -1, 1),
        Among(u"atrice", -1, 2),
        Among(u"ance", -1, 1),
        Among(u"ence", -1, 5),
        Among(u"logie", -1, 3),
        Among(u"able", -1, 1),
        Among(u"isme", -1, 1),
        Among(u"euse", -1, 12),
        Among(u"iste", -1, 1),
        Among(u"ive", -1, 8),
        Among(u"if", -1, 8),
        Among(u"usion", -1, 4),
        Among(u"ation", -1, 2),
        Among(u"ution", -1, 4),
        Among(u"ateur", -1, 2),
        Among(u"iqUes", -1, 1),
        Among(u"atrices", -1, 2),
        Among(u"ances", -1, 1),
        Among(u"ences", -1, 5),
        Among(u"logies", -1, 3),
        Among(u"ables", -1, 1),
        Among(u"ismes", -1, 1),
        Among(u"euses", -1, 12),
        Among(u"istes", -1, 1),
        Among(u"ives", -1, 8),
        Among(u"ifs", -1, 8),
        Among(u"usions", -1, 4),
        Among(u"ations", -1, 2),
        Among(u"utions", -1, 4),
        Among(u"ateurs", -1, 2),
        Among(u"ments", -1, 16),
        Among(u"ements", 30, 6),
        Among(u"issements", 31, 13),
        Among(u"ités", -1, 7),
        Among(u"ment", -1, 16),
        Among(u"ement", 34, 6),
        Among(u"issement", 35, 13),
        Among(u"amment", 34, 14),
        Among(u"emment", 34, 15),
        Among(u"aux", -1, 10),
        Among(u"eaux", 39, 9),
        Among(u"eux", -1, 1),
        Among(u"oux", -1, 11),
        Among(u"ité", -1, 7)
    ]

    a_5 = [
        Among(u"ira", -1, 1),
        Among(u"ie", -1, 1),
        Among(u"isse", -1, 1),
        Among(u"issante", -1, 1),
        Among(u"i", -1, 1),
        Among(u"irai", 4, 1),
        Among(u"ir", -1, 1),
        Among(u"iras", -1, 1),
        Among(u"ies", -1, 1),
        Among(u"îmes", -1, 1),
        Among(u"isses", -1, 1),
        Among(u"issantes", -1, 1),
        Among(u"îtes", -1, 1),
        Among(u"is", -1, 1),
        Among(u"irais", 13, 1),
        Among(u"issais", 13, 1),
        Among(u"irions", -1, 1),
        Among(u"issions", -1, 1),
        Among(u"irons", -1, 1),
        Among(u"issons", -1, 1),
        Among(u"issants", -1, 1),
        Among(u"it", -1, 1),
        Among(u"irait", 21, 1),
        Among(u"issait", 21, 1),
        Among(u"issant", -1, 1),
        Among(u"iraIent", -1, 1),
        Among(u"issaIent", -1, 1),
        Among(u"irent", -1, 1),
        Among(u"issent", -1, 1),
        Among(u"iront", -1, 1),
        Among(u"ît", -1, 1),
        Among(u"iriez", -1, 1),
        Among(u"issiez", -1, 1),
        Among(u"irez", -1, 1),
        Among(u"issez", -1, 1)
    ]

    a_6 = [
        Among(u"al", -1, 1),
        Among(u"épl", -1, -1),
        Among(u"auv", -1, -1)
    ]

    a_7 = [
        Among(u"a", -1, 3),
        Among(u"era", 0, 2),
        Among(u"aise", -1, 4),
        Among(u"asse", -1, 3),
        Among(u"ante", -1, 3),
        Among(u"ée", -1, 2),
        Among(u"ai", -1, 3),
        Among(u"erai", 6, 2),
        Among(u"er", -1, 2),
        Among(u"as", -1, 3),
        Among(u"eras", 9, 2),
        Among(u"âmes", -1, 3),
        Among(u"aises", -1, 4),
        Among(u"asses", -1, 3),
        Among(u"antes", -1, 3),
        Among(u"âtes", -1, 3),
        Among(u"ées", -1, 2),
        Among(u"ais", -1, 4),
        Among(u"eais", 17, 2),
        Among(u"erais", 17, 2),
        Among(u"ions", -1, 1),
        Among(u"erions", 20, 2),
        Among(u"assions", 20, 3),
        Among(u"erons", -1, 2),
        Among(u"ants", -1, 3),
        Among(u"és", -1, 2),
        Among(u"ait", -1, 3),
        Among(u"erait", 26, 2),
        Among(u"ant", -1, 3),
        Among(u"aIent", -1, 3),
        Among(u"eraIent", 29, 2),
        Among(u"èrent", -1, 2),
        Among(u"assent", -1, 3),
        Among(u"eront", -1, 2),
        Among(u"ât", -1, 3),
        Among(u"ez", -1, 2),
        Among(u"iez", 35, 2),
        Among(u"eriez", 36, 2),
        Among(u"assiez", 36, 3),
        Among(u"erez", 35, 2),
        Among(u"é", -1, 2)
    ]

    a_8 = [
        Among(u"e", -1, 3),
        Among(u"Ière", 0, 2),
        Among(u"ière", 0, 2),
        Among(u"ion", -1, 1),
        Among(u"Ier", -1, 2),
        Among(u"ier", -1, 2)
    ]

    a_9 = [
        Among(u"ell", -1, -1),
        Among(u"eill", -1, -1),
        Among(u"enn", -1, -1),
        Among(u"onn", -1, -1),
        Among(u"ett", -1, -1)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass


class lab2(BaseException): pass


class lab3(BaseException): pass


class lab4(BaseException): pass


class lab5(BaseException): pass


class lab6(BaseException): pass


class lab7(BaseException): pass


class lab8(BaseException): pass


class lab9(BaseException): pass


class lab10(BaseException): pass


class lab11(BaseException): pass


class lab12(BaseException): pass


class lab13(BaseException): pass


class lab14(BaseException): pass


class lab15(BaseException): pass

#-*- coding: utf-8 -*-
# Generated from english.sbl by Snowball 3.0.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class EnglishStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from english.sbl by Snowball 3.0.1 - https://snowballstem.org/
    '''

    g_aeo = {u"a", u"e", u"o"}

    g_v = {u"a", u"e", u"i", u"o", u"u", u"y"}

    g_v_WXY = {u"a", u"e", u"i", u"o", u"u", u"y", u"w", u"x", u"Y"}

    g_valid_LI = {u"c", u"d", u"e", u"g", u"h", u"k", u"m", u"n", u"r", u"t"}

    B_Y_found = False
    I_p2 = 0
    I_p1 = 0

    def __r_prelude(self):
        self.B_Y_found = False
        v_1 = self.cursor
        try:
            self.bra = self.cursor
            if not self.eq_s(u"'"):
                raise lab0()
            self.ket = self.cursor
            if not self.slice_del():
                return False

        except lab0: pass
        self.cursor = v_1
        v_2 = self.cursor
        try:
            self.bra = self.cursor
            if not self.eq_s(u"y"):
                raise lab1()
            self.ket = self.cursor
            if not self.slice_from(u"Y"):
                return False
            self.B_Y_found = True
        except lab1: pass
        self.cursor = v_2
        v_3 = self.cursor
        try:
            while True:
                v_4 = self.cursor
                try:
                    try:
                        while True:
                            v_5 = self.cursor
                            try:
                                if not self.in_grouping(EnglishStemmer.g_v):
                                    raise lab5()
                                self.bra = self.cursor
                                if not self.eq_s(u"y"):
                                    raise lab5()
                                self.ket = self.cursor
                                self.cursor = v_5
                                raise lab4()
                            except lab5: pass
                            self.cursor = v_5
                            if self.cursor >= self.limit:
                                raise lab3()
                            self.cursor += 1
                    except lab4: pass
                    if not self.slice_from(u"Y"):
                        return False
                    self.B_Y_found = True
                    continue
                except lab3: pass
                self.cursor = v_4
                break
        except lab2: pass
        self.cursor = v_3
        return True

    def __r_mark_regions(self):
        self.I_p1 = self.limit
        self.I_p2 = self.limit
        v_1 = self.cursor
        try:
            try:
                v_2 = self.cursor
                try:
                    if self.find_among(EnglishStemmer.a_0) == 0:
                        raise lab2()
                    raise lab1()
                except lab2: pass
                self.cursor = v_2
                if not self.go_out_grouping(EnglishStemmer.g_v):
                    raise lab0()
                self.cursor += 1
                if not self.go_in_grouping(EnglishStemmer.g_v):
                    raise lab0()
                self.cursor += 1
            except lab1: pass
            self.I_p1 = self.cursor
            if not self.go_out_grouping(EnglishStemmer.g_v):
                raise lab0()
            self.cursor += 1
            if not self.go_in_grouping(EnglishStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_p2 = self.cursor
        except lab0: pass
        self.cursor = v_1
        return True

    def __r_shortv(self):
        try:
            v_1 = self.limit - self.cursor
            try:
                if not self.out_grouping_b(EnglishStemmer.g_v_WXY):
                    raise lab1()
                if not self.in_grouping_b(EnglishStemmer.g_v):
                    raise lab1()
                if not self.out_grouping_b(EnglishStemmer.g_v):
                    raise lab1()
                raise lab0()
            except lab1: pass
            self.cursor = self.limit - v_1
            try:
                if not self.out_grouping_b(EnglishStemmer.g_v):
                    raise lab2()
                if not self.in_grouping_b(EnglishStemmer.g_v):
                    raise lab2()
                if self.cursor > self.limit_backward:
                    raise lab2()
                raise lab0()
            except lab2: pass
            self.cursor = self.limit - v_1
            if not self.eq_s_b(u"past"):
                return False
        except lab0: pass
        return True

    def __r_R1(self):
        return self.I_p1 <= self.cursor

    def __r_R2(self):
        return self.I_p2 <= self.cursor

    def __r_Step_1a(self):
        v_1 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if self.find_among_b(EnglishStemmer.a_1) == 0:
                self.cursor = self.limit - v_1
                raise lab0()
            self.bra = self.cursor
            if not self.slice_del():
                return False

        except lab0: pass
        self.ket = self.cursor
        among_var = self.find_among_b(EnglishStemmer.a_2)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if not self.slice_from(u"ss"):
                return False
        elif among_var == 2:
            try:
                v_2 = self.limit - self.cursor
                try:
                    c = self.cursor - 2
                    if c < self.limit_backward:
                        raise lab2()
                    self.cursor = c
                    if not self.slice_from(u"i"):
                        return False
                    raise lab1()
                except lab2: pass
                self.cursor = self.limit - v_2
                if not self.slice_from(u"ie"):
                    return False
            except lab1: pass
        elif among_var == 3:
            if self.cursor <= self.limit_backward:
                return False
            self.cursor -= 1
            if not self.go_out_grouping_b(EnglishStemmer.g_v):
                return False
            self.cursor -= 1
            if not self.slice_del():
                return False

        return True

    def __r_Step_1b(self):
        self.ket = self.cursor
        among_var = self.find_among_b(EnglishStemmer.a_5)
        self.bra = self.cursor
        try:
            v_1 = self.limit - self.cursor
            try:
                if among_var == 1:
                    v_2 = self.limit - self.cursor
                    try:
                        try:
                            v_3 = self.limit - self.cursor
                            try:
                                if self.find_among_b(EnglishStemmer.a_3) == 0:
                                    raise lab4()
                                if self.cursor > self.limit_backward:
                                    raise lab4()
                                raise lab3()
                            except lab4: pass
                            self.cursor = self.limit - v_3
                            if not self.__r_R1():
                                raise lab2()
                            if not self.slice_from(u"ee"):
                                return False
                        except lab3: pass
                    except lab2: pass
                    self.cursor = self.limit - v_2
                elif among_var == 2:
                    raise lab1()
                elif among_var == 3:
                    among_var = self.find_among_b(EnglishStemmer.a_4)
                    if among_var == 0:
                        raise lab1()
                    if among_var == 1:
                        v_4 = self.limit - self.cursor
                        if not self.out_grouping_b(EnglishStemmer.g_v):
                            raise lab1()
                        if self.cursor > self.limit_backward:
                            raise lab1()
                        self.cursor = self.limit - v_4
                        self.bra = self.cursor
                        if not self.slice_from(u"ie"):
                            return False
                    else:
                        if self.cursor > self.limit_backward:
                            raise lab1()
                raise lab0()
            except lab1: pass
            self.cursor = self.limit - v_1
            v_5 = self.limit - self.cursor
            if not self.go_out_grouping_b(EnglishStemmer.g_v):
                return False
            self.cursor -= 1
            self.cursor = self.limit - v_5
            if not self.slice_del():
                return False

            self.ket = self.cursor
            self.bra = self.cursor
            v_6 = self.limit - self.cursor
            among_var = self.find_among_b(EnglishStemmer.a_6)
            if among_var == 1:
                if not self.slice_from(u"e"):
                    return False
                return False
            elif among_var == 2:
                v_7 = self.limit - self.cursor
                try:
                    if not self.in_grouping_b(EnglishStemmer.g_aeo):
                        raise lab5()
                    if self.cursor > self.limit_backward:
                        raise lab5()
                    return False
                except lab5: pass
                self.cursor = self.limit - v_7
            else:
                if self.cursor != self.I_p1:
                    return False
                v_8 = self.limit - self.cursor
                if not self.__r_shortv():
                    return False
                self.cursor = self.limit - v_8
                if not self.slice_from(u"e"):
                    return False
                return False
            self.cursor = self.limit - v_6
            self.ket = self.cursor
            if self.cursor <= self.limit_backward:
                return False
            self.cursor -= 1
            self.bra = self.cursor
            if not self.slice_del():
                return False

        except lab0: pass
        return True

    def __r_Step_1c(self):
        self.ket = self.cursor
        try:
            v_1 = self.limit - self.cursor
            try:
                if not self.eq_s_b(u"y"):
                    raise lab1()
                raise lab0()
            except lab1: pass
            self.cursor = self.limit - v_1
            if not self.eq_s_b(u"Y"):
                return False
        except lab0: pass
        self.bra = self.cursor
        if not self.out_grouping_b(EnglishStemmer.g_v):
            return False
        try:
            if self.cursor > self.limit_backward:
                raise lab2()
            return False
        except lab2: pass
        if not self.slice_from(u"i"):
            return False
        return True

    def __r_Step_2(self):
        self.ket = self.cursor
        among_var = self.find_among_b(EnglishStemmer.a_7)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R1():
            return False
        if among_var == 1:
            if not self.slice_from(u"tion"):
                return False
        elif among_var == 2:
            if not self.slice_from(u"ence"):
                return False
        elif among_var == 3:
            if not self.slice_from(u"ance"):
                return False
        elif among_var == 4:
            if not self.slice_from(u"able"):
                return False
        elif among_var == 5:
            if not self.slice_from(u"ent"):
                return False
        elif among_var == 6:
            if not self.slice_from(u"ize"):
                return False
        elif among_var == 7:
            if not self.slice_from(u"ate"):
                return False
        elif among_var == 8:
            if not self.slice_from(u"al"):
                return False
        elif among_var == 9:
            if not self.slice_from(u"ful"):
                return False
        elif among_var == 10:
            if not self.slice_from(u"ous"):
                return False
        elif among_var == 11:
            if not self.slice_from(u"ive"):
                return False
        elif among_var == 12:
            if not self.slice_from(u"ble"):
                return False
        elif among_var == 13:
            if not self.slice_from(u"og"):
                return False
        elif among_var == 14:
            if not self.eq_s_b(u"l"):
                return False
            if not self.slice_from(u"og"):
                return False
        elif among_var == 15:
            if not self.slice_from(u"less"):
                return False
        else:
            if not self.in_grouping_b(EnglishStemmer.g_valid_LI):
                return False
            if not self.slice_del():
                return False

        return True

    def __r_Step_3(self):
        self.ket = self.cursor
        among_var = self.find_among_b(EnglishStemmer.a_8)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R1():
            return False
        if among_var == 1:
            if not self.slice_from(u"tion"):
                return False
        elif among_var == 2:
            if not self.slice_from(u"ate"):
                return False
        elif among_var == 3:
            if not self.slice_from(u"al"):
                return False
        elif among_var == 4:
            if not self.slice_from(u"ic"):
                return False
        elif among_var == 5:
            if not self.slice_del():
                return False

        else:
            if not self.__r_R2():
                return False
            if not self.slice_del():
                return False

        return True

    def __r_Step_4(self):
        self.ket = self.cursor
        among_var = self.find_among_b(EnglishStemmer.a_9)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R2():
            return False
        if among_var == 1:
            if not self.slice_del():
                return False

        else:
            try:
                v_1 = self.limit - self.cursor
                try:
                    if not self.eq_s_b(u"s"):
                        raise lab1()
                    raise lab0()
                except lab1: pass
                self.cursor = self.limit - v_1
                if not self.eq_s_b(u"t"):
                    return False
            except lab0: pass
            if not self.slice_del():
                return False

        return True

    def __r_Step_5(self):
        self.ket = self.cursor
        among_var = self.find_among_b(EnglishStemmer.a_10)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            try:
                try:
                    if not self.__r_R2():
                        raise lab1()
                    raise lab0()
                except lab1: pass
                if not self.__r_R1():
                    return False
                v_2 = self.limit - self.cursor
                try:
                    if not self.__r_shortv():
                        raise lab2()
                    return False
                except lab2: pass
                self.cursor = self.limit - v_2
            except lab0: pass
            if not self.slice_del():
                return False

        else:
            if not self.__r_R2():
                return False
            if not self.eq_s_b(u"l"):
                return False
            if not self.slice_del():
                return False

        return True

    def __r_exception1(self):
        self.bra = self.cursor
        among_var = self.find_among(EnglishStemmer.a_11)
        if among_var == 0:
            return False
        self.ket = self.cursor
        if self.cursor < self.limit:
            return False
        if among_var == 1:
            if not self.slice_from(u"sky"):
                return False
        elif among_var == 2:
            if not self.slice_from(u"idl"):
                return False
        elif among_var == 3:
            if not self.slice_from(u"gentl"):
                return False
        elif among_var == 4:
            if not self.slice_from(u"ugli"):
                return False
        elif among_var == 5:
            if not self.slice_from(u"earli"):
                return False
        elif among_var == 6:
            if not self.slice_from(u"onli"):
                return False
        elif among_var == 7:
            if not self.slice_from(u"singl"):
                return False
        return True

    def __r_postlude(self):
        if not self.B_Y_found:
            return False
        while True:
            v_1 = self.cursor
            try:
                try:
                    while True:
                        v_2 = self.cursor
                        try:
                            self.bra = self.cursor
                            if not self.eq_s(u"Y"):
                                raise lab2()
                            self.ket = self.cursor
                            self.cursor = v_2
                            raise lab1()
                        except lab2: pass
                        self.cursor = v_2
                        if self.cursor >= self.limit:
                            raise lab0()
                        self.cursor += 1
                except lab1: pass
                if not self.slice_from(u"y"):
                    return False
                continue
            except lab0: pass
            self.cursor = v_1
            break
        return True

    def _stem(self):
        try:
            v_1 = self.cursor
            try:
                if not self.__r_exception1():
                    raise lab1()
                raise lab0()
            except lab1: pass
            self.cursor = v_1
            try:
                v_2 = self.cursor
                try:
                    c = self.cursor + 3
                    if c > self.limit:
                        raise lab3()
                    self.cursor = c
                    raise lab2()
                except lab3: pass
                self.cursor = v_2
                raise lab0()
            except lab2: pass
            self.cursor = v_1
            self.__r_prelude()
            self.__r_mark_regions()
            self.limit_backward = self.cursor
            self.cursor = self.limit
            v_5 = self.limit - self.cursor
            self.__r_Step_1a()
            self.cursor = self.limit - v_5
            v_6 = self.limit - self.cursor
            self.__r_Step_1b()
            self.cursor = self.limit - v_6
            v_7 = self.limit - self.cursor
            self.__r_Step_1c()
            self.cursor = self.limit - v_7
            v_8 = self.limit - self.cursor
            self.__r_Step_2()
            self.cursor = self.limit - v_8
            v_9 = self.limit - self.cursor
            self.__r_Step_3()
            self.cursor = self.limit - v_9
            v_10 = self.limit - self.cursor
            self.__r_Step_4()
            self.cursor = self.limit - v_10
            v_11 = self.limit - self.cursor
            self.__r_Step_5()
            self.cursor = self.limit - v_11
            self.cursor = self.limit_backward
            v_12 = self.cursor
            self.__r_postlude()
            self.cursor = v_12
        except lab0: pass
        return True

    a_0 = [
        Among(u"arsen", -1, -1),
        Among(u"commun", -1, -1),
        Among(u"emerg", -1, -1),
        Among(u"gener", -1, -1),
        Among(u"later", -1, -1),
        Among(u"organ", -1, -1),
        Among(u"past", -1, -1),
        Among(u"univers", -1, -1)
    ]

    a_1 = [
        Among(u"'", -1, 1),
        Among(u"'s'", 0, 1),
        Among(u"'s", -1, 1)
    ]

    a_2 = [
        Among(u"ied", -1, 2),
        Among(u"s", -1, 3),
        Among(u"ies", 1, 2),
        Among(u"sses", 1, 1),
        Among(u"ss", 1, -1),
        Among(u"us", 1, -1)
    ]

    a_3 = [
        Among(u"succ", -1, 1),
        Among(u"proc", -1, 1),
        Among(u"exc", -1, 1)
    ]

    a_4 = [
        Among(u"even", -1, 2),
        Among(u"cann", -1, 2),
        Among(u"inn", -1, 2),
        Among(u"earr", -1, 2),
        Among(u"herr", -1, 2),
        Among(u"out", -1, 2),
        Among(u"y", -1, 1)
    ]

    a_5 = [
        Among(u"", -1, -1),
        Among(u"ed", 0, 2),
        Among(u"eed", 1, 1),
        Among(u"ing", 0, 3),
        Among(u"edly", 0, 2),
        Among(u"eedly", 4, 1),
        Among(u"ingly", 0, 2)
    ]

    a_6 = [
        Among(u"", -1, 3),
        Among(u"bb", 0, 2),
        Among(u"dd", 0, 2),
        Among(u"ff", 0, 2),
        Among(u"gg", 0, 2),
        Among(u"bl", 0, 1),
        Among(u"mm", 0, 2),
        Among(u"nn", 0, 2),
        Among(u"pp", 0, 2),
        Among(u"rr", 0, 2),
        Among(u"at", 0, 1),
        Among(u"tt", 0, 2),
        Among(u"iz", 0, 1)
    ]

    a_7 = [
        Among(u"anci", -1, 3),
        Among(u"enci", -1, 2),
        Among(u"ogi", -1, 14),
        Among(u"li", -1, 16),
        Among(u"bli", 3, 12),
        Among(u"abli", 4, 4),
        Among(u"alli", 3, 8),
        Among(u"fulli", 3, 9),
        Among(u"lessli", 3, 15),
        Among(u"ousli", 3, 10),
        Among(u"entli", 3, 5),
        Among(u"aliti", -1, 8),
        Among(u"biliti", -1, 12),
        Among(u"iviti", -1, 11),
        Among(u"tional", -1, 1),
        Among(u"ational", 14, 7),
        Among(u"alism", -1, 8),
        Among(u"ation", -1, 7),
        Among(u"ization", 17, 6),
        Among(u"izer", -1, 6),
        Among(u"ator", -1, 7),
        Among(u"iveness", -1, 11),
        Among(u"fulness", -1, 9),
        Among(u"ousness", -1, 10),
        Among(u"ogist", -1, 13)
    ]

    a_8 = [
        Among(u"icate", -1, 4),
        Among(u"ative", -1, 6),
        Among(u"alize", -1, 3),
        Among(u"iciti", -1, 4),
        Among(u"ical", -1, 4),
        Among(u"tional", -1, 1),
        Among(u"ational", 5, 2),
        Among(u"ful", -1, 5),
        Among(u"ness", -1, 5)
    ]

    a_9 = [
        Among(u"ic", -1, 1),
        Among(u"ance", -1, 1),
        Among(u"ence", -1, 1),
        Among(u"able", -1, 1),
        Among(u"ible", -1, 1),
        Among(u"ate", -1, 1),
        Among(u"ive", -1, 1),
        Among(u"ize", -1, 1),
        Among(u"iti", -1, 1),
        Among(u"al", -1, 1),
        Among(u"ism", -1, 1),
        Among(u"ion", -1, 2),
        Among(u"er", -1, 1),
        Among(u"ous", -1, 1),
        Among(u"ant", -1, 1),
        Among(u"ent", -1, 1),
        Among(u"ment", 15, 1),
        Among(u"ement", 16, 1)
    ]

    a_10 = [
        Among(u"e", -1, 1),
        Among(u"l", -1, 2)
    ]

    a_11 = [
        Among(u"andes", -1, -1),
        Among(u"atlas", -1, -1),
        Among(u"bias", -1, -1),
        Among(u"cosmos", -1, -1),
        Among(u"early", -1, 5),
        Among(u"gently", -1, 3),
        Among(u"howe", -1, -1),
        Among(u"idly", -1, 2),
        Among(u"news", -1, -1),
        Among(u"only", -1, 6),
        Among(u"singly", -1, 7),
        Among(u"skies", -1, 1),
        Among(u"sky", -1, -1),
        Among(u"ugly", -1, 4)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass


class lab2(BaseException): pass


class lab3(BaseException): pass


class lab4(BaseException): pass


class lab5(BaseException): pass

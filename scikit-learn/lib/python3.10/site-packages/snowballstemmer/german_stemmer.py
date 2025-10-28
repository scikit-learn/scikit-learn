#-*- coding: utf-8 -*-
# Generated from german.sbl by Snowball 3.0.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class GermanStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from german.sbl by Snowball 3.0.1 - https://snowballstem.org/
    '''

    g_v = {u"a", u"e", u"i", u"o", u"u", u"y", u"ä", u"ö", u"ü"}

    g_et_ending = {u"d", u"f", u"g", u"k", u"l", u"m", u"n", u"r", u"s", u"t", u"U", u"z", u"ä"}

    g_s_ending = {u"b", u"d", u"f", u"g", u"h", u"k", u"l", u"m", u"n", u"r", u"t"}

    g_st_ending = {u"b", u"d", u"f", u"g", u"h", u"k", u"l", u"m", u"n", u"t"}

    I_x = 0
    I_p2 = 0
    I_p1 = 0

    def __r_prelude(self):
        v_1 = self.cursor
        while True:
            v_2 = self.cursor
            try:
                try:
                    while True:
                        v_3 = self.cursor
                        try:
                            if not self.in_grouping(GermanStemmer.g_v):
                                raise lab2()
                            self.bra = self.cursor
                            try:
                                v_4 = self.cursor
                                try:
                                    if not self.eq_s(u"u"):
                                        raise lab4()
                                    self.ket = self.cursor
                                    if not self.in_grouping(GermanStemmer.g_v):
                                        raise lab4()
                                    if not self.slice_from(u"U"):
                                        return False
                                    raise lab3()
                                except lab4: pass
                                self.cursor = v_4
                                if not self.eq_s(u"y"):
                                    raise lab2()
                                self.ket = self.cursor
                                if not self.in_grouping(GermanStemmer.g_v):
                                    raise lab2()
                                if not self.slice_from(u"Y"):
                                    return False
                            except lab3: pass
                            self.cursor = v_3
                            raise lab1()
                        except lab2: pass
                        self.cursor = v_3
                        if self.cursor >= self.limit:
                            raise lab0()
                        self.cursor += 1
                except lab1: pass
                continue
            except lab0: pass
            self.cursor = v_2
            break
        self.cursor = v_1
        while True:
            v_5 = self.cursor
            try:
                self.bra = self.cursor
                among_var = self.find_among(GermanStemmer.a_0)
                self.ket = self.cursor
                if among_var == 1:
                    if not self.slice_from(u"ss"):
                        return False
                elif among_var == 2:
                    if not self.slice_from(u"ä"):
                        return False
                elif among_var == 3:
                    if not self.slice_from(u"ö"):
                        return False
                elif among_var == 4:
                    if not self.slice_from(u"ü"):
                        return False
                elif among_var == 5:
                    if self.cursor >= self.limit:
                        raise lab5()
                    self.cursor += 1
                continue
            except lab5: pass
            self.cursor = v_5
            break
        return True

    def __r_mark_regions(self):
        self.I_p1 = self.limit
        self.I_p2 = self.limit
        v_1 = self.cursor
        c = self.cursor + 3
        if c > self.limit:
            return False
        self.cursor = c
        self.I_x = self.cursor
        self.cursor = v_1
        if not self.go_out_grouping(GermanStemmer.g_v):
            return False
        self.cursor += 1
        if not self.go_in_grouping(GermanStemmer.g_v):
            return False
        self.cursor += 1
        self.I_p1 = self.cursor
        try:
            if self.I_p1 >= self.I_x:
                raise lab0()
            self.I_p1 = self.I_x
        except lab0: pass
        if not self.go_out_grouping(GermanStemmer.g_v):
            return False
        self.cursor += 1
        if not self.go_in_grouping(GermanStemmer.g_v):
            return False
        self.cursor += 1
        self.I_p2 = self.cursor
        return True

    def __r_postlude(self):
        while True:
            v_1 = self.cursor
            try:
                self.bra = self.cursor
                among_var = self.find_among(GermanStemmer.a_1)
                self.ket = self.cursor
                if among_var == 1:
                    if not self.slice_from(u"y"):
                        return False
                elif among_var == 2:
                    if not self.slice_from(u"u"):
                        return False
                elif among_var == 3:
                    if not self.slice_from(u"a"):
                        return False
                elif among_var == 4:
                    if not self.slice_from(u"o"):
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

    def __r_R1(self):
        return self.I_p1 <= self.cursor

    def __r_R2(self):
        return self.I_p2 <= self.cursor

    def __r_standard_suffix(self):
        v_1 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            among_var = self.find_among_b(GermanStemmer.a_2)
            if among_var == 0:
                raise lab0()
            self.bra = self.cursor
            if not self.__r_R1():
                raise lab0()
            if among_var == 1:
                v_2 = self.limit - self.cursor
                try:
                    if not self.eq_s_b(u"syst"):
                        raise lab1()
                    raise lab0()
                except lab1: pass
                self.cursor = self.limit - v_2
                if not self.slice_del():
                    return False

            elif among_var == 2:
                if not self.slice_del():
                    return False

            elif among_var == 3:
                if not self.slice_del():
                    return False

                v_3 = self.limit - self.cursor
                try:
                    self.ket = self.cursor
                    if not self.eq_s_b(u"s"):
                        self.cursor = self.limit - v_3
                        raise lab2()
                    self.bra = self.cursor
                    if not self.eq_s_b(u"nis"):
                        self.cursor = self.limit - v_3
                        raise lab2()
                    if not self.slice_del():
                        return False

                except lab2: pass
            elif among_var == 4:
                if not self.in_grouping_b(GermanStemmer.g_s_ending):
                    raise lab0()
                if not self.slice_del():
                    return False

            else:
                if not self.slice_from(u"l"):
                    return False
        except lab0: pass
        self.cursor = self.limit - v_1
        v_4 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            among_var = self.find_among_b(GermanStemmer.a_4)
            if among_var == 0:
                raise lab3()
            self.bra = self.cursor
            if not self.__r_R1():
                raise lab3()
            if among_var == 1:
                if not self.slice_del():
                    return False

            elif among_var == 2:
                if not self.in_grouping_b(GermanStemmer.g_st_ending):
                    raise lab3()
                c = self.cursor - 3
                if c < self.limit_backward:
                    raise lab3()
                self.cursor = c
                if not self.slice_del():
                    return False

            else:
                v_5 = self.limit - self.cursor
                if not self.in_grouping_b(GermanStemmer.g_et_ending):
                    raise lab3()
                self.cursor = self.limit - v_5
                v_6 = self.limit - self.cursor
                try:
                    if self.find_among_b(GermanStemmer.a_3) == 0:
                        raise lab4()
                    raise lab3()
                except lab4: pass
                self.cursor = self.limit - v_6
                if not self.slice_del():
                    return False

        except lab3: pass
        self.cursor = self.limit - v_4
        v_7 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            among_var = self.find_among_b(GermanStemmer.a_6)
            if among_var == 0:
                raise lab5()
            self.bra = self.cursor
            if not self.__r_R2():
                raise lab5()
            if among_var == 1:
                if not self.slice_del():
                    return False

                v_8 = self.limit - self.cursor
                try:
                    self.ket = self.cursor
                    if not self.eq_s_b(u"ig"):
                        self.cursor = self.limit - v_8
                        raise lab6()
                    self.bra = self.cursor
                    v_9 = self.limit - self.cursor
                    try:
                        if not self.eq_s_b(u"e"):
                            raise lab7()
                        self.cursor = self.limit - v_8
                        raise lab6()
                    except lab7: pass
                    self.cursor = self.limit - v_9
                    if not self.__r_R2():
                        self.cursor = self.limit - v_8
                        raise lab6()
                    if not self.slice_del():
                        return False

                except lab6: pass
            elif among_var == 2:
                v_10 = self.limit - self.cursor
                try:
                    if not self.eq_s_b(u"e"):
                        raise lab8()
                    raise lab5()
                except lab8: pass
                self.cursor = self.limit - v_10
                if not self.slice_del():
                    return False

            elif among_var == 3:
                if not self.slice_del():
                    return False

                v_11 = self.limit - self.cursor
                try:
                    self.ket = self.cursor
                    try:
                        v_12 = self.limit - self.cursor
                        try:
                            if not self.eq_s_b(u"er"):
                                raise lab11()
                            raise lab10()
                        except lab11: pass
                        self.cursor = self.limit - v_12
                        if not self.eq_s_b(u"en"):
                            self.cursor = self.limit - v_11
                            raise lab9()
                    except lab10: pass
                    self.bra = self.cursor
                    if not self.__r_R1():
                        self.cursor = self.limit - v_11
                        raise lab9()
                    if not self.slice_del():
                        return False

                except lab9: pass
            else:
                if not self.slice_del():
                    return False

                v_13 = self.limit - self.cursor
                try:
                    self.ket = self.cursor
                    if self.find_among_b(GermanStemmer.a_5) == 0:
                        self.cursor = self.limit - v_13
                        raise lab12()
                    self.bra = self.cursor
                    if not self.__r_R2():
                        self.cursor = self.limit - v_13
                        raise lab12()
                    if not self.slice_del():
                        return False

                except lab12: pass
        except lab5: pass
        self.cursor = self.limit - v_7
        return True

    def _stem(self):
        v_1 = self.cursor
        self.__r_prelude()
        self.cursor = v_1
        v_2 = self.cursor
        self.__r_mark_regions()
        self.cursor = v_2
        self.limit_backward = self.cursor
        self.cursor = self.limit
        self.__r_standard_suffix()
        self.cursor = self.limit_backward
        v_4 = self.cursor
        self.__r_postlude()
        self.cursor = v_4
        return True

    a_0 = [
        Among(u"", -1, 5),
        Among(u"ae", 0, 2),
        Among(u"oe", 0, 3),
        Among(u"qu", 0, -1),
        Among(u"ue", 0, 4),
        Among(u"ß", 0, 1)
    ]

    a_1 = [
        Among(u"", -1, 5),
        Among(u"U", 0, 2),
        Among(u"Y", 0, 1),
        Among(u"ä", 0, 3),
        Among(u"ö", 0, 4),
        Among(u"ü", 0, 2)
    ]

    a_2 = [
        Among(u"e", -1, 3),
        Among(u"em", -1, 1),
        Among(u"en", -1, 3),
        Among(u"erinnen", 2, 2),
        Among(u"erin", -1, 2),
        Among(u"ln", -1, 5),
        Among(u"ern", -1, 2),
        Among(u"er", -1, 2),
        Among(u"s", -1, 4),
        Among(u"es", 8, 3),
        Among(u"lns", 8, 5)
    ]

    a_3 = [
        Among(u"tick", -1, -1),
        Among(u"plan", -1, -1),
        Among(u"geordn", -1, -1),
        Among(u"intern", -1, -1),
        Among(u"tr", -1, -1)
    ]

    a_4 = [
        Among(u"en", -1, 1),
        Among(u"er", -1, 1),
        Among(u"et", -1, 3),
        Among(u"st", -1, 2),
        Among(u"est", 3, 1)
    ]

    a_5 = [
        Among(u"ig", -1, 1),
        Among(u"lich", -1, 1)
    ]

    a_6 = [
        Among(u"end", -1, 1),
        Among(u"ig", -1, 2),
        Among(u"ung", -1, 1),
        Among(u"lich", -1, 3),
        Among(u"isch", -1, 2),
        Among(u"ik", -1, 2),
        Among(u"heit", -1, 3),
        Among(u"keit", -1, 4)
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

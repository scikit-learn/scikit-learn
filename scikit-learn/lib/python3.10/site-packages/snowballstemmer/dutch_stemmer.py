#-*- coding: utf-8 -*-
# Generated from dutch.sbl by Snowball 3.0.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class DutchStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from dutch.sbl by Snowball 3.0.1 - https://snowballstem.org/
    '''

    g_E = {u"e", u"ë", u"é", u"è", u"ê"}

    g_AIOU = {u"a", u"ä", u"á", u"à", u"â", u"i", u"ï", u"í", u"ì", u"î", u"o", u"ö", u"ó", u"ò", u"ô", u"u", u"ü", u"ú", u"ù", u"û"}

    g_AEIOU = {u"a", u"ä", u"á", u"à", u"â", u"e", u"ë", u"é", u"è", u"ê", u"i", u"ï", u"í", u"ì", u"î", u"o", u"ö", u"ó", u"ò", u"ô", u"u", u"ü", u"ú", u"ù", u"û"}

    g_v = {u"a", u"ä", u"á", u"à", u"â", u"e", u"ë", u"é", u"è", u"ê", u"i", u"ï", u"í", u"ì", u"î", u"o", u"ö", u"ó", u"ò", u"ô", u"u", u"ü", u"ú", u"ù", u"û", u"y"}

    g_v_WX = {u"a", u"ä", u"á", u"à", u"â", u"e", u"ë", u"é", u"è", u"ê", u"i", u"ï", u"í", u"ì", u"î", u"o", u"ö", u"ó", u"ò", u"ô", u"u", u"ü", u"ú", u"ù", u"û", u"y", u"w", u"x"}

    B_GE_removed = False
    B_stemmed = False
    I_p2 = 0
    I_p1 = 0
    S_ch = ""

    def __r_R1(self):
        return self.I_p1 <= self.cursor

    def __r_R2(self):
        return self.I_p2 <= self.cursor

    def __r_V(self):
        v_1 = self.limit - self.cursor
        try:
            v_2 = self.limit - self.cursor
            try:
                if not self.in_grouping_b(DutchStemmer.g_v):
                    raise lab1()
                raise lab0()
            except lab1: pass
            self.cursor = self.limit - v_2
            if not self.eq_s_b(u"ij"):
                return False
        except lab0: pass
        self.cursor = self.limit - v_1
        return True

    def __r_VX(self):
        v_1 = self.limit - self.cursor
        if self.cursor <= self.limit_backward:
            return False
        self.cursor -= 1
        try:
            v_2 = self.limit - self.cursor
            try:
                if not self.in_grouping_b(DutchStemmer.g_v):
                    raise lab1()
                raise lab0()
            except lab1: pass
            self.cursor = self.limit - v_2
            if not self.eq_s_b(u"ij"):
                return False
        except lab0: pass
        self.cursor = self.limit - v_1
        return True

    def __r_C(self):
        v_1 = self.limit - self.cursor
        v_2 = self.limit - self.cursor
        try:
            if not self.eq_s_b(u"ij"):
                raise lab0()
            return False
        except lab0: pass
        self.cursor = self.limit - v_2
        if not self.out_grouping_b(DutchStemmer.g_v):
            return False
        self.cursor = self.limit - v_1
        return True

    def __r_lengthen_V(self):
        v_1 = self.limit - self.cursor
        try:
            if not self.out_grouping_b(DutchStemmer.g_v_WX):
                raise lab0()
            self.ket = self.cursor
            among_var = self.find_among_b(DutchStemmer.a_0)
            if among_var == 0:
                raise lab0()
            self.bra = self.cursor
            if among_var == 1:
                v_2 = self.limit - self.cursor
                try:
                    v_3 = self.limit - self.cursor
                    try:
                        if not self.out_grouping_b(DutchStemmer.g_AEIOU):
                            raise lab2()
                        raise lab1()
                    except lab2: pass
                    self.cursor = self.limit - v_3
                    if self.cursor > self.limit_backward:
                        raise lab0()
                except lab1: pass
                self.cursor = self.limit - v_2
                self.S_ch = self.slice_to()
                if self.S_ch == '':
                    return False
                c = self.cursor
                self.insert(self.cursor, self.cursor, self.S_ch)
                self.cursor = c
            elif among_var == 2:
                v_4 = self.limit - self.cursor
                try:
                    v_5 = self.limit - self.cursor
                    try:
                        if not self.out_grouping_b(DutchStemmer.g_AEIOU):
                            raise lab4()
                        raise lab3()
                    except lab4: pass
                    self.cursor = self.limit - v_5
                    if self.cursor > self.limit_backward:
                        raise lab0()
                except lab3: pass
                v_6 = self.limit - self.cursor
                try:
                    try:
                        v_7 = self.limit - self.cursor
                        try:
                            if not self.in_grouping_b(DutchStemmer.g_AIOU):
                                raise lab7()
                            raise lab6()
                        except lab7: pass
                        self.cursor = self.limit - v_7
                        if not self.in_grouping_b(DutchStemmer.g_E):
                            raise lab5()
                        if self.cursor > self.limit_backward:
                            raise lab5()
                    except lab6: pass
                    raise lab0()
                except lab5: pass
                self.cursor = self.limit - v_6
                v_8 = self.limit - self.cursor
                try:
                    if self.cursor <= self.limit_backward:
                        raise lab8()
                    self.cursor -= 1
                    if not self.in_grouping_b(DutchStemmer.g_AIOU):
                        raise lab8()
                    if not self.out_grouping_b(DutchStemmer.g_AEIOU):
                        raise lab8()
                    raise lab0()
                except lab8: pass
                self.cursor = self.limit - v_8
                self.cursor = self.limit - v_4
                self.S_ch = self.slice_to()
                if self.S_ch == '':
                    return False
                c = self.cursor
                self.insert(self.cursor, self.cursor, self.S_ch)
                self.cursor = c
            elif among_var == 3:
                if not self.slice_from(u"eëe"):
                    return False
            else:
                if not self.slice_from(u"iee"):
                    return False
        except lab0: pass
        self.cursor = self.limit - v_1
        return True

    def __r_Step_1(self):
        self.ket = self.cursor
        among_var = self.find_among_b(DutchStemmer.a_1)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if not self.slice_del():
                return False

        elif among_var == 2:
            if not self.__r_R1():
                return False
            v_1 = self.limit - self.cursor
            try:
                if not self.eq_s_b(u"t"):
                    raise lab0()
                if not self.__r_R1():
                    raise lab0()
                return False
            except lab0: pass
            self.cursor = self.limit - v_1
            if not self.__r_C():
                return False
            if not self.slice_del():
                return False

        elif among_var == 3:
            if not self.__r_R1():
                return False
            if not self.slice_from(u"ie"):
                return False
        elif among_var == 4:
            try:
                v_2 = self.limit - self.cursor
                try:
                    v_3 = self.limit - self.cursor
                    if not self.eq_s_b(u"ar"):
                        raise lab2()
                    if not self.__r_R1():
                        raise lab2()
                    if not self.__r_C():
                        raise lab2()
                    self.cursor = self.limit - v_3
                    if not self.slice_del():
                        return False

                    self.__r_lengthen_V()
                    raise lab1()
                except lab2: pass
                self.cursor = self.limit - v_2
                try:
                    v_4 = self.limit - self.cursor
                    if not self.eq_s_b(u"er"):
                        raise lab3()
                    if not self.__r_R1():
                        raise lab3()
                    if not self.__r_C():
                        raise lab3()
                    self.cursor = self.limit - v_4
                    if not self.slice_del():
                        return False

                    raise lab1()
                except lab3: pass
                self.cursor = self.limit - v_2
                if not self.__r_R1():
                    return False
                if not self.__r_C():
                    return False
                if not self.slice_from(u"e"):
                    return False
            except lab1: pass
        elif among_var == 5:
            if not self.__r_R1():
                return False
            if not self.slice_from(u"é"):
                return False
        elif among_var == 6:
            if not self.__r_R1():
                return False
            if not self.__r_V():
                return False
            if not self.slice_from(u"au"):
                return False
        elif among_var == 7:
            try:
                v_5 = self.limit - self.cursor
                try:
                    if not self.eq_s_b(u"hed"):
                        raise lab5()
                    if not self.__r_R1():
                        raise lab5()
                    self.bra = self.cursor
                    if not self.slice_from(u"heid"):
                        return False
                    raise lab4()
                except lab5: pass
                self.cursor = self.limit - v_5
                try:
                    if not self.eq_s_b(u"nd"):
                        raise lab6()
                    if not self.slice_del():
                        return False

                    raise lab4()
                except lab6: pass
                self.cursor = self.limit - v_5
                try:
                    if not self.eq_s_b(u"d"):
                        raise lab7()
                    if not self.__r_R1():
                        raise lab7()
                    if not self.__r_C():
                        raise lab7()
                    self.bra = self.cursor
                    if not self.slice_del():
                        return False

                    raise lab4()
                except lab7: pass
                self.cursor = self.limit - v_5
                try:
                    try:
                        v_6 = self.limit - self.cursor
                        try:
                            if not self.eq_s_b(u"i"):
                                raise lab10()
                            raise lab9()
                        except lab10: pass
                        self.cursor = self.limit - v_6
                        if not self.eq_s_b(u"j"):
                            raise lab8()
                    except lab9: pass
                    if not self.__r_V():
                        raise lab8()
                    if not self.slice_del():
                        return False

                    raise lab4()
                except lab8: pass
                self.cursor = self.limit - v_5
                if not self.__r_R1():
                    return False
                if not self.__r_C():
                    return False
                if not self.slice_del():
                    return False

                self.__r_lengthen_V()
            except lab4: pass
        else:
            if not self.slice_from(u"nd"):
                return False
        return True

    def __r_Step_2(self):
        self.ket = self.cursor
        among_var = self.find_among_b(DutchStemmer.a_2)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            try:
                v_1 = self.limit - self.cursor
                try:
                    if not self.eq_s_b(u"'t"):
                        raise lab1()
                    self.bra = self.cursor
                    if not self.slice_del():
                        return False

                    raise lab0()
                except lab1: pass
                self.cursor = self.limit - v_1
                try:
                    if not self.eq_s_b(u"et"):
                        raise lab2()
                    self.bra = self.cursor
                    if not self.__r_R1():
                        raise lab2()
                    if not self.__r_C():
                        raise lab2()
                    if not self.slice_del():
                        return False

                    raise lab0()
                except lab2: pass
                self.cursor = self.limit - v_1
                try:
                    if not self.eq_s_b(u"rnt"):
                        raise lab3()
                    self.bra = self.cursor
                    if not self.slice_from(u"rn"):
                        return False
                    raise lab0()
                except lab3: pass
                self.cursor = self.limit - v_1
                try:
                    if not self.eq_s_b(u"t"):
                        raise lab4()
                    self.bra = self.cursor
                    if not self.__r_R1():
                        raise lab4()
                    if not self.__r_VX():
                        raise lab4()
                    if not self.slice_del():
                        return False

                    raise lab0()
                except lab4: pass
                self.cursor = self.limit - v_1
                try:
                    if not self.eq_s_b(u"ink"):
                        raise lab5()
                    self.bra = self.cursor
                    if not self.slice_from(u"ing"):
                        return False
                    raise lab0()
                except lab5: pass
                self.cursor = self.limit - v_1
                try:
                    if not self.eq_s_b(u"mp"):
                        raise lab6()
                    self.bra = self.cursor
                    if not self.slice_from(u"m"):
                        return False
                    raise lab0()
                except lab6: pass
                self.cursor = self.limit - v_1
                try:
                    if not self.eq_s_b(u"'"):
                        raise lab7()
                    self.bra = self.cursor
                    if not self.__r_R1():
                        raise lab7()
                    if not self.slice_del():
                        return False

                    raise lab0()
                except lab7: pass
                self.cursor = self.limit - v_1
                self.bra = self.cursor
                if not self.__r_R1():
                    return False
                if not self.__r_C():
                    return False
                if not self.slice_del():
                    return False

            except lab0: pass
        elif among_var == 2:
            if not self.__r_R1():
                return False
            if not self.slice_from(u"g"):
                return False
        elif among_var == 3:
            if not self.__r_R1():
                return False
            if not self.slice_from(u"lijk"):
                return False
        elif among_var == 4:
            if not self.__r_R1():
                return False
            if not self.slice_from(u"isch"):
                return False
        elif among_var == 5:
            if not self.__r_R1():
                return False
            if not self.__r_C():
                return False
            if not self.slice_del():
                return False

        elif among_var == 6:
            if not self.__r_R1():
                return False
            if not self.slice_from(u"t"):
                return False
        elif among_var == 7:
            if not self.__r_R1():
                return False
            if not self.slice_from(u"s"):
                return False
        elif among_var == 8:
            if not self.__r_R1():
                return False
            if not self.slice_from(u"r"):
                return False
        elif among_var == 9:
            if not self.__r_R1():
                return False
            if not self.slice_del():
                return False

            self.insert(self.cursor, self.cursor, u"l")
            self.__r_lengthen_V()
        elif among_var == 10:
            if not self.__r_R1():
                return False
            if not self.__r_C():
                return False
            if not self.slice_del():
                return False

            self.insert(self.cursor, self.cursor, u"en")
            self.__r_lengthen_V()
        else:
            if not self.__r_R1():
                return False
            if not self.__r_C():
                return False
            if not self.slice_from(u"ief"):
                return False
        return True

    def __r_Step_3(self):
        self.ket = self.cursor
        among_var = self.find_among_b(DutchStemmer.a_3)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if not self.__r_R1():
                return False
            if not self.slice_from(u"eer"):
                return False
        elif among_var == 2:
            if not self.__r_R1():
                return False
            if not self.slice_del():
                return False

            self.__r_lengthen_V()
        elif among_var == 3:
            if not self.__r_R1():
                return False
            if not self.slice_del():
                return False

        elif among_var == 4:
            if not self.slice_from(u"r"):
                return False
        elif among_var == 5:
            try:
                v_1 = self.limit - self.cursor
                try:
                    if not self.eq_s_b(u"ild"):
                        raise lab1()
                    if not self.slice_from(u"er"):
                        return False
                    raise lab0()
                except lab1: pass
                self.cursor = self.limit - v_1
                if not self.__r_R1():
                    return False
                if not self.slice_del():
                    return False

                self.__r_lengthen_V()
            except lab0: pass
        elif among_var == 6:
            if not self.__r_R1():
                return False
            if not self.__r_C():
                return False
            if not self.slice_from(u"aar"):
                return False
        elif among_var == 7:
            if not self.__r_R2():
                return False
            if not self.slice_del():
                return False

            self.insert(self.cursor, self.cursor, u"f")
            self.__r_lengthen_V()
        elif among_var == 8:
            if not self.__r_R2():
                return False
            if not self.slice_del():
                return False

            self.insert(self.cursor, self.cursor, u"g")
            self.__r_lengthen_V()
        elif among_var == 9:
            if not self.__r_R1():
                return False
            if not self.__r_C():
                return False
            if not self.slice_from(u"t"):
                return False
        else:
            if not self.__r_R1():
                return False
            if not self.__r_C():
                return False
            if not self.slice_from(u"d"):
                return False
        return True

    def __r_Step_4(self):
        try:
            v_1 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                among_var = self.find_among_b(DutchStemmer.a_4)
                if among_var == 0:
                    raise lab1()
                self.bra = self.cursor
                if among_var == 1:
                    if not self.__r_R1():
                        raise lab1()
                    if not self.slice_from(u"ie"):
                        return False
                elif among_var == 2:
                    if not self.__r_R1():
                        raise lab1()
                    if not self.slice_from(u"eer"):
                        return False
                elif among_var == 3:
                    if not self.__r_R1():
                        raise lab1()
                    if not self.slice_del():
                        return False

                elif among_var == 4:
                    if not self.__r_R1():
                        raise lab1()
                    if not self.__r_V():
                        raise lab1()
                    if not self.slice_from(u"n"):
                        return False
                elif among_var == 5:
                    if not self.__r_R1():
                        raise lab1()
                    if not self.__r_V():
                        raise lab1()
                    if not self.slice_from(u"l"):
                        return False
                elif among_var == 6:
                    if not self.__r_R1():
                        raise lab1()
                    if not self.__r_V():
                        raise lab1()
                    if not self.slice_from(u"r"):
                        return False
                elif among_var == 7:
                    if not self.__r_R1():
                        raise lab1()
                    if not self.slice_from(u"teer"):
                        return False
                elif among_var == 8:
                    if not self.__r_R1():
                        raise lab1()
                    if not self.slice_from(u"lijk"):
                        return False
                else:
                    if not self.__r_R1():
                        raise lab1()
                    if not self.__r_C():
                        raise lab1()
                    if not self.slice_del():
                        return False

                    self.__r_lengthen_V()
                raise lab0()
            except lab1: pass
            self.cursor = self.limit - v_1
            self.ket = self.cursor
            if self.find_among_b(DutchStemmer.a_5) == 0:
                return False
            self.bra = self.cursor
            if not self.__r_R1():
                return False
            v_2 = self.limit - self.cursor
            try:
                if not self.eq_s_b(u"inn"):
                    raise lab2()
                if self.cursor > self.limit_backward:
                    raise lab2()
                return False
            except lab2: pass
            self.cursor = self.limit - v_2
            if not self.__r_C():
                return False
            if not self.slice_del():
                return False

            self.__r_lengthen_V()
        except lab0: pass
        return True

    def __r_Step_7(self):
        self.ket = self.cursor
        among_var = self.find_among_b(DutchStemmer.a_6)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if not self.slice_from(u"k"):
                return False
        elif among_var == 2:
            if not self.slice_from(u"f"):
                return False
        else:
            if not self.slice_from(u"p"):
                return False
        return True

    def __r_Step_6(self):
        self.ket = self.cursor
        among_var = self.find_among_b(DutchStemmer.a_7)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if not self.slice_from(u"b"):
                return False
        elif among_var == 2:
            if not self.slice_from(u"c"):
                return False
        elif among_var == 3:
            if not self.slice_from(u"d"):
                return False
        elif among_var == 4:
            if not self.slice_from(u"f"):
                return False
        elif among_var == 5:
            if not self.slice_from(u"g"):
                return False
        elif among_var == 6:
            if not self.slice_from(u"h"):
                return False
        elif among_var == 7:
            if not self.slice_from(u"j"):
                return False
        elif among_var == 8:
            if not self.slice_from(u"k"):
                return False
        elif among_var == 9:
            if not self.slice_from(u"l"):
                return False
        elif among_var == 10:
            if not self.slice_from(u"m"):
                return False
        elif among_var == 11:
            v_1 = self.limit - self.cursor
            try:
                if not self.eq_s_b(u"i"):
                    raise lab0()
                if self.cursor > self.limit_backward:
                    raise lab0()
                return False
            except lab0: pass
            self.cursor = self.limit - v_1
            if not self.slice_from(u"n"):
                return False
        elif among_var == 12:
            if not self.slice_from(u"p"):
                return False
        elif among_var == 13:
            if not self.slice_from(u"q"):
                return False
        elif among_var == 14:
            if not self.slice_from(u"r"):
                return False
        elif among_var == 15:
            if not self.slice_from(u"s"):
                return False
        elif among_var == 16:
            if not self.slice_from(u"t"):
                return False
        elif among_var == 17:
            if not self.slice_from(u"v"):
                return False
        elif among_var == 18:
            if not self.slice_from(u"w"):
                return False
        elif among_var == 19:
            if not self.slice_from(u"x"):
                return False
        else:
            if not self.slice_from(u"z"):
                return False
        return True

    def __r_Step_1c(self):
        self.ket = self.cursor
        among_var = self.find_among_b(DutchStemmer.a_8)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R1():
            return False
        if not self.__r_C():
            return False
        if among_var == 1:
            v_1 = self.limit - self.cursor
            try:
                if not self.eq_s_b(u"n"):
                    raise lab0()
                if not self.__r_R1():
                    raise lab0()
                return False
            except lab0: pass
            self.cursor = self.limit - v_1
            try:
                v_2 = self.limit - self.cursor
                try:
                    if not self.eq_s_b(u"in"):
                        raise lab2()
                    if self.cursor > self.limit_backward:
                        raise lab2()
                    if not self.slice_from(u"n"):
                        return False
                    raise lab1()
                except lab2: pass
                self.cursor = self.limit - v_2
                if not self.slice_del():
                    return False

            except lab1: pass
        else:
            v_3 = self.limit - self.cursor
            try:
                if not self.eq_s_b(u"h"):
                    raise lab3()
                if not self.__r_R1():
                    raise lab3()
                return False
            except lab3: pass
            self.cursor = self.limit - v_3
            v_4 = self.limit - self.cursor
            try:
                if not self.eq_s_b(u"en"):
                    raise lab4()
                if self.cursor > self.limit_backward:
                    raise lab4()
                return False
            except lab4: pass
            self.cursor = self.limit - v_4
            if not self.slice_del():
                return False

        return True

    def __r_Lose_prefix(self):
        self.bra = self.cursor
        if not self.eq_s(u"ge"):
            return False
        self.ket = self.cursor
        v_1 = self.cursor
        c = self.cursor + 3
        if c > self.limit:
            return False
        self.cursor = c
        self.cursor = v_1
        v_2 = self.cursor
        try:
            while True:
                v_3 = self.cursor
                try:
                    try:
                        v_4 = self.cursor
                        try:
                            if not self.eq_s(u"ij"):
                                raise lab3()
                            raise lab2()
                        except lab3: pass
                        self.cursor = v_4
                        if not self.in_grouping(DutchStemmer.g_v):
                            raise lab1()
                    except lab2: pass
                    raise lab0()
                except lab1: pass
                self.cursor = v_3
                if self.cursor >= self.limit:
                    return False
                self.cursor += 1
        except lab0: pass
        while True:
            v_5 = self.cursor
            try:
                try:
                    v_6 = self.cursor
                    try:
                        if not self.eq_s(u"ij"):
                            raise lab6()
                        raise lab5()
                    except lab6: pass
                    self.cursor = v_6
                    if not self.in_grouping(DutchStemmer.g_v):
                        raise lab4()
                except lab5: pass
                continue
            except lab4: pass
            self.cursor = v_5
            break
        try:
            if self.cursor < self.limit:
                raise lab7()
            return False
        except lab7: pass
        self.cursor = v_2
        among_var = self.find_among(DutchStemmer.a_9)
        if among_var == 1:
            return False
        self.B_GE_removed = True
        if not self.slice_del():
            return False

        v_8 = self.cursor
        try:
            self.bra = self.cursor
            among_var = self.find_among(DutchStemmer.a_10)
            if among_var == 0:
                raise lab8()
            self.ket = self.cursor
            if among_var == 1:
                if not self.slice_from(u"e"):
                    return False
            else:
                if not self.slice_from(u"i"):
                    return False
        except lab8: pass
        self.cursor = v_8
        return True

    def __r_Lose_infix(self):
        if self.cursor >= self.limit:
            return False
        self.cursor += 1
        try:
            while True:
                try:
                    self.bra = self.cursor
                    if not self.eq_s(u"ge"):
                        raise lab1()
                    self.ket = self.cursor
                    raise lab0()
                except lab1: pass
                if self.cursor >= self.limit:
                    return False
                self.cursor += 1
        except lab0: pass
        v_2 = self.cursor
        c = self.cursor + 3
        if c > self.limit:
            return False
        self.cursor = c
        self.cursor = v_2
        v_3 = self.cursor
        try:
            while True:
                v_4 = self.cursor
                try:
                    try:
                        v_5 = self.cursor
                        try:
                            if not self.eq_s(u"ij"):
                                raise lab5()
                            raise lab4()
                        except lab5: pass
                        self.cursor = v_5
                        if not self.in_grouping(DutchStemmer.g_v):
                            raise lab3()
                    except lab4: pass
                    raise lab2()
                except lab3: pass
                self.cursor = v_4
                if self.cursor >= self.limit:
                    return False
                self.cursor += 1
        except lab2: pass
        while True:
            v_6 = self.cursor
            try:
                try:
                    v_7 = self.cursor
                    try:
                        if not self.eq_s(u"ij"):
                            raise lab8()
                        raise lab7()
                    except lab8: pass
                    self.cursor = v_7
                    if not self.in_grouping(DutchStemmer.g_v):
                        raise lab6()
                except lab7: pass
                continue
            except lab6: pass
            self.cursor = v_6
            break
        try:
            if self.cursor < self.limit:
                raise lab9()
            return False
        except lab9: pass
        self.cursor = v_3
        self.B_GE_removed = True
        if not self.slice_del():
            return False

        v_9 = self.cursor
        try:
            self.bra = self.cursor
            among_var = self.find_among(DutchStemmer.a_11)
            if among_var == 0:
                raise lab10()
            self.ket = self.cursor
            if among_var == 1:
                if not self.slice_from(u"e"):
                    return False
            else:
                if not self.slice_from(u"i"):
                    return False
        except lab10: pass
        self.cursor = v_9
        return True

    def __r_measure(self):
        self.I_p1 = self.limit
        self.I_p2 = self.limit
        v_1 = self.cursor
        try:
            while True:
                try:
                    if not self.out_grouping(DutchStemmer.g_v):
                        raise lab1()
                    continue
                except lab1: pass
                break
            v_3 = 1
            while True:
                v_4 = self.cursor
                try:
                    try:
                        v_5 = self.cursor
                        try:
                            if not self.eq_s(u"ij"):
                                raise lab4()
                            raise lab3()
                        except lab4: pass
                        self.cursor = v_5
                        if not self.in_grouping(DutchStemmer.g_v):
                            raise lab2()
                    except lab3: pass
                    v_3 -= 1
                    continue
                except lab2: pass
                self.cursor = v_4
                break
            if v_3 > 0:
                raise lab0()
            if not self.out_grouping(DutchStemmer.g_v):
                raise lab0()
            self.I_p1 = self.cursor
            while True:
                try:
                    if not self.out_grouping(DutchStemmer.g_v):
                        raise lab5()
                    continue
                except lab5: pass
                break
            v_7 = 1
            while True:
                v_8 = self.cursor
                try:
                    try:
                        v_9 = self.cursor
                        try:
                            if not self.eq_s(u"ij"):
                                raise lab8()
                            raise lab7()
                        except lab8: pass
                        self.cursor = v_9
                        if not self.in_grouping(DutchStemmer.g_v):
                            raise lab6()
                    except lab7: pass
                    v_7 -= 1
                    continue
                except lab6: pass
                self.cursor = v_8
                break
            if v_7 > 0:
                raise lab0()
            if not self.out_grouping(DutchStemmer.g_v):
                raise lab0()
            self.I_p2 = self.cursor
        except lab0: pass
        self.cursor = v_1
        return True

    def _stem(self):
        self.B_stemmed = False
        self.__r_measure()
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_1 = self.limit - self.cursor
        try:
            if not self.__r_Step_1():
                raise lab0()
            self.B_stemmed = True
        except lab0: pass
        self.cursor = self.limit - v_1
        v_2 = self.limit - self.cursor
        try:
            if not self.__r_Step_2():
                raise lab1()
            self.B_stemmed = True
        except lab1: pass
        self.cursor = self.limit - v_2
        v_3 = self.limit - self.cursor
        try:
            if not self.__r_Step_3():
                raise lab2()
            self.B_stemmed = True
        except lab2: pass
        self.cursor = self.limit - v_3
        v_4 = self.limit - self.cursor
        try:
            if not self.__r_Step_4():
                raise lab3()
            self.B_stemmed = True
        except lab3: pass
        self.cursor = self.limit - v_4
        self.cursor = self.limit_backward
        self.B_GE_removed = False
        v_5 = self.cursor
        try:
            v_6 = self.cursor
            if not self.__r_Lose_prefix():
                raise lab4()
            self.cursor = v_6
            self.__r_measure()
        except lab4: pass
        self.cursor = v_5
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_7 = self.limit - self.cursor
        try:
            if not self.B_GE_removed:
                raise lab5()
            self.B_stemmed = True
            if not self.__r_Step_1c():
                raise lab5()
        except lab5: pass
        self.cursor = self.limit - v_7
        self.cursor = self.limit_backward
        self.B_GE_removed = False
        v_8 = self.cursor
        try:
            v_9 = self.cursor
            if not self.__r_Lose_infix():
                raise lab6()
            self.cursor = v_9
            self.__r_measure()
        except lab6: pass
        self.cursor = v_8
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_10 = self.limit - self.cursor
        try:
            if not self.B_GE_removed:
                raise lab7()
            self.B_stemmed = True
            if not self.__r_Step_1c():
                raise lab7()
        except lab7: pass
        self.cursor = self.limit - v_10
        self.cursor = self.limit_backward
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_11 = self.limit - self.cursor
        try:
            if not self.__r_Step_7():
                raise lab8()
            self.B_stemmed = True
        except lab8: pass
        self.cursor = self.limit - v_11
        v_12 = self.limit - self.cursor
        try:
            if not self.B_stemmed:
                raise lab9()
            if not self.__r_Step_6():
                raise lab9()
        except lab9: pass
        self.cursor = self.limit - v_12
        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among(u"a", -1, 1),
        Among(u"e", -1, 2),
        Among(u"o", -1, 1),
        Among(u"u", -1, 1),
        Among(u"à", -1, 1),
        Among(u"á", -1, 1),
        Among(u"â", -1, 1),
        Among(u"ä", -1, 1),
        Among(u"è", -1, 2),
        Among(u"é", -1, 2),
        Among(u"ê", -1, 2),
        Among(u"eë", -1, 3),
        Among(u"ië", -1, 4),
        Among(u"ò", -1, 1),
        Among(u"ó", -1, 1),
        Among(u"ô", -1, 1),
        Among(u"ö", -1, 1),
        Among(u"ù", -1, 1),
        Among(u"ú", -1, 1),
        Among(u"û", -1, 1),
        Among(u"ü", -1, 1)
    ]

    a_1 = [
        Among(u"nde", -1, 8),
        Among(u"en", -1, 7),
        Among(u"s", -1, 2),
        Among(u"'s", 2, 1),
        Among(u"es", 2, 4),
        Among(u"ies", 4, 3),
        Among(u"aus", 2, 6),
        Among(u"és", 2, 5)
    ]

    a_2 = [
        Among(u"de", -1, 5),
        Among(u"ge", -1, 2),
        Among(u"ische", -1, 4),
        Among(u"je", -1, 1),
        Among(u"lijke", -1, 3),
        Among(u"le", -1, 9),
        Among(u"ene", -1, 10),
        Among(u"re", -1, 8),
        Among(u"se", -1, 7),
        Among(u"te", -1, 6),
        Among(u"ieve", -1, 11)
    ]

    a_3 = [
        Among(u"heid", -1, 3),
        Among(u"fie", -1, 7),
        Among(u"gie", -1, 8),
        Among(u"atie", -1, 1),
        Among(u"isme", -1, 5),
        Among(u"ing", -1, 5),
        Among(u"arij", -1, 6),
        Among(u"erij", -1, 5),
        Among(u"sel", -1, 3),
        Among(u"rder", -1, 4),
        Among(u"ster", -1, 3),
        Among(u"iteit", -1, 2),
        Among(u"dst", -1, 10),
        Among(u"tst", -1, 9)
    ]

    a_4 = [
        Among(u"end", -1, 9),
        Among(u"atief", -1, 2),
        Among(u"erig", -1, 9),
        Among(u"achtig", -1, 3),
        Among(u"ioneel", -1, 1),
        Among(u"baar", -1, 3),
        Among(u"laar", -1, 5),
        Among(u"naar", -1, 4),
        Among(u"raar", -1, 6),
        Among(u"eriger", -1, 9),
        Among(u"achtiger", -1, 3),
        Among(u"lijker", -1, 8),
        Among(u"tant", -1, 7),
        Among(u"erigst", -1, 9),
        Among(u"achtigst", -1, 3),
        Among(u"lijkst", -1, 8)
    ]

    a_5 = [
        Among(u"ig", -1, 1),
        Among(u"iger", -1, 1),
        Among(u"igst", -1, 1)
    ]

    a_6 = [
        Among(u"ft", -1, 2),
        Among(u"kt", -1, 1),
        Among(u"pt", -1, 3)
    ]

    a_7 = [
        Among(u"bb", -1, 1),
        Among(u"cc", -1, 2),
        Among(u"dd", -1, 3),
        Among(u"ff", -1, 4),
        Among(u"gg", -1, 5),
        Among(u"hh", -1, 6),
        Among(u"jj", -1, 7),
        Among(u"kk", -1, 8),
        Among(u"ll", -1, 9),
        Among(u"mm", -1, 10),
        Among(u"nn", -1, 11),
        Among(u"pp", -1, 12),
        Among(u"qq", -1, 13),
        Among(u"rr", -1, 14),
        Among(u"ss", -1, 15),
        Among(u"tt", -1, 16),
        Among(u"v", -1, 4),
        Among(u"vv", 16, 17),
        Among(u"ww", -1, 18),
        Among(u"xx", -1, 19),
        Among(u"z", -1, 15),
        Among(u"zz", 20, 20)
    ]

    a_8 = [
        Among(u"d", -1, 1),
        Among(u"t", -1, 2)
    ]

    a_9 = [
        Among(u"", -1, -1),
        Among(u"eft", 0, 1),
        Among(u"vaa", 0, 1),
        Among(u"val", 0, 1),
        Among(u"vali", 3, -1),
        Among(u"vare", 0, 1)
    ]

    a_10 = [
        Among(u"ë", -1, 1),
        Among(u"ï", -1, 2)
    ]

    a_11 = [
        Among(u"ë", -1, 1),
        Among(u"ï", -1, 2)
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

#-*- coding: utf-8 -*-
# Generated from russian.sbl by Snowball 3.0.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class RussianStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from russian.sbl by Snowball 3.0.1 - https://snowballstem.org/
    '''

    g_v = {u"а", u"е", u"и", u"о", u"у", u"ы", u"э", u"ю", u"я"}

    I_p2 = 0
    I_pV = 0

    def __r_mark_regions(self):
        self.I_pV = self.limit
        self.I_p2 = self.limit
        v_1 = self.cursor
        try:
            if not self.go_out_grouping(RussianStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_pV = self.cursor
            if not self.go_in_grouping(RussianStemmer.g_v):
                raise lab0()
            self.cursor += 1
            if not self.go_out_grouping(RussianStemmer.g_v):
                raise lab0()
            self.cursor += 1
            if not self.go_in_grouping(RussianStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_p2 = self.cursor
        except lab0: pass
        self.cursor = v_1
        return True

    def __r_R2(self):
        return self.I_p2 <= self.cursor

    def __r_perfective_gerund(self):
        self.ket = self.cursor
        among_var = self.find_among_b(RussianStemmer.a_0)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            try:
                v_1 = self.limit - self.cursor
                try:
                    if not self.eq_s_b(u"а"):
                        raise lab1()
                    raise lab0()
                except lab1: pass
                self.cursor = self.limit - v_1
                if not self.eq_s_b(u"я"):
                    return False
            except lab0: pass
            if not self.slice_del():
                return False

        else:
            if not self.slice_del():
                return False

        return True

    def __r_adjective(self):
        self.ket = self.cursor
        if self.find_among_b(RussianStemmer.a_1) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        return True

    def __r_adjectival(self):
        if not self.__r_adjective():
            return False
        v_1 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            among_var = self.find_among_b(RussianStemmer.a_2)
            if among_var == 0:
                self.cursor = self.limit - v_1
                raise lab0()
            self.bra = self.cursor
            if among_var == 1:
                try:
                    v_2 = self.limit - self.cursor
                    try:
                        if not self.eq_s_b(u"а"):
                            raise lab2()
                        raise lab1()
                    except lab2: pass
                    self.cursor = self.limit - v_2
                    if not self.eq_s_b(u"я"):
                        self.cursor = self.limit - v_1
                        raise lab0()
                except lab1: pass
                if not self.slice_del():
                    return False

            else:
                if not self.slice_del():
                    return False

        except lab0: pass
        return True

    def __r_reflexive(self):
        self.ket = self.cursor
        if self.find_among_b(RussianStemmer.a_3) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        return True

    def __r_verb(self):
        self.ket = self.cursor
        among_var = self.find_among_b(RussianStemmer.a_4)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            try:
                v_1 = self.limit - self.cursor
                try:
                    if not self.eq_s_b(u"а"):
                        raise lab1()
                    raise lab0()
                except lab1: pass
                self.cursor = self.limit - v_1
                if not self.eq_s_b(u"я"):
                    return False
            except lab0: pass
            if not self.slice_del():
                return False

        else:
            if not self.slice_del():
                return False

        return True

    def __r_noun(self):
        self.ket = self.cursor
        if self.find_among_b(RussianStemmer.a_5) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        return True

    def __r_derivational(self):
        self.ket = self.cursor
        if self.find_among_b(RussianStemmer.a_6) == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R2():
            return False
        if not self.slice_del():
            return False

        return True

    def __r_tidy_up(self):
        self.ket = self.cursor
        among_var = self.find_among_b(RussianStemmer.a_7)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if not self.slice_del():
                return False

            self.ket = self.cursor
            if not self.eq_s_b(u"н"):
                return False
            self.bra = self.cursor
            if not self.eq_s_b(u"н"):
                return False
            if not self.slice_del():
                return False

        elif among_var == 2:
            if not self.eq_s_b(u"н"):
                return False
            if not self.slice_del():
                return False

        else:
            if not self.slice_del():
                return False

        return True

    def _stem(self):
        v_1 = self.cursor
        try:
            while True:
                v_2 = self.cursor
                try:
                    try:
                        while True:
                            v_3 = self.cursor
                            try:
                                self.bra = self.cursor
                                if not self.eq_s(u"ё"):
                                    raise lab3()
                                self.ket = self.cursor
                                self.cursor = v_3
                                raise lab2()
                            except lab3: pass
                            self.cursor = v_3
                            if self.cursor >= self.limit:
                                raise lab1()
                            self.cursor += 1
                    except lab2: pass
                    if not self.slice_from(u"е"):
                        return False
                    continue
                except lab1: pass
                self.cursor = v_2
                break
        except lab0: pass
        self.cursor = v_1
        self.__r_mark_regions()
        self.limit_backward = self.cursor
        self.cursor = self.limit
        if self.cursor < self.I_pV:
            return False
        v_6 = self.limit_backward
        self.limit_backward = self.I_pV
        v_7 = self.limit - self.cursor
        try:
            try:
                v_8 = self.limit - self.cursor
                try:
                    if not self.__r_perfective_gerund():
                        raise lab6()
                    raise lab5()
                except lab6: pass
                self.cursor = self.limit - v_8
                v_9 = self.limit - self.cursor
                try:
                    if not self.__r_reflexive():
                        self.cursor = self.limit - v_9
                        raise lab7()
                except lab7: pass
                try:
                    v_10 = self.limit - self.cursor
                    try:
                        if not self.__r_adjectival():
                            raise lab9()
                        raise lab8()
                    except lab9: pass
                    self.cursor = self.limit - v_10
                    try:
                        if not self.__r_verb():
                            raise lab10()
                        raise lab8()
                    except lab10: pass
                    self.cursor = self.limit - v_10
                    if not self.__r_noun():
                        raise lab4()
                except lab8: pass
            except lab5: pass
        except lab4: pass
        self.cursor = self.limit - v_7
        v_11 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if not self.eq_s_b(u"и"):
                self.cursor = self.limit - v_11
                raise lab11()
            self.bra = self.cursor
            if not self.slice_del():
                return False

        except lab11: pass
        v_12 = self.limit - self.cursor
        self.__r_derivational()
        self.cursor = self.limit - v_12
        v_13 = self.limit - self.cursor
        self.__r_tidy_up()
        self.cursor = self.limit - v_13
        self.limit_backward = v_6
        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among(u"в", -1, 1),
        Among(u"ив", 0, 2),
        Among(u"ыв", 0, 2),
        Among(u"вши", -1, 1),
        Among(u"ивши", 3, 2),
        Among(u"ывши", 3, 2),
        Among(u"вшись", -1, 1),
        Among(u"ившись", 6, 2),
        Among(u"ывшись", 6, 2)
    ]

    a_1 = [
        Among(u"ее", -1, 1),
        Among(u"ие", -1, 1),
        Among(u"ое", -1, 1),
        Among(u"ые", -1, 1),
        Among(u"ими", -1, 1),
        Among(u"ыми", -1, 1),
        Among(u"ей", -1, 1),
        Among(u"ий", -1, 1),
        Among(u"ой", -1, 1),
        Among(u"ый", -1, 1),
        Among(u"ем", -1, 1),
        Among(u"им", -1, 1),
        Among(u"ом", -1, 1),
        Among(u"ым", -1, 1),
        Among(u"его", -1, 1),
        Among(u"ого", -1, 1),
        Among(u"ему", -1, 1),
        Among(u"ому", -1, 1),
        Among(u"их", -1, 1),
        Among(u"ых", -1, 1),
        Among(u"ею", -1, 1),
        Among(u"ою", -1, 1),
        Among(u"ую", -1, 1),
        Among(u"юю", -1, 1),
        Among(u"ая", -1, 1),
        Among(u"яя", -1, 1)
    ]

    a_2 = [
        Among(u"ем", -1, 1),
        Among(u"нн", -1, 1),
        Among(u"вш", -1, 1),
        Among(u"ивш", 2, 2),
        Among(u"ывш", 2, 2),
        Among(u"щ", -1, 1),
        Among(u"ющ", 5, 1),
        Among(u"ующ", 6, 2)
    ]

    a_3 = [
        Among(u"сь", -1, 1),
        Among(u"ся", -1, 1)
    ]

    a_4 = [
        Among(u"ла", -1, 1),
        Among(u"ила", 0, 2),
        Among(u"ыла", 0, 2),
        Among(u"на", -1, 1),
        Among(u"ена", 3, 2),
        Among(u"ете", -1, 1),
        Among(u"ите", -1, 2),
        Among(u"йте", -1, 1),
        Among(u"ейте", 7, 2),
        Among(u"уйте", 7, 2),
        Among(u"ли", -1, 1),
        Among(u"или", 10, 2),
        Among(u"ыли", 10, 2),
        Among(u"й", -1, 1),
        Among(u"ей", 13, 2),
        Among(u"уй", 13, 2),
        Among(u"л", -1, 1),
        Among(u"ил", 16, 2),
        Among(u"ыл", 16, 2),
        Among(u"ем", -1, 1),
        Among(u"им", -1, 2),
        Among(u"ым", -1, 2),
        Among(u"н", -1, 1),
        Among(u"ен", 22, 2),
        Among(u"ло", -1, 1),
        Among(u"ило", 24, 2),
        Among(u"ыло", 24, 2),
        Among(u"но", -1, 1),
        Among(u"ено", 27, 2),
        Among(u"нно", 27, 1),
        Among(u"ет", -1, 1),
        Among(u"ует", 30, 2),
        Among(u"ит", -1, 2),
        Among(u"ыт", -1, 2),
        Among(u"ют", -1, 1),
        Among(u"уют", 34, 2),
        Among(u"ят", -1, 2),
        Among(u"ны", -1, 1),
        Among(u"ены", 37, 2),
        Among(u"ть", -1, 1),
        Among(u"ить", 39, 2),
        Among(u"ыть", 39, 2),
        Among(u"ешь", -1, 1),
        Among(u"ишь", -1, 2),
        Among(u"ю", -1, 2),
        Among(u"ую", 44, 2)
    ]

    a_5 = [
        Among(u"а", -1, 1),
        Among(u"ев", -1, 1),
        Among(u"ов", -1, 1),
        Among(u"е", -1, 1),
        Among(u"ие", 3, 1),
        Among(u"ье", 3, 1),
        Among(u"и", -1, 1),
        Among(u"еи", 6, 1),
        Among(u"ии", 6, 1),
        Among(u"ами", 6, 1),
        Among(u"ями", 6, 1),
        Among(u"иями", 10, 1),
        Among(u"й", -1, 1),
        Among(u"ей", 12, 1),
        Among(u"ией", 13, 1),
        Among(u"ий", 12, 1),
        Among(u"ой", 12, 1),
        Among(u"ам", -1, 1),
        Among(u"ем", -1, 1),
        Among(u"ием", 18, 1),
        Among(u"ом", -1, 1),
        Among(u"ям", -1, 1),
        Among(u"иям", 21, 1),
        Among(u"о", -1, 1),
        Among(u"у", -1, 1),
        Among(u"ах", -1, 1),
        Among(u"ях", -1, 1),
        Among(u"иях", 26, 1),
        Among(u"ы", -1, 1),
        Among(u"ь", -1, 1),
        Among(u"ю", -1, 1),
        Among(u"ию", 30, 1),
        Among(u"ью", 30, 1),
        Among(u"я", -1, 1),
        Among(u"ия", 33, 1),
        Among(u"ья", 33, 1)
    ]

    a_6 = [
        Among(u"ост", -1, 1),
        Among(u"ость", -1, 1)
    ]

    a_7 = [
        Among(u"ейше", -1, 1),
        Among(u"н", -1, 2),
        Among(u"ейш", -1, 1),
        Among(u"ь", -1, 3)
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

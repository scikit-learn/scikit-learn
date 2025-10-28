#-*- coding: utf-8 -*-
# Generated from romanian.sbl by Snowball 3.0.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class RomanianStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from romanian.sbl by Snowball 3.0.1 - https://snowballstem.org/
    '''

    g_v = {u"a", u"e", u"i", u"o", u"u", u"â", u"î", u"ă"}

    B_standard_suffix_removed = False
    I_p2 = 0
    I_p1 = 0
    I_pV = 0

    def __r_norm(self):
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
                                among_var = self.find_among(RomanianStemmer.a_0)
                                if among_var == 0:
                                    raise lab3()
                                self.ket = self.cursor
                                if among_var == 1:
                                    if not self.slice_from(u"ș"):
                                        return False
                                else:
                                    if not self.slice_from(u"ț"):
                                        return False
                                self.cursor = v_3
                                raise lab2()
                            except lab3: pass
                            self.cursor = v_3
                            if self.cursor >= self.limit:
                                raise lab1()
                            self.cursor += 1
                    except lab2: pass
                    continue
                except lab1: pass
                self.cursor = v_2
                break
        except lab0: pass
        self.cursor = v_1
        return True

    def __r_prelude(self):
        while True:
            v_1 = self.cursor
            try:
                try:
                    while True:
                        v_2 = self.cursor
                        try:
                            if not self.in_grouping(RomanianStemmer.g_v):
                                raise lab2()
                            self.bra = self.cursor
                            try:
                                v_3 = self.cursor
                                try:
                                    if not self.eq_s(u"u"):
                                        raise lab4()
                                    self.ket = self.cursor
                                    if not self.in_grouping(RomanianStemmer.g_v):
                                        raise lab4()
                                    if not self.slice_from(u"U"):
                                        return False
                                    raise lab3()
                                except lab4: pass
                                self.cursor = v_3
                                if not self.eq_s(u"i"):
                                    raise lab2()
                                self.ket = self.cursor
                                if not self.in_grouping(RomanianStemmer.g_v):
                                    raise lab2()
                                if not self.slice_from(u"I"):
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
                    if not self.in_grouping(RomanianStemmer.g_v):
                        raise lab2()
                    try:
                        v_3 = self.cursor
                        try:
                            if not self.out_grouping(RomanianStemmer.g_v):
                                raise lab4()
                            if not self.go_out_grouping(RomanianStemmer.g_v):
                                raise lab4()
                            self.cursor += 1
                            raise lab3()
                        except lab4: pass
                        self.cursor = v_3
                        if not self.in_grouping(RomanianStemmer.g_v):
                            raise lab2()
                        if not self.go_in_grouping(RomanianStemmer.g_v):
                            raise lab2()
                        self.cursor += 1
                    except lab3: pass
                    raise lab1()
                except lab2: pass
                self.cursor = v_2
                if not self.out_grouping(RomanianStemmer.g_v):
                    raise lab0()
                try:
                    v_4 = self.cursor
                    try:
                        if not self.out_grouping(RomanianStemmer.g_v):
                            raise lab6()
                        if not self.go_out_grouping(RomanianStemmer.g_v):
                            raise lab6()
                        self.cursor += 1
                        raise lab5()
                    except lab6: pass
                    self.cursor = v_4
                    if not self.in_grouping(RomanianStemmer.g_v):
                        raise lab0()
                    if self.cursor >= self.limit:
                        raise lab0()
                    self.cursor += 1
                except lab5: pass
            except lab1: pass
            self.I_pV = self.cursor
        except lab0: pass
        self.cursor = v_1
        v_5 = self.cursor
        try:
            if not self.go_out_grouping(RomanianStemmer.g_v):
                raise lab7()
            self.cursor += 1
            if not self.go_in_grouping(RomanianStemmer.g_v):
                raise lab7()
            self.cursor += 1
            self.I_p1 = self.cursor
            if not self.go_out_grouping(RomanianStemmer.g_v):
                raise lab7()
            self.cursor += 1
            if not self.go_in_grouping(RomanianStemmer.g_v):
                raise lab7()
            self.cursor += 1
            self.I_p2 = self.cursor
        except lab7: pass
        self.cursor = v_5
        return True

    def __r_postlude(self):
        while True:
            v_1 = self.cursor
            try:
                self.bra = self.cursor
                among_var = self.find_among(RomanianStemmer.a_1)
                self.ket = self.cursor
                if among_var == 1:
                    if not self.slice_from(u"i"):
                        return False
                elif among_var == 2:
                    if not self.slice_from(u"u"):
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

    def __r_step_0(self):
        self.ket = self.cursor
        among_var = self.find_among_b(RomanianStemmer.a_2)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R1():
            return False
        if among_var == 1:
            if not self.slice_del():
                return False

        elif among_var == 2:
            if not self.slice_from(u"a"):
                return False
        elif among_var == 3:
            if not self.slice_from(u"e"):
                return False
        elif among_var == 4:
            if not self.slice_from(u"i"):
                return False
        elif among_var == 5:
            v_1 = self.limit - self.cursor
            try:
                if not self.eq_s_b(u"ab"):
                    raise lab0()
                return False
            except lab0: pass
            self.cursor = self.limit - v_1
            if not self.slice_from(u"i"):
                return False
        elif among_var == 6:
            if not self.slice_from(u"at"):
                return False
        else:
            if not self.slice_from(u"ați"):
                return False
        return True

    def __r_combo_suffix(self):
        v_1 = self.limit - self.cursor
        self.ket = self.cursor
        among_var = self.find_among_b(RomanianStemmer.a_3)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R1():
            return False
        if among_var == 1:
            if not self.slice_from(u"abil"):
                return False
        elif among_var == 2:
            if not self.slice_from(u"ibil"):
                return False
        elif among_var == 3:
            if not self.slice_from(u"iv"):
                return False
        elif among_var == 4:
            if not self.slice_from(u"ic"):
                return False
        elif among_var == 5:
            if not self.slice_from(u"at"):
                return False
        else:
            if not self.slice_from(u"it"):
                return False
        self.B_standard_suffix_removed = True
        self.cursor = self.limit - v_1
        return True

    def __r_standard_suffix(self):
        self.B_standard_suffix_removed = False
        while True:
            v_1 = self.limit - self.cursor
            try:
                if not self.__r_combo_suffix():
                    raise lab0()
                continue
            except lab0: pass
            self.cursor = self.limit - v_1
            break
        self.ket = self.cursor
        among_var = self.find_among_b(RomanianStemmer.a_4)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R2():
            return False
        if among_var == 1:
            if not self.slice_del():
                return False

        elif among_var == 2:
            if not self.eq_s_b(u"ț"):
                return False
            self.bra = self.cursor
            if not self.slice_from(u"t"):
                return False
        else:
            if not self.slice_from(u"ist"):
                return False
        self.B_standard_suffix_removed = True
        return True

    def __r_verb_suffix(self):
        if self.cursor < self.I_pV:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_pV
        self.ket = self.cursor
        among_var = self.find_among_b(RomanianStemmer.a_5)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        if among_var == 1:
            try:
                v_3 = self.limit - self.cursor
                try:
                    if not self.out_grouping_b(RomanianStemmer.g_v):
                        raise lab1()
                    raise lab0()
                except lab1: pass
                self.cursor = self.limit - v_3
                if not self.eq_s_b(u"u"):
                    self.limit_backward = v_2
                    return False
            except lab0: pass
            if not self.slice_del():
                return False

        else:
            if not self.slice_del():
                return False

        self.limit_backward = v_2
        return True

    def __r_vowel_suffix(self):
        self.ket = self.cursor
        if self.find_among_b(RomanianStemmer.a_6) == 0:
            return False
        self.bra = self.cursor
        if not self.__r_RV():
            return False
        if not self.slice_del():
            return False

        return True

    def _stem(self):
        self.__r_norm()
        v_2 = self.cursor
        self.__r_prelude()
        self.cursor = v_2
        self.__r_mark_regions()
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_4 = self.limit - self.cursor
        self.__r_step_0()
        self.cursor = self.limit - v_4
        v_5 = self.limit - self.cursor
        self.__r_standard_suffix()
        self.cursor = self.limit - v_5
        v_6 = self.limit - self.cursor
        try:
            try:
                v_7 = self.limit - self.cursor
                try:
                    if not self.B_standard_suffix_removed:
                        raise lab2()
                    raise lab1()
                except lab2: pass
                self.cursor = self.limit - v_7
                if not self.__r_verb_suffix():
                    raise lab0()
            except lab1: pass
        except lab0: pass
        self.cursor = self.limit - v_6
        v_8 = self.limit - self.cursor
        self.__r_vowel_suffix()
        self.cursor = self.limit - v_8
        self.cursor = self.limit_backward
        v_9 = self.cursor
        self.__r_postlude()
        self.cursor = v_9
        return True

    a_0 = [
        Among(u"ş", -1, 1),
        Among(u"ţ", -1, 2)
    ]

    a_1 = [
        Among(u"", -1, 3),
        Among(u"I", 0, 1),
        Among(u"U", 0, 2)
    ]

    a_2 = [
        Among(u"ea", -1, 3),
        Among(u"ația", -1, 7),
        Among(u"aua", -1, 2),
        Among(u"iua", -1, 4),
        Among(u"ație", -1, 7),
        Among(u"ele", -1, 3),
        Among(u"ile", -1, 5),
        Among(u"iile", 6, 4),
        Among(u"iei", -1, 4),
        Among(u"atei", -1, 6),
        Among(u"ii", -1, 4),
        Among(u"ului", -1, 1),
        Among(u"ul", -1, 1),
        Among(u"elor", -1, 3),
        Among(u"ilor", -1, 4),
        Among(u"iilor", 14, 4)
    ]

    a_3 = [
        Among(u"icala", -1, 4),
        Among(u"iciva", -1, 4),
        Among(u"ativa", -1, 5),
        Among(u"itiva", -1, 6),
        Among(u"icale", -1, 4),
        Among(u"ațiune", -1, 5),
        Among(u"ițiune", -1, 6),
        Among(u"atoare", -1, 5),
        Among(u"itoare", -1, 6),
        Among(u"ătoare", -1, 5),
        Among(u"icitate", -1, 4),
        Among(u"abilitate", -1, 1),
        Among(u"ibilitate", -1, 2),
        Among(u"ivitate", -1, 3),
        Among(u"icive", -1, 4),
        Among(u"ative", -1, 5),
        Among(u"itive", -1, 6),
        Among(u"icali", -1, 4),
        Among(u"atori", -1, 5),
        Among(u"icatori", 18, 4),
        Among(u"itori", -1, 6),
        Among(u"ători", -1, 5),
        Among(u"icitati", -1, 4),
        Among(u"abilitati", -1, 1),
        Among(u"ivitati", -1, 3),
        Among(u"icivi", -1, 4),
        Among(u"ativi", -1, 5),
        Among(u"itivi", -1, 6),
        Among(u"icităi", -1, 4),
        Among(u"abilităi", -1, 1),
        Among(u"ivităi", -1, 3),
        Among(u"icități", -1, 4),
        Among(u"abilități", -1, 1),
        Among(u"ivități", -1, 3),
        Among(u"ical", -1, 4),
        Among(u"ator", -1, 5),
        Among(u"icator", 35, 4),
        Among(u"itor", -1, 6),
        Among(u"ător", -1, 5),
        Among(u"iciv", -1, 4),
        Among(u"ativ", -1, 5),
        Among(u"itiv", -1, 6),
        Among(u"icală", -1, 4),
        Among(u"icivă", -1, 4),
        Among(u"ativă", -1, 5),
        Among(u"itivă", -1, 6)
    ]

    a_4 = [
        Among(u"ica", -1, 1),
        Among(u"abila", -1, 1),
        Among(u"ibila", -1, 1),
        Among(u"oasa", -1, 1),
        Among(u"ata", -1, 1),
        Among(u"ita", -1, 1),
        Among(u"anta", -1, 1),
        Among(u"ista", -1, 3),
        Among(u"uta", -1, 1),
        Among(u"iva", -1, 1),
        Among(u"ic", -1, 1),
        Among(u"ice", -1, 1),
        Among(u"abile", -1, 1),
        Among(u"ibile", -1, 1),
        Among(u"isme", -1, 3),
        Among(u"iune", -1, 2),
        Among(u"oase", -1, 1),
        Among(u"ate", -1, 1),
        Among(u"itate", 17, 1),
        Among(u"ite", -1, 1),
        Among(u"ante", -1, 1),
        Among(u"iste", -1, 3),
        Among(u"ute", -1, 1),
        Among(u"ive", -1, 1),
        Among(u"ici", -1, 1),
        Among(u"abili", -1, 1),
        Among(u"ibili", -1, 1),
        Among(u"iuni", -1, 2),
        Among(u"atori", -1, 1),
        Among(u"osi", -1, 1),
        Among(u"ati", -1, 1),
        Among(u"itati", 30, 1),
        Among(u"iti", -1, 1),
        Among(u"anti", -1, 1),
        Among(u"isti", -1, 3),
        Among(u"uti", -1, 1),
        Among(u"iști", -1, 3),
        Among(u"ivi", -1, 1),
        Among(u"ităi", -1, 1),
        Among(u"oși", -1, 1),
        Among(u"ități", -1, 1),
        Among(u"abil", -1, 1),
        Among(u"ibil", -1, 1),
        Among(u"ism", -1, 3),
        Among(u"ator", -1, 1),
        Among(u"os", -1, 1),
        Among(u"at", -1, 1),
        Among(u"it", -1, 1),
        Among(u"ant", -1, 1),
        Among(u"ist", -1, 3),
        Among(u"ut", -1, 1),
        Among(u"iv", -1, 1),
        Among(u"ică", -1, 1),
        Among(u"abilă", -1, 1),
        Among(u"ibilă", -1, 1),
        Among(u"oasă", -1, 1),
        Among(u"ată", -1, 1),
        Among(u"ită", -1, 1),
        Among(u"antă", -1, 1),
        Among(u"istă", -1, 3),
        Among(u"ută", -1, 1),
        Among(u"ivă", -1, 1)
    ]

    a_5 = [
        Among(u"ea", -1, 1),
        Among(u"ia", -1, 1),
        Among(u"esc", -1, 1),
        Among(u"ăsc", -1, 1),
        Among(u"ind", -1, 1),
        Among(u"ând", -1, 1),
        Among(u"are", -1, 1),
        Among(u"ere", -1, 1),
        Among(u"ire", -1, 1),
        Among(u"âre", -1, 1),
        Among(u"se", -1, 2),
        Among(u"ase", 10, 1),
        Among(u"sese", 10, 2),
        Among(u"ise", 10, 1),
        Among(u"use", 10, 1),
        Among(u"âse", 10, 1),
        Among(u"ește", -1, 1),
        Among(u"ăște", -1, 1),
        Among(u"eze", -1, 1),
        Among(u"ai", -1, 1),
        Among(u"eai", 19, 1),
        Among(u"iai", 19, 1),
        Among(u"sei", -1, 2),
        Among(u"ești", -1, 1),
        Among(u"ăști", -1, 1),
        Among(u"ui", -1, 1),
        Among(u"ezi", -1, 1),
        Among(u"âi", -1, 1),
        Among(u"ași", -1, 1),
        Among(u"seși", -1, 2),
        Among(u"aseși", 29, 1),
        Among(u"seseși", 29, 2),
        Among(u"iseși", 29, 1),
        Among(u"useși", 29, 1),
        Among(u"âseși", 29, 1),
        Among(u"iși", -1, 1),
        Among(u"uși", -1, 1),
        Among(u"âși", -1, 1),
        Among(u"ați", -1, 2),
        Among(u"eați", 38, 1),
        Among(u"iați", 38, 1),
        Among(u"eți", -1, 2),
        Among(u"iți", -1, 2),
        Among(u"âți", -1, 2),
        Among(u"arăți", -1, 1),
        Among(u"serăți", -1, 2),
        Among(u"aserăți", 45, 1),
        Among(u"seserăți", 45, 2),
        Among(u"iserăți", 45, 1),
        Among(u"userăți", 45, 1),
        Among(u"âserăți", 45, 1),
        Among(u"irăți", -1, 1),
        Among(u"urăți", -1, 1),
        Among(u"ârăți", -1, 1),
        Among(u"am", -1, 1),
        Among(u"eam", 54, 1),
        Among(u"iam", 54, 1),
        Among(u"em", -1, 2),
        Among(u"asem", 57, 1),
        Among(u"sesem", 57, 2),
        Among(u"isem", 57, 1),
        Among(u"usem", 57, 1),
        Among(u"âsem", 57, 1),
        Among(u"im", -1, 2),
        Among(u"âm", -1, 2),
        Among(u"ăm", -1, 2),
        Among(u"arăm", 65, 1),
        Among(u"serăm", 65, 2),
        Among(u"aserăm", 67, 1),
        Among(u"seserăm", 67, 2),
        Among(u"iserăm", 67, 1),
        Among(u"userăm", 67, 1),
        Among(u"âserăm", 67, 1),
        Among(u"irăm", 65, 1),
        Among(u"urăm", 65, 1),
        Among(u"ârăm", 65, 1),
        Among(u"au", -1, 1),
        Among(u"eau", 76, 1),
        Among(u"iau", 76, 1),
        Among(u"indu", -1, 1),
        Among(u"ându", -1, 1),
        Among(u"ez", -1, 1),
        Among(u"ească", -1, 1),
        Among(u"ară", -1, 1),
        Among(u"seră", -1, 2),
        Among(u"aseră", 84, 1),
        Among(u"seseră", 84, 2),
        Among(u"iseră", 84, 1),
        Among(u"useră", 84, 1),
        Among(u"âseră", 84, 1),
        Among(u"iră", -1, 1),
        Among(u"ură", -1, 1),
        Among(u"âră", -1, 1),
        Among(u"ează", -1, 1)
    ]

    a_6 = [
        Among(u"a", -1, 1),
        Among(u"e", -1, 1),
        Among(u"ie", 1, 1),
        Among(u"i", -1, 1),
        Among(u"ă", -1, 1)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass


class lab2(BaseException): pass


class lab3(BaseException): pass


class lab4(BaseException): pass


class lab5(BaseException): pass


class lab6(BaseException): pass


class lab7(BaseException): pass

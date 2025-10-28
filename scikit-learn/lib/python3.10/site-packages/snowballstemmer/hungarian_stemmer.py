#-*- coding: utf-8 -*-
# Generated from hungarian.sbl by Snowball 3.0.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class HungarianStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from hungarian.sbl by Snowball 3.0.1 - https://snowballstem.org/
    '''

    g_v = {u"a", u"e", u"i", u"o", u"u", u"á", u"é", u"í", u"ó", u"ö", u"ő", u"ú", u"ü", u"ű"}

    I_p1 = 0

    def __r_mark_regions(self):
        self.I_p1 = self.limit
        try:
            v_1 = self.cursor
            try:
                if not self.in_grouping(HungarianStemmer.g_v):
                    raise lab1()
                v_2 = self.cursor
                try:
                    if not self.go_in_grouping(HungarianStemmer.g_v):
                        raise lab2()
                    self.cursor += 1
                    self.I_p1 = self.cursor
                except lab2: pass
                self.cursor = v_2
                raise lab0()
            except lab1: pass
            self.cursor = v_1
            if not self.go_out_grouping(HungarianStemmer.g_v):
                return False
            self.cursor += 1
            self.I_p1 = self.cursor
        except lab0: pass
        return True

    def __r_R1(self):
        return self.I_p1 <= self.cursor

    def __r_v_ending(self):
        self.ket = self.cursor
        among_var = self.find_among_b(HungarianStemmer.a_0)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R1():
            return False
        if among_var == 1:
            if not self.slice_from(u"a"):
                return False
        else:
            if not self.slice_from(u"e"):
                return False
        return True

    def __r_double(self):
        v_1 = self.limit - self.cursor
        if self.find_among_b(HungarianStemmer.a_1) == 0:
            return False
        self.cursor = self.limit - v_1
        return True

    def __r_undouble(self):
        if self.cursor <= self.limit_backward:
            return False
        self.cursor -= 1
        self.ket = self.cursor
        if self.cursor <= self.limit_backward:
            return False
        self.cursor -= 1
        self.bra = self.cursor
        if not self.slice_del():
            return False

        return True

    def __r_instrum(self):
        self.ket = self.cursor
        if self.find_among_b(HungarianStemmer.a_2) == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R1():
            return False
        if not self.__r_double():
            return False
        if not self.slice_del():
            return False

        if not self.__r_undouble():
            return False
        return True

    def __r_case(self):
        self.ket = self.cursor
        if self.find_among_b(HungarianStemmer.a_3) == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R1():
            return False
        if not self.slice_del():
            return False

        if not self.__r_v_ending():
            return False
        return True

    def __r_case_special(self):
        self.ket = self.cursor
        among_var = self.find_among_b(HungarianStemmer.a_4)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R1():
            return False
        if among_var == 1:
            if not self.slice_from(u"e"):
                return False
        else:
            if not self.slice_from(u"a"):
                return False
        return True

    def __r_case_other(self):
        self.ket = self.cursor
        among_var = self.find_among_b(HungarianStemmer.a_5)
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
        else:
            if not self.slice_from(u"e"):
                return False
        return True

    def __r_factive(self):
        self.ket = self.cursor
        if self.find_among_b(HungarianStemmer.a_6) == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R1():
            return False
        if not self.__r_double():
            return False
        if not self.slice_del():
            return False

        if not self.__r_undouble():
            return False
        return True

    def __r_plural(self):
        self.ket = self.cursor
        among_var = self.find_among_b(HungarianStemmer.a_7)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R1():
            return False
        if among_var == 1:
            if not self.slice_from(u"a"):
                return False
        elif among_var == 2:
            if not self.slice_from(u"e"):
                return False
        else:
            if not self.slice_del():
                return False

        return True

    def __r_owned(self):
        self.ket = self.cursor
        among_var = self.find_among_b(HungarianStemmer.a_8)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R1():
            return False
        if among_var == 1:
            if not self.slice_del():
                return False

        elif among_var == 2:
            if not self.slice_from(u"e"):
                return False
        else:
            if not self.slice_from(u"a"):
                return False
        return True

    def __r_sing_owner(self):
        self.ket = self.cursor
        among_var = self.find_among_b(HungarianStemmer.a_9)
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
        else:
            if not self.slice_from(u"e"):
                return False
        return True

    def __r_plur_owner(self):
        self.ket = self.cursor
        among_var = self.find_among_b(HungarianStemmer.a_10)
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
        else:
            if not self.slice_from(u"e"):
                return False
        return True

    def _stem(self):
        v_1 = self.cursor
        self.__r_mark_regions()
        self.cursor = v_1
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_2 = self.limit - self.cursor
        self.__r_instrum()
        self.cursor = self.limit - v_2
        v_3 = self.limit - self.cursor
        self.__r_case()
        self.cursor = self.limit - v_3
        v_4 = self.limit - self.cursor
        self.__r_case_special()
        self.cursor = self.limit - v_4
        v_5 = self.limit - self.cursor
        self.__r_case_other()
        self.cursor = self.limit - v_5
        v_6 = self.limit - self.cursor
        self.__r_factive()
        self.cursor = self.limit - v_6
        v_7 = self.limit - self.cursor
        self.__r_owned()
        self.cursor = self.limit - v_7
        v_8 = self.limit - self.cursor
        self.__r_sing_owner()
        self.cursor = self.limit - v_8
        v_9 = self.limit - self.cursor
        self.__r_plur_owner()
        self.cursor = self.limit - v_9
        v_10 = self.limit - self.cursor
        self.__r_plural()
        self.cursor = self.limit - v_10
        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among(u"á", -1, 1),
        Among(u"é", -1, 2)
    ]

    a_1 = [
        Among(u"bb", -1, -1),
        Among(u"cc", -1, -1),
        Among(u"dd", -1, -1),
        Among(u"ff", -1, -1),
        Among(u"gg", -1, -1),
        Among(u"jj", -1, -1),
        Among(u"kk", -1, -1),
        Among(u"ll", -1, -1),
        Among(u"mm", -1, -1),
        Among(u"nn", -1, -1),
        Among(u"pp", -1, -1),
        Among(u"rr", -1, -1),
        Among(u"ccs", -1, -1),
        Among(u"ss", -1, -1),
        Among(u"zzs", -1, -1),
        Among(u"tt", -1, -1),
        Among(u"vv", -1, -1),
        Among(u"ggy", -1, -1),
        Among(u"lly", -1, -1),
        Among(u"nny", -1, -1),
        Among(u"tty", -1, -1),
        Among(u"ssz", -1, -1),
        Among(u"zz", -1, -1)
    ]

    a_2 = [
        Among(u"al", -1, 1),
        Among(u"el", -1, 1)
    ]

    a_3 = [
        Among(u"ba", -1, -1),
        Among(u"ra", -1, -1),
        Among(u"be", -1, -1),
        Among(u"re", -1, -1),
        Among(u"ig", -1, -1),
        Among(u"nak", -1, -1),
        Among(u"nek", -1, -1),
        Among(u"val", -1, -1),
        Among(u"vel", -1, -1),
        Among(u"ul", -1, -1),
        Among(u"nál", -1, -1),
        Among(u"nél", -1, -1),
        Among(u"ból", -1, -1),
        Among(u"ról", -1, -1),
        Among(u"tól", -1, -1),
        Among(u"ül", -1, -1),
        Among(u"ből", -1, -1),
        Among(u"ről", -1, -1),
        Among(u"től", -1, -1),
        Among(u"n", -1, -1),
        Among(u"an", 19, -1),
        Among(u"ban", 20, -1),
        Among(u"en", 19, -1),
        Among(u"ben", 22, -1),
        Among(u"képpen", 22, -1),
        Among(u"on", 19, -1),
        Among(u"ön", 19, -1),
        Among(u"képp", -1, -1),
        Among(u"kor", -1, -1),
        Among(u"t", -1, -1),
        Among(u"at", 29, -1),
        Among(u"et", 29, -1),
        Among(u"ként", 29, -1),
        Among(u"anként", 32, -1),
        Among(u"enként", 32, -1),
        Among(u"onként", 32, -1),
        Among(u"ot", 29, -1),
        Among(u"ért", 29, -1),
        Among(u"öt", 29, -1),
        Among(u"hez", -1, -1),
        Among(u"hoz", -1, -1),
        Among(u"höz", -1, -1),
        Among(u"vá", -1, -1),
        Among(u"vé", -1, -1)
    ]

    a_4 = [
        Among(u"án", -1, 2),
        Among(u"én", -1, 1),
        Among(u"ánként", -1, 2)
    ]

    a_5 = [
        Among(u"stul", -1, 1),
        Among(u"astul", 0, 1),
        Among(u"ástul", 0, 2),
        Among(u"stül", -1, 1),
        Among(u"estül", 3, 1),
        Among(u"éstül", 3, 3)
    ]

    a_6 = [
        Among(u"á", -1, 1),
        Among(u"é", -1, 1)
    ]

    a_7 = [
        Among(u"k", -1, 3),
        Among(u"ak", 0, 3),
        Among(u"ek", 0, 3),
        Among(u"ok", 0, 3),
        Among(u"ák", 0, 1),
        Among(u"ék", 0, 2),
        Among(u"ök", 0, 3)
    ]

    a_8 = [
        Among(u"éi", -1, 1),
        Among(u"áéi", 0, 3),
        Among(u"ééi", 0, 2),
        Among(u"é", -1, 1),
        Among(u"ké", 3, 1),
        Among(u"aké", 4, 1),
        Among(u"eké", 4, 1),
        Among(u"oké", 4, 1),
        Among(u"áké", 4, 3),
        Among(u"éké", 4, 2),
        Among(u"öké", 4, 1),
        Among(u"éé", 3, 2)
    ]

    a_9 = [
        Among(u"a", -1, 1),
        Among(u"ja", 0, 1),
        Among(u"d", -1, 1),
        Among(u"ad", 2, 1),
        Among(u"ed", 2, 1),
        Among(u"od", 2, 1),
        Among(u"ád", 2, 2),
        Among(u"éd", 2, 3),
        Among(u"öd", 2, 1),
        Among(u"e", -1, 1),
        Among(u"je", 9, 1),
        Among(u"nk", -1, 1),
        Among(u"unk", 11, 1),
        Among(u"ánk", 11, 2),
        Among(u"énk", 11, 3),
        Among(u"ünk", 11, 1),
        Among(u"uk", -1, 1),
        Among(u"juk", 16, 1),
        Among(u"ájuk", 17, 2),
        Among(u"ük", -1, 1),
        Among(u"jük", 19, 1),
        Among(u"éjük", 20, 3),
        Among(u"m", -1, 1),
        Among(u"am", 22, 1),
        Among(u"em", 22, 1),
        Among(u"om", 22, 1),
        Among(u"ám", 22, 2),
        Among(u"ém", 22, 3),
        Among(u"o", -1, 1),
        Among(u"á", -1, 2),
        Among(u"é", -1, 3)
    ]

    a_10 = [
        Among(u"id", -1, 1),
        Among(u"aid", 0, 1),
        Among(u"jaid", 1, 1),
        Among(u"eid", 0, 1),
        Among(u"jeid", 3, 1),
        Among(u"áid", 0, 2),
        Among(u"éid", 0, 3),
        Among(u"i", -1, 1),
        Among(u"ai", 7, 1),
        Among(u"jai", 8, 1),
        Among(u"ei", 7, 1),
        Among(u"jei", 10, 1),
        Among(u"ái", 7, 2),
        Among(u"éi", 7, 3),
        Among(u"itek", -1, 1),
        Among(u"eitek", 14, 1),
        Among(u"jeitek", 15, 1),
        Among(u"éitek", 14, 3),
        Among(u"ik", -1, 1),
        Among(u"aik", 18, 1),
        Among(u"jaik", 19, 1),
        Among(u"eik", 18, 1),
        Among(u"jeik", 21, 1),
        Among(u"áik", 18, 2),
        Among(u"éik", 18, 3),
        Among(u"ink", -1, 1),
        Among(u"aink", 25, 1),
        Among(u"jaink", 26, 1),
        Among(u"eink", 25, 1),
        Among(u"jeink", 28, 1),
        Among(u"áink", 25, 2),
        Among(u"éink", 25, 3),
        Among(u"aitok", -1, 1),
        Among(u"jaitok", 32, 1),
        Among(u"áitok", -1, 2),
        Among(u"im", -1, 1),
        Among(u"aim", 35, 1),
        Among(u"jaim", 36, 1),
        Among(u"eim", 35, 1),
        Among(u"jeim", 38, 1),
        Among(u"áim", 35, 2),
        Among(u"éim", 35, 3)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass


class lab2(BaseException): pass

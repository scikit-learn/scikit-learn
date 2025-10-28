#-*- coding: utf-8 -*-
# Generated from estonian.sbl by Snowball 3.0.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class EstonianStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from estonian.sbl by Snowball 3.0.1 - https://snowballstem.org/
    '''

    g_V1 = {u"a", u"e", u"i", u"o", u"u", u"õ", u"ä", u"ö", u"ü"}

    g_RV = {u"a", u"e", u"i", u"u", u"o"}

    g_KI = {u"k", u"p", u"t", u"g", u"b", u"d", u"s", u"h", u"f", u"š", u"z", u"ž"}

    g_GI = {u"c", u"j", u"l", u"m", u"n", u"q", u"r", u"v", u"w", u"x", u"a", u"e", u"i", u"o", u"u", u"õ", u"ä", u"ö", u"ü"}

    I_p1 = 0

    def __r_mark_regions(self):
        self.I_p1 = self.limit
        if not self.go_out_grouping(EstonianStemmer.g_V1):
            return False
        self.cursor += 1
        if not self.go_in_grouping(EstonianStemmer.g_V1):
            return False
        self.cursor += 1
        self.I_p1 = self.cursor
        return True

    def __r_emphasis(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        among_var = self.find_among_b(EstonianStemmer.a_0)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        v_3 = self.limit - self.cursor
        c = self.cursor - 4
        if c < self.limit_backward:
            return False
        self.cursor = c
        self.cursor = self.limit - v_3
        if among_var == 1:
            v_4 = self.limit - self.cursor
            if not self.in_grouping_b(EstonianStemmer.g_GI):
                return False
            self.cursor = self.limit - v_4
            v_5 = self.limit - self.cursor
            try:
                if not self.__r_LONGV():
                    raise lab0()
                return False
            except lab0: pass
            self.cursor = self.limit - v_5
            if not self.slice_del():
                return False

        else:
            if not self.in_grouping_b(EstonianStemmer.g_KI):
                return False
            if not self.slice_del():
                return False

        return True

    def __r_verb(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        among_var = self.find_among_b(EstonianStemmer.a_1)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if among_var == 1:
            if not self.slice_del():
                return False

        elif among_var == 2:
            if not self.slice_from(u"a"):
                return False
        else:
            if not self.in_grouping_b(EstonianStemmer.g_V1):
                return False
            if not self.slice_del():
                return False

        return True

    def __r_LONGV(self):
        if self.find_among_b(EstonianStemmer.a_2) == 0:
            return False
        return True

    def __r_i_plural(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        if self.find_among_b(EstonianStemmer.a_3) == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if not self.in_grouping_b(EstonianStemmer.g_RV):
            return False
        if not self.slice_del():
            return False

        return True

    def __r_special_noun_endings(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        among_var = self.find_among_b(EstonianStemmer.a_4)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if among_var == 1:
            if not self.slice_from(u"lase"):
                return False
        elif among_var == 2:
            if not self.slice_from(u"mise"):
                return False
        else:
            if not self.slice_from(u"lise"):
                return False
        return True

    def __r_case_ending(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        among_var = self.find_among_b(EstonianStemmer.a_5)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if among_var == 1:
            try:
                v_3 = self.limit - self.cursor
                try:
                    if not self.in_grouping_b(EstonianStemmer.g_RV):
                        raise lab1()
                    raise lab0()
                except lab1: pass
                self.cursor = self.limit - v_3
                if not self.__r_LONGV():
                    return False
            except lab0: pass
        else:
            v_4 = self.limit - self.cursor
            c = self.cursor - 4
            if c < self.limit_backward:
                return False
            self.cursor = c
            self.cursor = self.limit - v_4
        if not self.slice_del():
            return False

        return True

    def __r_plural_three_first_cases(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        among_var = self.find_among_b(EstonianStemmer.a_7)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if among_var == 1:
            if not self.slice_from(u"iku"):
                return False
        elif among_var == 2:
            v_3 = self.limit - self.cursor
            try:
                if not self.__r_LONGV():
                    raise lab0()
                return False
            except lab0: pass
            self.cursor = self.limit - v_3
            if not self.slice_del():
                return False

        elif among_var == 3:
            try:
                v_4 = self.limit - self.cursor
                try:
                    v_5 = self.limit - self.cursor
                    c = self.cursor - 4
                    if c < self.limit_backward:
                        raise lab2()
                    self.cursor = c
                    self.cursor = self.limit - v_5
                    among_var = self.find_among_b(EstonianStemmer.a_6)
                    if among_var == 1:
                        if not self.slice_from(u"e"):
                            return False
                    elif among_var == 2:
                        if not self.slice_del():
                            return False

                    raise lab1()
                except lab2: pass
                self.cursor = self.limit - v_4
                if not self.slice_from(u"t"):
                    return False
            except lab1: pass
        else:
            try:
                v_6 = self.limit - self.cursor
                try:
                    if not self.in_grouping_b(EstonianStemmer.g_RV):
                        raise lab4()
                    raise lab3()
                except lab4: pass
                self.cursor = self.limit - v_6
                if not self.__r_LONGV():
                    return False
            except lab3: pass
            if not self.slice_del():
                return False

        return True

    def __r_nu(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        if self.find_among_b(EstonianStemmer.a_8) == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if not self.slice_del():
            return False

        return True

    def __r_undouble_kpt(self):
        if not self.in_grouping_b(EstonianStemmer.g_V1):
            return False
        if self.I_p1 > self.cursor:
            return False
        self.ket = self.cursor
        among_var = self.find_among_b(EstonianStemmer.a_9)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if not self.slice_from(u"k"):
                return False
        elif among_var == 2:
            if not self.slice_from(u"p"):
                return False
        else:
            if not self.slice_from(u"t"):
                return False
        return True

    def __r_degrees(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        among_var = self.find_among_b(EstonianStemmer.a_10)
        if among_var == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if among_var == 1:
            if not self.in_grouping_b(EstonianStemmer.g_RV):
                return False
            if not self.slice_del():
                return False

        else:
            if not self.slice_del():
                return False

        return True

    def __r_substantive(self):
        v_1 = self.limit - self.cursor
        self.__r_special_noun_endings()
        self.cursor = self.limit - v_1
        v_2 = self.limit - self.cursor
        self.__r_case_ending()
        self.cursor = self.limit - v_2
        v_3 = self.limit - self.cursor
        self.__r_plural_three_first_cases()
        self.cursor = self.limit - v_3
        v_4 = self.limit - self.cursor
        self.__r_degrees()
        self.cursor = self.limit - v_4
        v_5 = self.limit - self.cursor
        self.__r_i_plural()
        self.cursor = self.limit - v_5
        v_6 = self.limit - self.cursor
        self.__r_nu()
        self.cursor = self.limit - v_6
        return True

    def __r_verb_exceptions(self):
        self.bra = self.cursor
        among_var = self.find_among(EstonianStemmer.a_11)
        if among_var == 0:
            return False
        self.ket = self.cursor
        if self.cursor < self.limit:
            return False
        if among_var == 1:
            if not self.slice_from(u"joo"):
                return False
        elif among_var == 2:
            if not self.slice_from(u"saa"):
                return False
        elif among_var == 3:
            if not self.slice_from(u"viima"):
                return False
        elif among_var == 4:
            if not self.slice_from(u"keesi"):
                return False
        elif among_var == 5:
            if not self.slice_from(u"löö"):
                return False
        elif among_var == 6:
            if not self.slice_from(u"lõi"):
                return False
        elif among_var == 7:
            if not self.slice_from(u"loo"):
                return False
        elif among_var == 8:
            if not self.slice_from(u"käisi"):
                return False
        elif among_var == 9:
            if not self.slice_from(u"söö"):
                return False
        elif among_var == 10:
            if not self.slice_from(u"too"):
                return False
        elif among_var == 11:
            if not self.slice_from(u"võisi"):
                return False
        elif among_var == 12:
            if not self.slice_from(u"jääma"):
                return False
        elif among_var == 13:
            if not self.slice_from(u"müüsi"):
                return False
        elif among_var == 14:
            if not self.slice_from(u"luge"):
                return False
        elif among_var == 15:
            if not self.slice_from(u"põde"):
                return False
        elif among_var == 16:
            if not self.slice_from(u"ladu"):
                return False
        elif among_var == 17:
            if not self.slice_from(u"tegi"):
                return False
        else:
            if not self.slice_from(u"nägi"):
                return False
        return True

    def _stem(self):
        v_1 = self.cursor
        try:
            if not self.__r_verb_exceptions():
                raise lab0()
            return False
        except lab0: pass
        self.cursor = v_1
        v_2 = self.cursor
        self.__r_mark_regions()
        self.cursor = v_2
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_3 = self.limit - self.cursor
        self.__r_emphasis()
        self.cursor = self.limit - v_3
        v_4 = self.limit - self.cursor
        try:
            try:
                v_5 = self.limit - self.cursor
                try:
                    if not self.__r_verb():
                        raise lab3()
                    raise lab2()
                except lab3: pass
                self.cursor = self.limit - v_5
                self.__r_substantive()
            except lab2: pass
        except lab1: pass
        self.cursor = self.limit - v_4
        v_6 = self.limit - self.cursor
        self.__r_undouble_kpt()
        self.cursor = self.limit - v_6
        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among(u"gi", -1, 1),
        Among(u"ki", -1, 2)
    ]

    a_1 = [
        Among(u"da", -1, 3),
        Among(u"mata", -1, 1),
        Among(u"b", -1, 3),
        Among(u"ksid", -1, 1),
        Among(u"nuksid", 3, 1),
        Among(u"me", -1, 3),
        Among(u"sime", 5, 1),
        Among(u"ksime", 6, 1),
        Among(u"nuksime", 7, 1),
        Among(u"akse", -1, 2),
        Among(u"dakse", 9, 1),
        Among(u"takse", 9, 1),
        Among(u"site", -1, 1),
        Among(u"ksite", 12, 1),
        Among(u"nuksite", 13, 1),
        Among(u"n", -1, 3),
        Among(u"sin", 15, 1),
        Among(u"ksin", 16, 1),
        Among(u"nuksin", 17, 1),
        Among(u"daks", -1, 1),
        Among(u"taks", -1, 1)
    ]

    a_2 = [
        Among(u"aa", -1, -1),
        Among(u"ee", -1, -1),
        Among(u"ii", -1, -1),
        Among(u"oo", -1, -1),
        Among(u"uu", -1, -1),
        Among(u"ää", -1, -1),
        Among(u"õõ", -1, -1),
        Among(u"öö", -1, -1),
        Among(u"üü", -1, -1)
    ]

    a_3 = [
        Among(u"i", -1, 1)
    ]

    a_4 = [
        Among(u"lane", -1, 1),
        Among(u"line", -1, 3),
        Among(u"mine", -1, 2),
        Among(u"lasse", -1, 1),
        Among(u"lisse", -1, 3),
        Among(u"misse", -1, 2),
        Among(u"lasi", -1, 1),
        Among(u"lisi", -1, 3),
        Among(u"misi", -1, 2),
        Among(u"last", -1, 1),
        Among(u"list", -1, 3),
        Among(u"mist", -1, 2)
    ]

    a_5 = [
        Among(u"ga", -1, 1),
        Among(u"ta", -1, 1),
        Among(u"le", -1, 1),
        Among(u"sse", -1, 1),
        Among(u"l", -1, 1),
        Among(u"s", -1, 1),
        Among(u"ks", 5, 1),
        Among(u"t", -1, 2),
        Among(u"lt", 7, 1),
        Among(u"st", 7, 1)
    ]

    a_6 = [
        Among(u"", -1, 2),
        Among(u"las", 0, 1),
        Among(u"lis", 0, 1),
        Among(u"mis", 0, 1),
        Among(u"t", 0, -1)
    ]

    a_7 = [
        Among(u"d", -1, 4),
        Among(u"sid", 0, 2),
        Among(u"de", -1, 4),
        Among(u"ikkude", 2, 1),
        Among(u"ike", -1, 1),
        Among(u"ikke", -1, 1),
        Among(u"te", -1, 3)
    ]

    a_8 = [
        Among(u"va", -1, -1),
        Among(u"du", -1, -1),
        Among(u"nu", -1, -1),
        Among(u"tu", -1, -1)
    ]

    a_9 = [
        Among(u"kk", -1, 1),
        Among(u"pp", -1, 2),
        Among(u"tt", -1, 3)
    ]

    a_10 = [
        Among(u"ma", -1, 2),
        Among(u"mai", -1, 1),
        Among(u"m", -1, 1)
    ]

    a_11 = [
        Among(u"joob", -1, 1),
        Among(u"jood", -1, 1),
        Among(u"joodakse", 1, 1),
        Among(u"jooma", -1, 1),
        Among(u"joomata", 3, 1),
        Among(u"joome", -1, 1),
        Among(u"joon", -1, 1),
        Among(u"joote", -1, 1),
        Among(u"joovad", -1, 1),
        Among(u"juua", -1, 1),
        Among(u"juuakse", 9, 1),
        Among(u"jäi", -1, 12),
        Among(u"jäid", 11, 12),
        Among(u"jäime", 11, 12),
        Among(u"jäin", 11, 12),
        Among(u"jäite", 11, 12),
        Among(u"jääb", -1, 12),
        Among(u"jääd", -1, 12),
        Among(u"jääda", 17, 12),
        Among(u"jäädakse", 18, 12),
        Among(u"jäädi", 17, 12),
        Among(u"jääks", -1, 12),
        Among(u"jääksid", 21, 12),
        Among(u"jääksime", 21, 12),
        Among(u"jääksin", 21, 12),
        Among(u"jääksite", 21, 12),
        Among(u"jääma", -1, 12),
        Among(u"jäämata", 26, 12),
        Among(u"jääme", -1, 12),
        Among(u"jään", -1, 12),
        Among(u"jääte", -1, 12),
        Among(u"jäävad", -1, 12),
        Among(u"jõi", -1, 1),
        Among(u"jõid", 32, 1),
        Among(u"jõime", 32, 1),
        Among(u"jõin", 32, 1),
        Among(u"jõite", 32, 1),
        Among(u"keeb", -1, 4),
        Among(u"keed", -1, 4),
        Among(u"keedakse", 38, 4),
        Among(u"keeks", -1, 4),
        Among(u"keeksid", 40, 4),
        Among(u"keeksime", 40, 4),
        Among(u"keeksin", 40, 4),
        Among(u"keeksite", 40, 4),
        Among(u"keema", -1, 4),
        Among(u"keemata", 45, 4),
        Among(u"keeme", -1, 4),
        Among(u"keen", -1, 4),
        Among(u"kees", -1, 4),
        Among(u"keeta", -1, 4),
        Among(u"keete", -1, 4),
        Among(u"keevad", -1, 4),
        Among(u"käia", -1, 8),
        Among(u"käiakse", 53, 8),
        Among(u"käib", -1, 8),
        Among(u"käid", -1, 8),
        Among(u"käidi", 56, 8),
        Among(u"käiks", -1, 8),
        Among(u"käiksid", 58, 8),
        Among(u"käiksime", 58, 8),
        Among(u"käiksin", 58, 8),
        Among(u"käiksite", 58, 8),
        Among(u"käima", -1, 8),
        Among(u"käimata", 63, 8),
        Among(u"käime", -1, 8),
        Among(u"käin", -1, 8),
        Among(u"käis", -1, 8),
        Among(u"käite", -1, 8),
        Among(u"käivad", -1, 8),
        Among(u"laob", -1, 16),
        Among(u"laod", -1, 16),
        Among(u"laoks", -1, 16),
        Among(u"laoksid", 72, 16),
        Among(u"laoksime", 72, 16),
        Among(u"laoksin", 72, 16),
        Among(u"laoksite", 72, 16),
        Among(u"laome", -1, 16),
        Among(u"laon", -1, 16),
        Among(u"laote", -1, 16),
        Among(u"laovad", -1, 16),
        Among(u"loeb", -1, 14),
        Among(u"loed", -1, 14),
        Among(u"loeks", -1, 14),
        Among(u"loeksid", 83, 14),
        Among(u"loeksime", 83, 14),
        Among(u"loeksin", 83, 14),
        Among(u"loeksite", 83, 14),
        Among(u"loeme", -1, 14),
        Among(u"loen", -1, 14),
        Among(u"loete", -1, 14),
        Among(u"loevad", -1, 14),
        Among(u"loob", -1, 7),
        Among(u"lood", -1, 7),
        Among(u"loodi", 93, 7),
        Among(u"looks", -1, 7),
        Among(u"looksid", 95, 7),
        Among(u"looksime", 95, 7),
        Among(u"looksin", 95, 7),
        Among(u"looksite", 95, 7),
        Among(u"looma", -1, 7),
        Among(u"loomata", 100, 7),
        Among(u"loome", -1, 7),
        Among(u"loon", -1, 7),
        Among(u"loote", -1, 7),
        Among(u"loovad", -1, 7),
        Among(u"luua", -1, 7),
        Among(u"luuakse", 106, 7),
        Among(u"lõi", -1, 6),
        Among(u"lõid", 108, 6),
        Among(u"lõime", 108, 6),
        Among(u"lõin", 108, 6),
        Among(u"lõite", 108, 6),
        Among(u"lööb", -1, 5),
        Among(u"lööd", -1, 5),
        Among(u"löödakse", 114, 5),
        Among(u"löödi", 114, 5),
        Among(u"lööks", -1, 5),
        Among(u"lööksid", 117, 5),
        Among(u"lööksime", 117, 5),
        Among(u"lööksin", 117, 5),
        Among(u"lööksite", 117, 5),
        Among(u"lööma", -1, 5),
        Among(u"löömata", 122, 5),
        Among(u"lööme", -1, 5),
        Among(u"löön", -1, 5),
        Among(u"lööte", -1, 5),
        Among(u"löövad", -1, 5),
        Among(u"lüüa", -1, 5),
        Among(u"lüüakse", 128, 5),
        Among(u"müüa", -1, 13),
        Among(u"müüakse", 130, 13),
        Among(u"müüb", -1, 13),
        Among(u"müüd", -1, 13),
        Among(u"müüdi", 133, 13),
        Among(u"müüks", -1, 13),
        Among(u"müüksid", 135, 13),
        Among(u"müüksime", 135, 13),
        Among(u"müüksin", 135, 13),
        Among(u"müüksite", 135, 13),
        Among(u"müüma", -1, 13),
        Among(u"müümata", 140, 13),
        Among(u"müüme", -1, 13),
        Among(u"müün", -1, 13),
        Among(u"müüs", -1, 13),
        Among(u"müüte", -1, 13),
        Among(u"müüvad", -1, 13),
        Among(u"näeb", -1, 18),
        Among(u"näed", -1, 18),
        Among(u"näeks", -1, 18),
        Among(u"näeksid", 149, 18),
        Among(u"näeksime", 149, 18),
        Among(u"näeksin", 149, 18),
        Among(u"näeksite", 149, 18),
        Among(u"näeme", -1, 18),
        Among(u"näen", -1, 18),
        Among(u"näete", -1, 18),
        Among(u"näevad", -1, 18),
        Among(u"nägema", -1, 18),
        Among(u"nägemata", 158, 18),
        Among(u"näha", -1, 18),
        Among(u"nähakse", 160, 18),
        Among(u"nähti", -1, 18),
        Among(u"põeb", -1, 15),
        Among(u"põed", -1, 15),
        Among(u"põeks", -1, 15),
        Among(u"põeksid", 165, 15),
        Among(u"põeksime", 165, 15),
        Among(u"põeksin", 165, 15),
        Among(u"põeksite", 165, 15),
        Among(u"põeme", -1, 15),
        Among(u"põen", -1, 15),
        Among(u"põete", -1, 15),
        Among(u"põevad", -1, 15),
        Among(u"saab", -1, 2),
        Among(u"saad", -1, 2),
        Among(u"saada", 175, 2),
        Among(u"saadakse", 176, 2),
        Among(u"saadi", 175, 2),
        Among(u"saaks", -1, 2),
        Among(u"saaksid", 179, 2),
        Among(u"saaksime", 179, 2),
        Among(u"saaksin", 179, 2),
        Among(u"saaksite", 179, 2),
        Among(u"saama", -1, 2),
        Among(u"saamata", 184, 2),
        Among(u"saame", -1, 2),
        Among(u"saan", -1, 2),
        Among(u"saate", -1, 2),
        Among(u"saavad", -1, 2),
        Among(u"sai", -1, 2),
        Among(u"said", 190, 2),
        Among(u"saime", 190, 2),
        Among(u"sain", 190, 2),
        Among(u"saite", 190, 2),
        Among(u"sõi", -1, 9),
        Among(u"sõid", 195, 9),
        Among(u"sõime", 195, 9),
        Among(u"sõin", 195, 9),
        Among(u"sõite", 195, 9),
        Among(u"sööb", -1, 9),
        Among(u"sööd", -1, 9),
        Among(u"söödakse", 201, 9),
        Among(u"söödi", 201, 9),
        Among(u"sööks", -1, 9),
        Among(u"sööksid", 204, 9),
        Among(u"sööksime", 204, 9),
        Among(u"sööksin", 204, 9),
        Among(u"sööksite", 204, 9),
        Among(u"sööma", -1, 9),
        Among(u"söömata", 209, 9),
        Among(u"sööme", -1, 9),
        Among(u"söön", -1, 9),
        Among(u"sööte", -1, 9),
        Among(u"söövad", -1, 9),
        Among(u"süüa", -1, 9),
        Among(u"süüakse", 215, 9),
        Among(u"teeb", -1, 17),
        Among(u"teed", -1, 17),
        Among(u"teeks", -1, 17),
        Among(u"teeksid", 219, 17),
        Among(u"teeksime", 219, 17),
        Among(u"teeksin", 219, 17),
        Among(u"teeksite", 219, 17),
        Among(u"teeme", -1, 17),
        Among(u"teen", -1, 17),
        Among(u"teete", -1, 17),
        Among(u"teevad", -1, 17),
        Among(u"tegema", -1, 17),
        Among(u"tegemata", 228, 17),
        Among(u"teha", -1, 17),
        Among(u"tehakse", 230, 17),
        Among(u"tehti", -1, 17),
        Among(u"toob", -1, 10),
        Among(u"tood", -1, 10),
        Among(u"toodi", 234, 10),
        Among(u"tooks", -1, 10),
        Among(u"tooksid", 236, 10),
        Among(u"tooksime", 236, 10),
        Among(u"tooksin", 236, 10),
        Among(u"tooksite", 236, 10),
        Among(u"tooma", -1, 10),
        Among(u"toomata", 241, 10),
        Among(u"toome", -1, 10),
        Among(u"toon", -1, 10),
        Among(u"toote", -1, 10),
        Among(u"toovad", -1, 10),
        Among(u"tuua", -1, 10),
        Among(u"tuuakse", 247, 10),
        Among(u"tõi", -1, 10),
        Among(u"tõid", 249, 10),
        Among(u"tõime", 249, 10),
        Among(u"tõin", 249, 10),
        Among(u"tõite", 249, 10),
        Among(u"viia", -1, 3),
        Among(u"viiakse", 254, 3),
        Among(u"viib", -1, 3),
        Among(u"viid", -1, 3),
        Among(u"viidi", 257, 3),
        Among(u"viiks", -1, 3),
        Among(u"viiksid", 259, 3),
        Among(u"viiksime", 259, 3),
        Among(u"viiksin", 259, 3),
        Among(u"viiksite", 259, 3),
        Among(u"viima", -1, 3),
        Among(u"viimata", 264, 3),
        Among(u"viime", -1, 3),
        Among(u"viin", -1, 3),
        Among(u"viisime", -1, 3),
        Among(u"viisin", -1, 3),
        Among(u"viisite", -1, 3),
        Among(u"viite", -1, 3),
        Among(u"viivad", -1, 3),
        Among(u"võib", -1, 11),
        Among(u"võid", -1, 11),
        Among(u"võida", 274, 11),
        Among(u"võidakse", 275, 11),
        Among(u"võidi", 274, 11),
        Among(u"võiks", -1, 11),
        Among(u"võiksid", 278, 11),
        Among(u"võiksime", 278, 11),
        Among(u"võiksin", 278, 11),
        Among(u"võiksite", 278, 11),
        Among(u"võima", -1, 11),
        Among(u"võimata", 283, 11),
        Among(u"võime", -1, 11),
        Among(u"võin", -1, 11),
        Among(u"võis", -1, 11),
        Among(u"võite", -1, 11),
        Among(u"võivad", -1, 11)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass


class lab2(BaseException): pass


class lab3(BaseException): pass


class lab4(BaseException): pass

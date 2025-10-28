#-*- coding: utf-8 -*-
# Generated from lithuanian.sbl by Snowball 3.0.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class LithuanianStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from lithuanian.sbl by Snowball 3.0.1 - https://snowballstem.org/
    '''

    g_v = {u"a", u"e", u"i", u"y", u"o", u"u", u"ą", u"ę", u"į", u"ų", u"ė", u"ū"}

    I_p1 = 0

    def __r_step1(self):
        if self.cursor < self.I_p1:
            return False
        v_2 = self.limit_backward
        self.limit_backward = self.I_p1
        self.ket = self.cursor
        if self.find_among_b(LithuanianStemmer.a_0) == 0:
            self.limit_backward = v_2
            return False
        self.bra = self.cursor
        self.limit_backward = v_2
        if not self.slice_del():
            return False

        return True

    def __r_step2(self):
        while True:
            v_1 = self.limit - self.cursor
            try:
                if self.cursor < self.I_p1:
                    raise lab0()
                v_3 = self.limit_backward
                self.limit_backward = self.I_p1
                self.ket = self.cursor
                if self.find_among_b(LithuanianStemmer.a_1) == 0:
                    self.limit_backward = v_3
                    raise lab0()
                self.bra = self.cursor
                self.limit_backward = v_3
                if not self.slice_del():
                    return False

                continue
            except lab0: pass
            self.cursor = self.limit - v_1
            break
        return True

    def __r_fix_conflicts(self):
        self.ket = self.cursor
        among_var = self.find_among_b(LithuanianStemmer.a_2)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if not self.slice_from(u"aitė"):
                return False
        elif among_var == 2:
            if not self.slice_from(u"uotė"):
                return False
        elif among_var == 3:
            if not self.slice_from(u"ėjimas"):
                return False
        elif among_var == 4:
            if not self.slice_from(u"esys"):
                return False
        elif among_var == 5:
            if not self.slice_from(u"asys"):
                return False
        elif among_var == 6:
            if not self.slice_from(u"avimas"):
                return False
        elif among_var == 7:
            if not self.slice_from(u"ojimas"):
                return False
        else:
            if not self.slice_from(u"okatė"):
                return False
        return True

    def __r_fix_chdz(self):
        self.ket = self.cursor
        among_var = self.find_among_b(LithuanianStemmer.a_3)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if not self.slice_from(u"t"):
                return False
        else:
            if not self.slice_from(u"d"):
                return False
        return True

    def __r_fix_gd(self):
        self.ket = self.cursor
        if self.find_among_b(LithuanianStemmer.a_4) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_from(u"g"):
            return False
        return True

    def _stem(self):
        self.I_p1 = self.limit
        v_1 = self.cursor
        try:
            v_2 = self.cursor
            try:
                v_3 = self.cursor
                if not self.eq_s(u"a"):
                    self.cursor = v_2
                    raise lab1()
                self.cursor = v_3
                if len(self.current) <= 6:
                    self.cursor = v_2
                    raise lab1()
                if self.cursor >= self.limit:
                    self.cursor = v_2
                    raise lab1()
                self.cursor += 1
            except lab1: pass
            if not self.go_out_grouping(LithuanianStemmer.g_v):
                raise lab0()
            self.cursor += 1
            if not self.go_in_grouping(LithuanianStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_p1 = self.cursor
        except lab0: pass
        self.cursor = v_1
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_4 = self.limit - self.cursor
        self.__r_fix_conflicts()
        self.cursor = self.limit - v_4
        v_5 = self.limit - self.cursor
        self.__r_step1()
        self.cursor = self.limit - v_5
        v_6 = self.limit - self.cursor
        self.__r_fix_chdz()
        self.cursor = self.limit - v_6
        v_7 = self.limit - self.cursor
        self.__r_step2()
        self.cursor = self.limit - v_7
        v_8 = self.limit - self.cursor
        self.__r_fix_chdz()
        self.cursor = self.limit - v_8
        v_9 = self.limit - self.cursor
        self.__r_fix_gd()
        self.cursor = self.limit - v_9
        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among(u"a", -1, -1),
        Among(u"ia", 0, -1),
        Among(u"eria", 1, -1),
        Among(u"osna", 0, -1),
        Among(u"iosna", 3, -1),
        Among(u"uosna", 3, -1),
        Among(u"iuosna", 5, -1),
        Among(u"ysna", 0, -1),
        Among(u"ėsna", 0, -1),
        Among(u"e", -1, -1),
        Among(u"ie", 9, -1),
        Among(u"enie", 10, -1),
        Among(u"erie", 10, -1),
        Among(u"oje", 9, -1),
        Among(u"ioje", 13, -1),
        Among(u"uje", 9, -1),
        Among(u"iuje", 15, -1),
        Among(u"yje", 9, -1),
        Among(u"enyje", 17, -1),
        Among(u"eryje", 17, -1),
        Among(u"ėje", 9, -1),
        Among(u"ame", 9, -1),
        Among(u"iame", 21, -1),
        Among(u"sime", 9, -1),
        Among(u"ome", 9, -1),
        Among(u"ėme", 9, -1),
        Among(u"tumėme", 25, -1),
        Among(u"ose", 9, -1),
        Among(u"iose", 27, -1),
        Among(u"uose", 27, -1),
        Among(u"iuose", 29, -1),
        Among(u"yse", 9, -1),
        Among(u"enyse", 31, -1),
        Among(u"eryse", 31, -1),
        Among(u"ėse", 9, -1),
        Among(u"ate", 9, -1),
        Among(u"iate", 35, -1),
        Among(u"ite", 9, -1),
        Among(u"kite", 37, -1),
        Among(u"site", 37, -1),
        Among(u"ote", 9, -1),
        Among(u"tute", 9, -1),
        Among(u"ėte", 9, -1),
        Among(u"tumėte", 42, -1),
        Among(u"i", -1, -1),
        Among(u"ai", 44, -1),
        Among(u"iai", 45, -1),
        Among(u"eriai", 46, -1),
        Among(u"ei", 44, -1),
        Among(u"tumei", 48, -1),
        Among(u"ki", 44, -1),
        Among(u"imi", 44, -1),
        Among(u"erimi", 51, -1),
        Among(u"umi", 44, -1),
        Among(u"iumi", 53, -1),
        Among(u"si", 44, -1),
        Among(u"asi", 55, -1),
        Among(u"iasi", 56, -1),
        Among(u"esi", 55, -1),
        Among(u"iesi", 58, -1),
        Among(u"siesi", 59, -1),
        Among(u"isi", 55, -1),
        Among(u"aisi", 61, -1),
        Among(u"eisi", 61, -1),
        Among(u"tumeisi", 63, -1),
        Among(u"uisi", 61, -1),
        Among(u"osi", 55, -1),
        Among(u"ėjosi", 66, -1),
        Among(u"uosi", 66, -1),
        Among(u"iuosi", 68, -1),
        Among(u"siuosi", 69, -1),
        Among(u"usi", 55, -1),
        Among(u"ausi", 71, -1),
        Among(u"čiausi", 72, -1),
        Among(u"ąsi", 55, -1),
        Among(u"ėsi", 55, -1),
        Among(u"ųsi", 55, -1),
        Among(u"tųsi", 76, -1),
        Among(u"ti", 44, -1),
        Among(u"enti", 78, -1),
        Among(u"inti", 78, -1),
        Among(u"oti", 78, -1),
        Among(u"ioti", 81, -1),
        Among(u"uoti", 81, -1),
        Among(u"iuoti", 83, -1),
        Among(u"auti", 78, -1),
        Among(u"iauti", 85, -1),
        Among(u"yti", 78, -1),
        Among(u"ėti", 78, -1),
        Among(u"telėti", 88, -1),
        Among(u"inėti", 88, -1),
        Among(u"terėti", 88, -1),
        Among(u"ui", 44, -1),
        Among(u"iui", 92, -1),
        Among(u"eniui", 93, -1),
        Among(u"oj", -1, -1),
        Among(u"ėj", -1, -1),
        Among(u"k", -1, -1),
        Among(u"am", -1, -1),
        Among(u"iam", 98, -1),
        Among(u"iem", -1, -1),
        Among(u"im", -1, -1),
        Among(u"sim", 101, -1),
        Among(u"om", -1, -1),
        Among(u"tum", -1, -1),
        Among(u"ėm", -1, -1),
        Among(u"tumėm", 105, -1),
        Among(u"an", -1, -1),
        Among(u"on", -1, -1),
        Among(u"ion", 108, -1),
        Among(u"un", -1, -1),
        Among(u"iun", 110, -1),
        Among(u"ėn", -1, -1),
        Among(u"o", -1, -1),
        Among(u"io", 113, -1),
        Among(u"enio", 114, -1),
        Among(u"ėjo", 113, -1),
        Among(u"uo", 113, -1),
        Among(u"s", -1, -1),
        Among(u"as", 118, -1),
        Among(u"ias", 119, -1),
        Among(u"es", 118, -1),
        Among(u"ies", 121, -1),
        Among(u"is", 118, -1),
        Among(u"ais", 123, -1),
        Among(u"iais", 124, -1),
        Among(u"tumeis", 123, -1),
        Among(u"imis", 123, -1),
        Among(u"enimis", 127, -1),
        Among(u"omis", 123, -1),
        Among(u"iomis", 129, -1),
        Among(u"umis", 123, -1),
        Among(u"ėmis", 123, -1),
        Among(u"enis", 123, -1),
        Among(u"asis", 123, -1),
        Among(u"ysis", 123, -1),
        Among(u"ams", 118, -1),
        Among(u"iams", 136, -1),
        Among(u"iems", 118, -1),
        Among(u"ims", 118, -1),
        Among(u"enims", 139, -1),
        Among(u"erims", 139, -1),
        Among(u"oms", 118, -1),
        Among(u"ioms", 142, -1),
        Among(u"ums", 118, -1),
        Among(u"ėms", 118, -1),
        Among(u"ens", 118, -1),
        Among(u"os", 118, -1),
        Among(u"ios", 147, -1),
        Among(u"uos", 147, -1),
        Among(u"iuos", 149, -1),
        Among(u"ers", 118, -1),
        Among(u"us", 118, -1),
        Among(u"aus", 152, -1),
        Among(u"iaus", 153, -1),
        Among(u"ius", 152, -1),
        Among(u"ys", 118, -1),
        Among(u"enys", 156, -1),
        Among(u"erys", 156, -1),
        Among(u"ąs", 118, -1),
        Among(u"iąs", 159, -1),
        Among(u"ės", 118, -1),
        Among(u"amės", 161, -1),
        Among(u"iamės", 162, -1),
        Among(u"imės", 161, -1),
        Among(u"kimės", 164, -1),
        Among(u"simės", 164, -1),
        Among(u"omės", 161, -1),
        Among(u"ėmės", 161, -1),
        Among(u"tumėmės", 168, -1),
        Among(u"atės", 161, -1),
        Among(u"iatės", 170, -1),
        Among(u"sitės", 161, -1),
        Among(u"otės", 161, -1),
        Among(u"ėtės", 161, -1),
        Among(u"tumėtės", 174, -1),
        Among(u"įs", 118, -1),
        Among(u"ūs", 118, -1),
        Among(u"tųs", 118, -1),
        Among(u"at", -1, -1),
        Among(u"iat", 179, -1),
        Among(u"it", -1, -1),
        Among(u"sit", 181, -1),
        Among(u"ot", -1, -1),
        Among(u"ėt", -1, -1),
        Among(u"tumėt", 184, -1),
        Among(u"u", -1, -1),
        Among(u"au", 186, -1),
        Among(u"iau", 187, -1),
        Among(u"čiau", 188, -1),
        Among(u"iu", 186, -1),
        Among(u"eniu", 190, -1),
        Among(u"siu", 190, -1),
        Among(u"y", -1, -1),
        Among(u"ą", -1, -1),
        Among(u"ią", 194, -1),
        Among(u"ė", -1, -1),
        Among(u"ę", -1, -1),
        Among(u"į", -1, -1),
        Among(u"enį", 198, -1),
        Among(u"erį", 198, -1),
        Among(u"ų", -1, -1),
        Among(u"ių", 201, -1),
        Among(u"erų", 201, -1)
    ]

    a_1 = [
        Among(u"ing", -1, -1),
        Among(u"aj", -1, -1),
        Among(u"iaj", 1, -1),
        Among(u"iej", -1, -1),
        Among(u"oj", -1, -1),
        Among(u"ioj", 4, -1),
        Among(u"uoj", 4, -1),
        Among(u"iuoj", 6, -1),
        Among(u"auj", -1, -1),
        Among(u"ąj", -1, -1),
        Among(u"iąj", 9, -1),
        Among(u"ėj", -1, -1),
        Among(u"ųj", -1, -1),
        Among(u"iųj", 12, -1),
        Among(u"ok", -1, -1),
        Among(u"iok", 14, -1),
        Among(u"iuk", -1, -1),
        Among(u"uliuk", 16, -1),
        Among(u"učiuk", 16, -1),
        Among(u"išk", -1, -1),
        Among(u"iul", -1, -1),
        Among(u"yl", -1, -1),
        Among(u"ėl", -1, -1),
        Among(u"am", -1, -1),
        Among(u"dam", 23, -1),
        Among(u"jam", 23, -1),
        Among(u"zgan", -1, -1),
        Among(u"ain", -1, -1),
        Among(u"esn", -1, -1),
        Among(u"op", -1, -1),
        Among(u"iop", 29, -1),
        Among(u"ias", -1, -1),
        Among(u"ies", -1, -1),
        Among(u"ais", -1, -1),
        Among(u"iais", 33, -1),
        Among(u"os", -1, -1),
        Among(u"ios", 35, -1),
        Among(u"uos", 35, -1),
        Among(u"iuos", 37, -1),
        Among(u"aus", -1, -1),
        Among(u"iaus", 39, -1),
        Among(u"ąs", -1, -1),
        Among(u"iąs", 41, -1),
        Among(u"ęs", -1, -1),
        Among(u"utėait", -1, -1),
        Among(u"ant", -1, -1),
        Among(u"iant", 45, -1),
        Among(u"siant", 46, -1),
        Among(u"int", -1, -1),
        Among(u"ot", -1, -1),
        Among(u"uot", 49, -1),
        Among(u"iuot", 50, -1),
        Among(u"yt", -1, -1),
        Among(u"ėt", -1, -1),
        Among(u"ykšt", -1, -1),
        Among(u"iau", -1, -1),
        Among(u"dav", -1, -1),
        Among(u"sv", -1, -1),
        Among(u"šv", -1, -1),
        Among(u"ykšč", -1, -1),
        Among(u"ę", -1, -1),
        Among(u"ėję", 60, -1)
    ]

    a_2 = [
        Among(u"ojime", -1, 7),
        Among(u"ėjime", -1, 3),
        Among(u"avime", -1, 6),
        Among(u"okate", -1, 8),
        Among(u"aite", -1, 1),
        Among(u"uote", -1, 2),
        Among(u"asius", -1, 5),
        Among(u"okatės", -1, 8),
        Among(u"aitės", -1, 1),
        Among(u"uotės", -1, 2),
        Among(u"esiu", -1, 4)
    ]

    a_3 = [
        Among(u"č", -1, 1),
        Among(u"dž", -1, 2)
    ]

    a_4 = [
        Among(u"gd", -1, 1)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass

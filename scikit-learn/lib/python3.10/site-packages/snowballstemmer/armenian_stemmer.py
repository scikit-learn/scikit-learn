#-*- coding: utf-8 -*-
# Generated from armenian.sbl by Snowball 3.0.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class ArmenianStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from armenian.sbl by Snowball 3.0.1 - https://snowballstem.org/
    '''

    g_v = {u"ա", u"է", u"ի", u"օ", u"ւ", u"ե", u"ո", u"ը"}

    I_p2 = 0
    I_pV = 0

    def __r_mark_regions(self):
        self.I_pV = self.limit
        self.I_p2 = self.limit
        v_1 = self.cursor
        try:
            if not self.go_out_grouping(ArmenianStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_pV = self.cursor
            if not self.go_in_grouping(ArmenianStemmer.g_v):
                raise lab0()
            self.cursor += 1
            if not self.go_out_grouping(ArmenianStemmer.g_v):
                raise lab0()
            self.cursor += 1
            if not self.go_in_grouping(ArmenianStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_p2 = self.cursor
        except lab0: pass
        self.cursor = v_1
        return True

    def __r_R2(self):
        return self.I_p2 <= self.cursor

    def __r_adjective(self):
        self.ket = self.cursor
        if self.find_among_b(ArmenianStemmer.a_0) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        return True

    def __r_verb(self):
        self.ket = self.cursor
        if self.find_among_b(ArmenianStemmer.a_1) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        return True

    def __r_noun(self):
        self.ket = self.cursor
        if self.find_among_b(ArmenianStemmer.a_2) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        return True

    def __r_ending(self):
        self.ket = self.cursor
        if self.find_among_b(ArmenianStemmer.a_3) == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R2():
            return False
        if not self.slice_del():
            return False

        return True

    def _stem(self):
        self.__r_mark_regions()
        self.limit_backward = self.cursor
        self.cursor = self.limit
        if self.cursor < self.I_pV:
            return False
        v_3 = self.limit_backward
        self.limit_backward = self.I_pV
        v_4 = self.limit - self.cursor
        self.__r_ending()
        self.cursor = self.limit - v_4
        v_5 = self.limit - self.cursor
        self.__r_verb()
        self.cursor = self.limit - v_5
        v_6 = self.limit - self.cursor
        self.__r_adjective()
        self.cursor = self.limit - v_6
        v_7 = self.limit - self.cursor
        self.__r_noun()
        self.cursor = self.limit - v_7
        self.limit_backward = v_3
        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among(u"րորդ", -1, 1),
        Among(u"երորդ", 0, 1),
        Among(u"ալի", -1, 1),
        Among(u"ակի", -1, 1),
        Among(u"որակ", -1, 1),
        Among(u"եղ", -1, 1),
        Among(u"ական", -1, 1),
        Among(u"արան", -1, 1),
        Among(u"են", -1, 1),
        Among(u"եկեն", 8, 1),
        Among(u"երեն", 8, 1),
        Among(u"որէն", -1, 1),
        Among(u"ին", -1, 1),
        Among(u"գին", 12, 1),
        Among(u"ովին", 12, 1),
        Among(u"լայն", -1, 1),
        Among(u"վուն", -1, 1),
        Among(u"պես", -1, 1),
        Among(u"իվ", -1, 1),
        Among(u"ատ", -1, 1),
        Among(u"ավետ", -1, 1),
        Among(u"կոտ", -1, 1),
        Among(u"բար", -1, 1)
    ]

    a_1 = [
        Among(u"ա", -1, 1),
        Among(u"ացա", 0, 1),
        Among(u"եցա", 0, 1),
        Among(u"վե", -1, 1),
        Among(u"ացրի", -1, 1),
        Among(u"ացի", -1, 1),
        Among(u"եցի", -1, 1),
        Among(u"վեցի", 6, 1),
        Among(u"ալ", -1, 1),
        Among(u"ըալ", 8, 1),
        Among(u"անալ", 8, 1),
        Among(u"ենալ", 8, 1),
        Among(u"ացնալ", 8, 1),
        Among(u"ել", -1, 1),
        Among(u"ըել", 13, 1),
        Among(u"նել", 13, 1),
        Among(u"ցնել", 15, 1),
        Among(u"եցնել", 16, 1),
        Among(u"չել", 13, 1),
        Among(u"վել", 13, 1),
        Among(u"ացվել", 19, 1),
        Among(u"եցվել", 19, 1),
        Among(u"տել", 13, 1),
        Among(u"ատել", 22, 1),
        Among(u"ոտել", 22, 1),
        Among(u"կոտել", 24, 1),
        Among(u"ված", -1, 1),
        Among(u"ում", -1, 1),
        Among(u"վում", 27, 1),
        Among(u"ան", -1, 1),
        Among(u"ցան", 29, 1),
        Among(u"ացան", 30, 1),
        Among(u"ացրին", -1, 1),
        Among(u"ացին", -1, 1),
        Among(u"եցին", -1, 1),
        Among(u"վեցին", 34, 1),
        Among(u"ալիս", -1, 1),
        Among(u"ելիս", -1, 1),
        Among(u"ավ", -1, 1),
        Among(u"ացավ", 38, 1),
        Among(u"եցավ", 38, 1),
        Among(u"ալով", -1, 1),
        Among(u"ելով", -1, 1),
        Among(u"ար", -1, 1),
        Among(u"ացար", 43, 1),
        Among(u"եցար", 43, 1),
        Among(u"ացրիր", -1, 1),
        Among(u"ացիր", -1, 1),
        Among(u"եցիր", -1, 1),
        Among(u"վեցիր", 48, 1),
        Among(u"աց", -1, 1),
        Among(u"եց", -1, 1),
        Among(u"ացրեց", 51, 1),
        Among(u"ալուց", -1, 1),
        Among(u"ելուց", -1, 1),
        Among(u"ալու", -1, 1),
        Among(u"ելու", -1, 1),
        Among(u"աք", -1, 1),
        Among(u"ցաք", 57, 1),
        Among(u"ացաք", 58, 1),
        Among(u"ացրիք", -1, 1),
        Among(u"ացիք", -1, 1),
        Among(u"եցիք", -1, 1),
        Among(u"վեցիք", 62, 1),
        Among(u"անք", -1, 1),
        Among(u"ցանք", 64, 1),
        Among(u"ացանք", 65, 1),
        Among(u"ացրինք", -1, 1),
        Among(u"ացինք", -1, 1),
        Among(u"եցինք", -1, 1),
        Among(u"վեցինք", 69, 1)
    ]

    a_2 = [
        Among(u"որդ", -1, 1),
        Among(u"ույթ", -1, 1),
        Among(u"ուհի", -1, 1),
        Among(u"ցի", -1, 1),
        Among(u"իլ", -1, 1),
        Among(u"ակ", -1, 1),
        Among(u"յակ", 5, 1),
        Among(u"անակ", 5, 1),
        Among(u"իկ", -1, 1),
        Among(u"ուկ", -1, 1),
        Among(u"ան", -1, 1),
        Among(u"պան", 10, 1),
        Among(u"ստան", 10, 1),
        Among(u"արան", 10, 1),
        Among(u"եղէն", -1, 1),
        Among(u"յուն", -1, 1),
        Among(u"ություն", 15, 1),
        Among(u"ածո", -1, 1),
        Among(u"իչ", -1, 1),
        Among(u"ուս", -1, 1),
        Among(u"ուստ", -1, 1),
        Among(u"գար", -1, 1),
        Among(u"վոր", -1, 1),
        Among(u"ավոր", 22, 1),
        Among(u"ոց", -1, 1),
        Among(u"անօց", -1, 1),
        Among(u"ու", -1, 1),
        Among(u"ք", -1, 1),
        Among(u"չեք", 27, 1),
        Among(u"իք", 27, 1),
        Among(u"ալիք", 29, 1),
        Among(u"անիք", 29, 1),
        Among(u"վածք", 27, 1),
        Among(u"ույք", 27, 1),
        Among(u"ենք", 27, 1),
        Among(u"ոնք", 27, 1),
        Among(u"ունք", 27, 1),
        Among(u"մունք", 36, 1),
        Among(u"իչք", 27, 1),
        Among(u"արք", 27, 1)
    ]

    a_3 = [
        Among(u"սա", -1, 1),
        Among(u"վա", -1, 1),
        Among(u"ամբ", -1, 1),
        Among(u"դ", -1, 1),
        Among(u"անդ", 3, 1),
        Among(u"ությանդ", 4, 1),
        Among(u"վանդ", 4, 1),
        Among(u"ոջդ", 3, 1),
        Among(u"երդ", 3, 1),
        Among(u"ներդ", 8, 1),
        Among(u"ուդ", 3, 1),
        Among(u"ը", -1, 1),
        Among(u"անը", 11, 1),
        Among(u"ությանը", 12, 1),
        Among(u"վանը", 12, 1),
        Among(u"ոջը", 11, 1),
        Among(u"երը", 11, 1),
        Among(u"ները", 16, 1),
        Among(u"ի", -1, 1),
        Among(u"վի", 18, 1),
        Among(u"երի", 18, 1),
        Among(u"ների", 20, 1),
        Among(u"անում", -1, 1),
        Among(u"երում", -1, 1),
        Among(u"ներում", 23, 1),
        Among(u"ն", -1, 1),
        Among(u"ան", 25, 1),
        Among(u"ության", 26, 1),
        Among(u"վան", 26, 1),
        Among(u"ին", 25, 1),
        Among(u"երին", 29, 1),
        Among(u"ներին", 30, 1),
        Among(u"ությանն", 25, 1),
        Among(u"երն", 25, 1),
        Among(u"ներն", 33, 1),
        Among(u"ուն", 25, 1),
        Among(u"ոջ", -1, 1),
        Among(u"ությանս", -1, 1),
        Among(u"վանս", -1, 1),
        Among(u"ոջս", -1, 1),
        Among(u"ով", -1, 1),
        Among(u"անով", 40, 1),
        Among(u"վով", 40, 1),
        Among(u"երով", 40, 1),
        Among(u"ներով", 43, 1),
        Among(u"եր", -1, 1),
        Among(u"ներ", 45, 1),
        Among(u"ց", -1, 1),
        Among(u"ից", 47, 1),
        Among(u"վանից", 48, 1),
        Among(u"ոջից", 48, 1),
        Among(u"վից", 48, 1),
        Among(u"երից", 48, 1),
        Among(u"ներից", 52, 1),
        Among(u"ցից", 48, 1),
        Among(u"ոց", 47, 1),
        Among(u"ուց", 47, 1)
    ]


class lab0(BaseException): pass

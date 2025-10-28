#-*- coding: utf-8 -*-
# Generated from catalan.sbl by Snowball 3.0.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class CatalanStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from catalan.sbl by Snowball 3.0.1 - https://snowballstem.org/
    '''

    g_v = {u"a", u"e", u"i", u"o", u"u", u"á", u"à", u"é", u"è", u"í", u"ï", u"ó", u"ò", u"ú", u"ü"}

    I_p2 = 0
    I_p1 = 0

    def __r_mark_regions(self):
        self.I_p1 = self.limit
        self.I_p2 = self.limit
        v_1 = self.cursor
        try:
            if not self.go_out_grouping(CatalanStemmer.g_v):
                raise lab0()
            self.cursor += 1
            if not self.go_in_grouping(CatalanStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_p1 = self.cursor
            if not self.go_out_grouping(CatalanStemmer.g_v):
                raise lab0()
            self.cursor += 1
            if not self.go_in_grouping(CatalanStemmer.g_v):
                raise lab0()
            self.cursor += 1
            self.I_p2 = self.cursor
        except lab0: pass
        self.cursor = v_1
        return True

    def __r_cleaning(self):
        while True:
            v_1 = self.cursor
            try:
                self.bra = self.cursor
                among_var = self.find_among(CatalanStemmer.a_0)
                self.ket = self.cursor
                if among_var == 1:
                    if not self.slice_from(u"a"):
                        return False
                elif among_var == 2:
                    if not self.slice_from(u"e"):
                        return False
                elif among_var == 3:
                    if not self.slice_from(u"i"):
                        return False
                elif among_var == 4:
                    if not self.slice_from(u"o"):
                        return False
                elif among_var == 5:
                    if not self.slice_from(u"u"):
                        return False
                elif among_var == 6:
                    if not self.slice_from(u"."):
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

    def __r_attached_pronoun(self):
        self.ket = self.cursor
        if self.find_among_b(CatalanStemmer.a_1) == 0:
            return False
        self.bra = self.cursor
        if not self.__r_R1():
            return False
        if not self.slice_del():
            return False

        return True

    def __r_standard_suffix(self):
        self.ket = self.cursor
        among_var = self.find_among_b(CatalanStemmer.a_2)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if not self.__r_R1():
                return False
            if not self.slice_del():
                return False

        elif among_var == 2:
            if not self.__r_R2():
                return False
            if not self.slice_del():
                return False

        elif among_var == 3:
            if not self.__r_R2():
                return False
            if not self.slice_from(u"log"):
                return False
        elif among_var == 4:
            if not self.__r_R2():
                return False
            if not self.slice_from(u"ic"):
                return False
        else:
            if not self.__r_R1():
                return False
            if not self.slice_from(u"c"):
                return False
        return True

    def __r_verb_suffix(self):
        self.ket = self.cursor
        among_var = self.find_among_b(CatalanStemmer.a_3)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if not self.__r_R1():
                return False
            if not self.slice_del():
                return False

        else:
            if not self.__r_R2():
                return False
            if not self.slice_del():
                return False

        return True

    def __r_residual_suffix(self):
        self.ket = self.cursor
        among_var = self.find_among_b(CatalanStemmer.a_4)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if not self.__r_R1():
                return False
            if not self.slice_del():
                return False

        else:
            if not self.__r_R1():
                return False
            if not self.slice_from(u"ic"):
                return False
        return True

    def _stem(self):
        self.__r_mark_regions()
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_2 = self.limit - self.cursor
        self.__r_attached_pronoun()
        self.cursor = self.limit - v_2
        v_3 = self.limit - self.cursor
        try:
            try:
                v_4 = self.limit - self.cursor
                try:
                    if not self.__r_standard_suffix():
                        raise lab2()
                    raise lab1()
                except lab2: pass
                self.cursor = self.limit - v_4
                if not self.__r_verb_suffix():
                    raise lab0()
            except lab1: pass
        except lab0: pass
        self.cursor = self.limit - v_3
        v_5 = self.limit - self.cursor
        self.__r_residual_suffix()
        self.cursor = self.limit - v_5
        self.cursor = self.limit_backward
        v_6 = self.cursor
        self.__r_cleaning()
        self.cursor = v_6
        return True

    a_0 = [
        Among(u"", -1, 7),
        Among(u"·", 0, 6),
        Among(u"à", 0, 1),
        Among(u"á", 0, 1),
        Among(u"è", 0, 2),
        Among(u"é", 0, 2),
        Among(u"ì", 0, 3),
        Among(u"í", 0, 3),
        Among(u"ï", 0, 3),
        Among(u"ò", 0, 4),
        Among(u"ó", 0, 4),
        Among(u"ú", 0, 5),
        Among(u"ü", 0, 5)
    ]

    a_1 = [
        Among(u"la", -1, 1),
        Among(u"-la", 0, 1),
        Among(u"sela", 0, 1),
        Among(u"le", -1, 1),
        Among(u"me", -1, 1),
        Among(u"-me", 4, 1),
        Among(u"se", -1, 1),
        Among(u"-te", -1, 1),
        Among(u"hi", -1, 1),
        Among(u"'hi", 8, 1),
        Among(u"li", -1, 1),
        Among(u"-li", 10, 1),
        Among(u"'l", -1, 1),
        Among(u"'m", -1, 1),
        Among(u"-m", -1, 1),
        Among(u"'n", -1, 1),
        Among(u"-n", -1, 1),
        Among(u"ho", -1, 1),
        Among(u"'ho", 17, 1),
        Among(u"lo", -1, 1),
        Among(u"selo", 19, 1),
        Among(u"'s", -1, 1),
        Among(u"las", -1, 1),
        Among(u"selas", 22, 1),
        Among(u"les", -1, 1),
        Among(u"-les", 24, 1),
        Among(u"'ls", -1, 1),
        Among(u"-ls", -1, 1),
        Among(u"'ns", -1, 1),
        Among(u"-ns", -1, 1),
        Among(u"ens", -1, 1),
        Among(u"los", -1, 1),
        Among(u"selos", 31, 1),
        Among(u"nos", -1, 1),
        Among(u"-nos", 33, 1),
        Among(u"vos", -1, 1),
        Among(u"us", -1, 1),
        Among(u"-us", 36, 1),
        Among(u"'t", -1, 1)
    ]

    a_2 = [
        Among(u"ica", -1, 4),
        Among(u"lógica", 0, 3),
        Among(u"enca", -1, 1),
        Among(u"ada", -1, 2),
        Among(u"ancia", -1, 1),
        Among(u"encia", -1, 1),
        Among(u"ència", -1, 1),
        Among(u"ícia", -1, 1),
        Among(u"logia", -1, 3),
        Among(u"inia", -1, 1),
        Among(u"íinia", 9, 1),
        Among(u"eria", -1, 1),
        Among(u"ària", -1, 1),
        Among(u"atòria", -1, 1),
        Among(u"alla", -1, 1),
        Among(u"ella", -1, 1),
        Among(u"ívola", -1, 1),
        Among(u"ima", -1, 1),
        Among(u"íssima", 17, 1),
        Among(u"quíssima", 18, 5),
        Among(u"ana", -1, 1),
        Among(u"ina", -1, 1),
        Among(u"era", -1, 1),
        Among(u"sfera", 22, 1),
        Among(u"ora", -1, 1),
        Among(u"dora", 24, 1),
        Among(u"adora", 25, 1),
        Among(u"adura", -1, 1),
        Among(u"esa", -1, 1),
        Among(u"osa", -1, 1),
        Among(u"assa", -1, 1),
        Among(u"essa", -1, 1),
        Among(u"issa", -1, 1),
        Among(u"eta", -1, 1),
        Among(u"ita", -1, 1),
        Among(u"ota", -1, 1),
        Among(u"ista", -1, 1),
        Among(u"ialista", 36, 1),
        Among(u"ionista", 36, 1),
        Among(u"iva", -1, 1),
        Among(u"ativa", 39, 1),
        Among(u"nça", -1, 1),
        Among(u"logía", -1, 3),
        Among(u"ic", -1, 4),
        Among(u"ístic", 43, 1),
        Among(u"enc", -1, 1),
        Among(u"esc", -1, 1),
        Among(u"ud", -1, 1),
        Among(u"atge", -1, 1),
        Among(u"ble", -1, 1),
        Among(u"able", 49, 1),
        Among(u"ible", 49, 1),
        Among(u"isme", -1, 1),
        Among(u"ialisme", 52, 1),
        Among(u"ionisme", 52, 1),
        Among(u"ivisme", 52, 1),
        Among(u"aire", -1, 1),
        Among(u"icte", -1, 1),
        Among(u"iste", -1, 1),
        Among(u"ici", -1, 1),
        Among(u"íci", -1, 1),
        Among(u"logi", -1, 3),
        Among(u"ari", -1, 1),
        Among(u"tori", -1, 1),
        Among(u"al", -1, 1),
        Among(u"il", -1, 1),
        Among(u"all", -1, 1),
        Among(u"ell", -1, 1),
        Among(u"ívol", -1, 1),
        Among(u"isam", -1, 1),
        Among(u"issem", -1, 1),
        Among(u"ìssem", -1, 1),
        Among(u"íssem", -1, 1),
        Among(u"íssim", -1, 1),
        Among(u"quíssim", 73, 5),
        Among(u"amen", -1, 1),
        Among(u"ìssin", -1, 1),
        Among(u"ar", -1, 1),
        Among(u"ificar", 77, 1),
        Among(u"egar", 77, 1),
        Among(u"ejar", 77, 1),
        Among(u"itar", 77, 1),
        Among(u"itzar", 77, 1),
        Among(u"fer", -1, 1),
        Among(u"or", -1, 1),
        Among(u"dor", 84, 1),
        Among(u"dur", -1, 1),
        Among(u"doras", -1, 1),
        Among(u"ics", -1, 4),
        Among(u"lógics", 88, 3),
        Among(u"uds", -1, 1),
        Among(u"nces", -1, 1),
        Among(u"ades", -1, 2),
        Among(u"ancies", -1, 1),
        Among(u"encies", -1, 1),
        Among(u"ències", -1, 1),
        Among(u"ícies", -1, 1),
        Among(u"logies", -1, 3),
        Among(u"inies", -1, 1),
        Among(u"ínies", -1, 1),
        Among(u"eries", -1, 1),
        Among(u"àries", -1, 1),
        Among(u"atòries", -1, 1),
        Among(u"bles", -1, 1),
        Among(u"ables", 103, 1),
        Among(u"ibles", 103, 1),
        Among(u"imes", -1, 1),
        Among(u"íssimes", 106, 1),
        Among(u"quíssimes", 107, 5),
        Among(u"formes", -1, 1),
        Among(u"ismes", -1, 1),
        Among(u"ialismes", 110, 1),
        Among(u"ines", -1, 1),
        Among(u"eres", -1, 1),
        Among(u"ores", -1, 1),
        Among(u"dores", 114, 1),
        Among(u"idores", 115, 1),
        Among(u"dures", -1, 1),
        Among(u"eses", -1, 1),
        Among(u"oses", -1, 1),
        Among(u"asses", -1, 1),
        Among(u"ictes", -1, 1),
        Among(u"ites", -1, 1),
        Among(u"otes", -1, 1),
        Among(u"istes", -1, 1),
        Among(u"ialistes", 124, 1),
        Among(u"ionistes", 124, 1),
        Among(u"iques", -1, 4),
        Among(u"lógiques", 127, 3),
        Among(u"ives", -1, 1),
        Among(u"atives", 129, 1),
        Among(u"logíes", -1, 3),
        Among(u"allengües", -1, 1),
        Among(u"icis", -1, 1),
        Among(u"ícis", -1, 1),
        Among(u"logis", -1, 3),
        Among(u"aris", -1, 1),
        Among(u"toris", -1, 1),
        Among(u"ls", -1, 1),
        Among(u"als", 138, 1),
        Among(u"ells", 138, 1),
        Among(u"ims", -1, 1),
        Among(u"íssims", 141, 1),
        Among(u"quíssims", 142, 5),
        Among(u"ions", -1, 1),
        Among(u"cions", 144, 1),
        Among(u"acions", 145, 2),
        Among(u"esos", -1, 1),
        Among(u"osos", -1, 1),
        Among(u"assos", -1, 1),
        Among(u"issos", -1, 1),
        Among(u"ers", -1, 1),
        Among(u"ors", -1, 1),
        Among(u"dors", 152, 1),
        Among(u"adors", 153, 1),
        Among(u"idors", 153, 1),
        Among(u"ats", -1, 1),
        Among(u"itats", 156, 1),
        Among(u"bilitats", 157, 1),
        Among(u"ivitats", 157, 1),
        Among(u"ativitats", 159, 1),
        Among(u"ïtats", 156, 1),
        Among(u"ets", -1, 1),
        Among(u"ants", -1, 1),
        Among(u"ents", -1, 1),
        Among(u"ments", 164, 1),
        Among(u"aments", 165, 1),
        Among(u"ots", -1, 1),
        Among(u"uts", -1, 1),
        Among(u"ius", -1, 1),
        Among(u"trius", 169, 1),
        Among(u"atius", 169, 1),
        Among(u"ès", -1, 1),
        Among(u"és", -1, 1),
        Among(u"ís", -1, 1),
        Among(u"dís", 174, 1),
        Among(u"ós", -1, 1),
        Among(u"itat", -1, 1),
        Among(u"bilitat", 177, 1),
        Among(u"ivitat", 177, 1),
        Among(u"ativitat", 179, 1),
        Among(u"ïtat", -1, 1),
        Among(u"et", -1, 1),
        Among(u"ant", -1, 1),
        Among(u"ent", -1, 1),
        Among(u"ient", 184, 1),
        Among(u"ment", 184, 1),
        Among(u"ament", 186, 1),
        Among(u"isament", 187, 1),
        Among(u"ot", -1, 1),
        Among(u"isseu", -1, 1),
        Among(u"ìsseu", -1, 1),
        Among(u"ísseu", -1, 1),
        Among(u"triu", -1, 1),
        Among(u"íssiu", -1, 1),
        Among(u"atiu", -1, 1),
        Among(u"ó", -1, 1),
        Among(u"ió", 196, 1),
        Among(u"ció", 197, 1),
        Among(u"ació", 198, 1)
    ]

    a_3 = [
        Among(u"aba", -1, 1),
        Among(u"esca", -1, 1),
        Among(u"isca", -1, 1),
        Among(u"ïsca", -1, 1),
        Among(u"ada", -1, 1),
        Among(u"ida", -1, 1),
        Among(u"uda", -1, 1),
        Among(u"ïda", -1, 1),
        Among(u"ia", -1, 1),
        Among(u"aria", 8, 1),
        Among(u"iria", 8, 1),
        Among(u"ara", -1, 1),
        Among(u"iera", -1, 1),
        Among(u"ira", -1, 1),
        Among(u"adora", -1, 1),
        Among(u"ïra", -1, 1),
        Among(u"ava", -1, 1),
        Among(u"ixa", -1, 1),
        Among(u"itza", -1, 1),
        Among(u"ía", -1, 1),
        Among(u"aría", 19, 1),
        Among(u"ería", 19, 1),
        Among(u"iría", 19, 1),
        Among(u"ïa", -1, 1),
        Among(u"isc", -1, 1),
        Among(u"ïsc", -1, 1),
        Among(u"ad", -1, 1),
        Among(u"ed", -1, 1),
        Among(u"id", -1, 1),
        Among(u"ie", -1, 1),
        Among(u"re", -1, 1),
        Among(u"dre", 30, 1),
        Among(u"ase", -1, 1),
        Among(u"iese", -1, 1),
        Among(u"aste", -1, 1),
        Among(u"iste", -1, 1),
        Among(u"ii", -1, 1),
        Among(u"ini", -1, 1),
        Among(u"esqui", -1, 1),
        Among(u"eixi", -1, 1),
        Among(u"itzi", -1, 1),
        Among(u"am", -1, 1),
        Among(u"em", -1, 1),
        Among(u"arem", 42, 1),
        Among(u"irem", 42, 1),
        Among(u"àrem", 42, 1),
        Among(u"írem", 42, 1),
        Among(u"àssem", 42, 1),
        Among(u"éssem", 42, 1),
        Among(u"iguem", 42, 1),
        Among(u"ïguem", 42, 1),
        Among(u"avem", 42, 1),
        Among(u"àvem", 42, 1),
        Among(u"ávem", 42, 1),
        Among(u"irìem", 42, 1),
        Among(u"íem", 42, 1),
        Among(u"aríem", 55, 1),
        Among(u"iríem", 55, 1),
        Among(u"assim", -1, 1),
        Among(u"essim", -1, 1),
        Among(u"issim", -1, 1),
        Among(u"àssim", -1, 1),
        Among(u"èssim", -1, 1),
        Among(u"éssim", -1, 1),
        Among(u"íssim", -1, 1),
        Among(u"ïm", -1, 1),
        Among(u"an", -1, 1),
        Among(u"aban", 66, 1),
        Among(u"arian", 66, 1),
        Among(u"aran", 66, 1),
        Among(u"ieran", 66, 1),
        Among(u"iran", 66, 1),
        Among(u"ían", 66, 1),
        Among(u"arían", 72, 1),
        Among(u"erían", 72, 1),
        Among(u"irían", 72, 1),
        Among(u"en", -1, 1),
        Among(u"ien", 76, 1),
        Among(u"arien", 77, 1),
        Among(u"irien", 77, 1),
        Among(u"aren", 76, 1),
        Among(u"eren", 76, 1),
        Among(u"iren", 76, 1),
        Among(u"àren", 76, 1),
        Among(u"ïren", 76, 1),
        Among(u"asen", 76, 1),
        Among(u"iesen", 76, 1),
        Among(u"assen", 76, 1),
        Among(u"essen", 76, 1),
        Among(u"issen", 76, 1),
        Among(u"éssen", 76, 1),
        Among(u"ïssen", 76, 1),
        Among(u"esquen", 76, 1),
        Among(u"isquen", 76, 1),
        Among(u"ïsquen", 76, 1),
        Among(u"aven", 76, 1),
        Among(u"ixen", 76, 1),
        Among(u"eixen", 96, 1),
        Among(u"ïxen", 76, 1),
        Among(u"ïen", 76, 1),
        Among(u"in", -1, 1),
        Among(u"inin", 100, 1),
        Among(u"sin", 100, 1),
        Among(u"isin", 102, 1),
        Among(u"assin", 102, 1),
        Among(u"essin", 102, 1),
        Among(u"issin", 102, 1),
        Among(u"ïssin", 102, 1),
        Among(u"esquin", 100, 1),
        Among(u"eixin", 100, 1),
        Among(u"aron", -1, 1),
        Among(u"ieron", -1, 1),
        Among(u"arán", -1, 1),
        Among(u"erán", -1, 1),
        Among(u"irán", -1, 1),
        Among(u"iïn", -1, 1),
        Among(u"ado", -1, 1),
        Among(u"ido", -1, 1),
        Among(u"ando", -1, 2),
        Among(u"iendo", -1, 1),
        Among(u"io", -1, 1),
        Among(u"ixo", -1, 1),
        Among(u"eixo", 121, 1),
        Among(u"ïxo", -1, 1),
        Among(u"itzo", -1, 1),
        Among(u"ar", -1, 1),
        Among(u"tzar", 125, 1),
        Among(u"er", -1, 1),
        Among(u"eixer", 127, 1),
        Among(u"ir", -1, 1),
        Among(u"ador", -1, 1),
        Among(u"as", -1, 1),
        Among(u"abas", 131, 1),
        Among(u"adas", 131, 1),
        Among(u"idas", 131, 1),
        Among(u"aras", 131, 1),
        Among(u"ieras", 131, 1),
        Among(u"ías", 131, 1),
        Among(u"arías", 137, 1),
        Among(u"erías", 137, 1),
        Among(u"irías", 137, 1),
        Among(u"ids", -1, 1),
        Among(u"es", -1, 1),
        Among(u"ades", 142, 1),
        Among(u"ides", 142, 1),
        Among(u"udes", 142, 1),
        Among(u"ïdes", 142, 1),
        Among(u"atges", 142, 1),
        Among(u"ies", 142, 1),
        Among(u"aries", 148, 1),
        Among(u"iries", 148, 1),
        Among(u"ares", 142, 1),
        Among(u"ires", 142, 1),
        Among(u"adores", 142, 1),
        Among(u"ïres", 142, 1),
        Among(u"ases", 142, 1),
        Among(u"ieses", 142, 1),
        Among(u"asses", 142, 1),
        Among(u"esses", 142, 1),
        Among(u"isses", 142, 1),
        Among(u"ïsses", 142, 1),
        Among(u"ques", 142, 1),
        Among(u"esques", 161, 1),
        Among(u"ïsques", 161, 1),
        Among(u"aves", 142, 1),
        Among(u"ixes", 142, 1),
        Among(u"eixes", 165, 1),
        Among(u"ïxes", 142, 1),
        Among(u"ïes", 142, 1),
        Among(u"abais", -1, 1),
        Among(u"arais", -1, 1),
        Among(u"ierais", -1, 1),
        Among(u"íais", -1, 1),
        Among(u"aríais", 172, 1),
        Among(u"eríais", 172, 1),
        Among(u"iríais", 172, 1),
        Among(u"aseis", -1, 1),
        Among(u"ieseis", -1, 1),
        Among(u"asteis", -1, 1),
        Among(u"isteis", -1, 1),
        Among(u"inis", -1, 1),
        Among(u"sis", -1, 1),
        Among(u"isis", 181, 1),
        Among(u"assis", 181, 1),
        Among(u"essis", 181, 1),
        Among(u"issis", 181, 1),
        Among(u"ïssis", 181, 1),
        Among(u"esquis", -1, 1),
        Among(u"eixis", -1, 1),
        Among(u"itzis", -1, 1),
        Among(u"áis", -1, 1),
        Among(u"aréis", -1, 1),
        Among(u"eréis", -1, 1),
        Among(u"iréis", -1, 1),
        Among(u"ams", -1, 1),
        Among(u"ados", -1, 1),
        Among(u"idos", -1, 1),
        Among(u"amos", -1, 1),
        Among(u"ábamos", 197, 1),
        Among(u"áramos", 197, 1),
        Among(u"iéramos", 197, 1),
        Among(u"íamos", 197, 1),
        Among(u"aríamos", 201, 1),
        Among(u"eríamos", 201, 1),
        Among(u"iríamos", 201, 1),
        Among(u"aremos", -1, 1),
        Among(u"eremos", -1, 1),
        Among(u"iremos", -1, 1),
        Among(u"ásemos", -1, 1),
        Among(u"iésemos", -1, 1),
        Among(u"imos", -1, 1),
        Among(u"adors", -1, 1),
        Among(u"ass", -1, 1),
        Among(u"erass", 212, 1),
        Among(u"ess", -1, 1),
        Among(u"ats", -1, 1),
        Among(u"its", -1, 1),
        Among(u"ents", -1, 1),
        Among(u"às", -1, 1),
        Among(u"aràs", 218, 1),
        Among(u"iràs", 218, 1),
        Among(u"arás", -1, 1),
        Among(u"erás", -1, 1),
        Among(u"irás", -1, 1),
        Among(u"és", -1, 1),
        Among(u"arés", 224, 1),
        Among(u"ís", -1, 1),
        Among(u"iïs", -1, 1),
        Among(u"at", -1, 1),
        Among(u"it", -1, 1),
        Among(u"ant", -1, 1),
        Among(u"ent", -1, 1),
        Among(u"int", -1, 1),
        Among(u"ut", -1, 1),
        Among(u"ït", -1, 1),
        Among(u"au", -1, 1),
        Among(u"erau", 235, 1),
        Among(u"ieu", -1, 1),
        Among(u"ineu", -1, 1),
        Among(u"areu", -1, 1),
        Among(u"ireu", -1, 1),
        Among(u"àreu", -1, 1),
        Among(u"íreu", -1, 1),
        Among(u"asseu", -1, 1),
        Among(u"esseu", -1, 1),
        Among(u"eresseu", 244, 1),
        Among(u"àsseu", -1, 1),
        Among(u"ésseu", -1, 1),
        Among(u"igueu", -1, 1),
        Among(u"ïgueu", -1, 1),
        Among(u"àveu", -1, 1),
        Among(u"áveu", -1, 1),
        Among(u"itzeu", -1, 1),
        Among(u"ìeu", -1, 1),
        Among(u"irìeu", 253, 1),
        Among(u"íeu", -1, 1),
        Among(u"aríeu", 255, 1),
        Among(u"iríeu", 255, 1),
        Among(u"assiu", -1, 1),
        Among(u"issiu", -1, 1),
        Among(u"àssiu", -1, 1),
        Among(u"èssiu", -1, 1),
        Among(u"éssiu", -1, 1),
        Among(u"íssiu", -1, 1),
        Among(u"ïu", -1, 1),
        Among(u"ix", -1, 1),
        Among(u"eix", 265, 1),
        Among(u"ïx", -1, 1),
        Among(u"itz", -1, 1),
        Among(u"ià", -1, 1),
        Among(u"arà", -1, 1),
        Among(u"irà", -1, 1),
        Among(u"itzà", -1, 1),
        Among(u"ará", -1, 1),
        Among(u"erá", -1, 1),
        Among(u"irá", -1, 1),
        Among(u"irè", -1, 1),
        Among(u"aré", -1, 1),
        Among(u"eré", -1, 1),
        Among(u"iré", -1, 1),
        Among(u"í", -1, 1),
        Among(u"iï", -1, 1),
        Among(u"ió", -1, 1)
    ]

    a_4 = [
        Among(u"a", -1, 1),
        Among(u"e", -1, 1),
        Among(u"i", -1, 1),
        Among(u"ïn", -1, 1),
        Among(u"o", -1, 1),
        Among(u"ir", -1, 1),
        Among(u"s", -1, 1),
        Among(u"is", 6, 1),
        Among(u"os", 6, 1),
        Among(u"ïs", 6, 1),
        Among(u"it", -1, 1),
        Among(u"eu", -1, 1),
        Among(u"iu", -1, 1),
        Among(u"iqu", -1, 2),
        Among(u"itz", -1, 1),
        Among(u"à", -1, 1),
        Among(u"á", -1, 1),
        Among(u"é", -1, 1),
        Among(u"ì", -1, 1),
        Among(u"í", -1, 1),
        Among(u"ï", -1, 1),
        Among(u"ó", -1, 1)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass


class lab2(BaseException): pass

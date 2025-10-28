#-*- coding: utf-8 -*-
# Generated from tamil.sbl by Snowball 3.0.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class TamilStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from tamil.sbl by Snowball 3.0.1 - https://snowballstem.org/
    '''

    B_found_vetrumai_urupu = False
    B_found_a_match = False

    def __r_has_min_length(self):
        return len(self.current) > 4

    def __r_fix_va_start(self):
        self.bra = self.cursor
        among_var = self.find_among(TamilStemmer.a_0)
        if among_var == 0:
            return False
        self.ket = self.cursor
        if among_var == 1:
            if not self.slice_from(u"\u0B93"):
                return False
        elif among_var == 2:
            if not self.slice_from(u"\u0B92"):
                return False
        elif among_var == 3:
            if not self.slice_from(u"\u0B89"):
                return False
        else:
            if not self.slice_from(u"\u0B8A"):
                return False
        return True

    def __r_fix_endings(self):
        v_1 = self.cursor
        try:
            while True:
                v_2 = self.cursor
                try:
                    if not self.__r_fix_ending():
                        raise lab1()
                    continue
                except lab1: pass
                self.cursor = v_2
                break
        except lab0: pass
        self.cursor = v_1
        return True

    def __r_remove_question_prefixes(self):
        self.bra = self.cursor
        if not self.eq_s(u"\u0B8E"):
            return False
        if self.find_among(TamilStemmer.a_1) == 0:
            return False
        if not self.eq_s(u"\u0BCD"):
            return False
        self.ket = self.cursor
        if not self.slice_del():
            return False

        v_1 = self.cursor
        self.__r_fix_va_start()
        self.cursor = v_1
        return True

    def __r_fix_ending(self):
        if len(self.current) <= 3:
            return False
        self.limit_backward = self.cursor
        self.cursor = self.limit
        try:
            v_1 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                among_var = self.find_among_b(TamilStemmer.a_5)
                if among_var == 0:
                    raise lab1()
                self.bra = self.cursor
                if among_var == 1:
                    if not self.slice_del():
                        return False

                elif among_var == 2:
                    v_2 = self.limit - self.cursor
                    if self.find_among_b(TamilStemmer.a_2) == 0:
                        raise lab1()
                    self.cursor = self.limit - v_2
                    if not self.slice_del():
                        return False

                elif among_var == 3:
                    if not self.slice_from(u"\u0BB3\u0BCD"):
                        return False
                elif among_var == 4:
                    if not self.slice_from(u"\u0BB2\u0BCD"):
                        return False
                elif among_var == 5:
                    if not self.slice_from(u"\u0B9F\u0BC1"):
                        return False
                elif among_var == 6:
                    if not self.B_found_vetrumai_urupu:
                        raise lab1()
                    v_3 = self.limit - self.cursor
                    try:
                        if not self.eq_s_b(u"\u0BC8"):
                            raise lab2()
                        raise lab1()
                    except lab2: pass
                    self.cursor = self.limit - v_3
                    if not self.slice_from(u"\u0BAE\u0BCD"):
                        return False
                elif among_var == 7:
                    if not self.slice_from(u"\u0BCD"):
                        return False
                elif among_var == 8:
                    v_4 = self.limit - self.cursor
                    try:
                        if self.find_among_b(TamilStemmer.a_3) == 0:
                            raise lab3()
                        raise lab1()
                    except lab3: pass
                    self.cursor = self.limit - v_4
                    if not self.slice_del():
                        return False

                else:
                    among_var = self.find_among_b(TamilStemmer.a_4)
                    if among_var == 1:
                        if not self.slice_del():
                            return False

                    else:
                        if not self.slice_from(u"\u0BAE\u0BCD"):
                            return False
                raise lab0()
            except lab1: pass
            self.cursor = self.limit - v_1
            self.ket = self.cursor
            if not self.eq_s_b(u"\u0BCD"):
                return False
            try:
                v_5 = self.limit - self.cursor
                try:
                    if self.find_among_b(TamilStemmer.a_6) == 0:
                        raise lab5()
                    v_6 = self.limit - self.cursor
                    try:
                        if not self.eq_s_b(u"\u0BCD"):
                            self.cursor = self.limit - v_6
                            raise lab6()
                        if self.find_among_b(TamilStemmer.a_7) == 0:
                            self.cursor = self.limit - v_6
                            raise lab6()
                    except lab6: pass
                    self.bra = self.cursor
                    if not self.slice_del():
                        return False

                    raise lab4()
                except lab5: pass
                self.cursor = self.limit - v_5
                try:
                    if self.find_among_b(TamilStemmer.a_8) == 0:
                        raise lab7()
                    self.bra = self.cursor
                    if not self.eq_s_b(u"\u0BCD"):
                        raise lab7()
                    if not self.slice_del():
                        return False

                    raise lab4()
                except lab7: pass
                self.cursor = self.limit - v_5
                v_7 = self.limit - self.cursor
                if self.find_among_b(TamilStemmer.a_9) == 0:
                    return False
                self.cursor = self.limit - v_7
                self.bra = self.cursor
                if not self.slice_del():
                    return False

            except lab4: pass
        except lab0: pass
        self.cursor = self.limit_backward
        return True

    def __r_remove_pronoun_prefixes(self):
        self.bra = self.cursor
        if self.find_among(TamilStemmer.a_10) == 0:
            return False
        if self.find_among(TamilStemmer.a_11) == 0:
            return False
        if not self.eq_s(u"\u0BCD"):
            return False
        self.ket = self.cursor
        if not self.slice_del():
            return False

        v_1 = self.cursor
        self.__r_fix_va_start()
        self.cursor = v_1
        return True

    def __r_remove_plural_suffix(self):
        self.limit_backward = self.cursor
        self.cursor = self.limit
        self.ket = self.cursor
        among_var = self.find_among_b(TamilStemmer.a_13)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            try:
                v_1 = self.limit - self.cursor
                try:
                    if self.find_among_b(TamilStemmer.a_12) == 0:
                        raise lab1()
                    if not self.slice_from(u"\u0BC1\u0B99\u0BCD"):
                        return False
                    raise lab0()
                except lab1: pass
                self.cursor = self.limit - v_1
                if not self.slice_from(u"\u0BCD"):
                    return False
            except lab0: pass
        elif among_var == 2:
            if not self.slice_from(u"\u0BB2\u0BCD"):
                return False
        elif among_var == 3:
            if not self.slice_from(u"\u0BB3\u0BCD"):
                return False
        else:
            if not self.slice_del():
                return False

        self.cursor = self.limit_backward
        return True

    def __r_remove_question_suffixes(self):
        if not self.__r_has_min_length():
            return False
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_1 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if self.find_among_b(TamilStemmer.a_14) == 0:
                raise lab0()
            self.bra = self.cursor
            if not self.slice_from(u"\u0BCD"):
                return False
        except lab0: pass
        self.cursor = self.limit - v_1
        self.cursor = self.limit_backward
        self.__r_fix_endings()
        return True

    def __r_remove_command_suffixes(self):
        if not self.__r_has_min_length():
            return False
        self.limit_backward = self.cursor
        self.cursor = self.limit
        self.ket = self.cursor
        if self.find_among_b(TamilStemmer.a_15) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.cursor = self.limit_backward
        return True

    def __r_remove_um(self):
        if not self.__r_has_min_length():
            return False
        self.limit_backward = self.cursor
        self.cursor = self.limit
        self.ket = self.cursor
        if not self.eq_s_b(u"\u0BC1\u0BAE\u0BCD"):
            return False
        self.bra = self.cursor
        if not self.slice_from(u"\u0BCD"):
            return False
        self.cursor = self.limit_backward
        v_1 = self.cursor
        self.__r_fix_ending()
        self.cursor = v_1
        return True

    def __r_remove_common_word_endings(self):
        if not self.__r_has_min_length():
            return False
        self.limit_backward = self.cursor
        self.cursor = self.limit
        self.ket = self.cursor
        among_var = self.find_among_b(TamilStemmer.a_17)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if not self.slice_from(u"\u0BCD"):
                return False
        elif among_var == 2:
            v_1 = self.limit - self.cursor
            try:
                if self.find_among_b(TamilStemmer.a_16) == 0:
                    raise lab0()
                return False
            except lab0: pass
            self.cursor = self.limit - v_1
            if not self.slice_from(u"\u0BCD"):
                return False
        else:
            if not self.slice_del():
                return False

        self.cursor = self.limit_backward
        self.__r_fix_endings()
        return True

    def __r_remove_vetrumai_urupukal(self):
        self.B_found_vetrumai_urupu = False
        if not self.__r_has_min_length():
            return False
        self.limit_backward = self.cursor
        self.cursor = self.limit
        try:
            v_1 = self.limit - self.cursor
            try:
                v_2 = self.limit - self.cursor
                self.ket = self.cursor
                among_var = self.find_among_b(TamilStemmer.a_20)
                if among_var == 0:
                    raise lab1()
                self.bra = self.cursor
                if among_var == 1:
                    if not self.slice_del():
                        return False

                elif among_var == 2:
                    if not self.slice_from(u"\u0BCD"):
                        return False
                elif among_var == 3:
                    v_3 = self.limit - self.cursor
                    try:
                        if not self.eq_s_b(u"\u0BAE"):
                            raise lab2()
                        raise lab1()
                    except lab2: pass
                    self.cursor = self.limit - v_3
                    if not self.slice_from(u"\u0BCD"):
                        return False
                elif among_var == 4:
                    if len(self.current) < 7:
                        raise lab1()
                    if not self.slice_from(u"\u0BCD"):
                        return False
                elif among_var == 5:
                    v_4 = self.limit - self.cursor
                    try:
                        if self.find_among_b(TamilStemmer.a_18) == 0:
                            raise lab3()
                        raise lab1()
                    except lab3: pass
                    self.cursor = self.limit - v_4
                    if not self.slice_from(u"\u0BCD"):
                        return False
                elif among_var == 6:
                    v_5 = self.limit - self.cursor
                    try:
                        if self.find_among_b(TamilStemmer.a_19) == 0:
                            raise lab4()
                        raise lab1()
                    except lab4: pass
                    self.cursor = self.limit - v_5
                    if not self.slice_del():
                        return False

                else:
                    if not self.slice_from(u"\u0BBF"):
                        return False
                self.cursor = self.limit - v_2
                raise lab0()
            except lab1: pass
            self.cursor = self.limit - v_1
            v_6 = self.limit - self.cursor
            self.ket = self.cursor
            if not self.eq_s_b(u"\u0BC8"):
                return False
            try:
                v_7 = self.limit - self.cursor
                try:
                    v_8 = self.limit - self.cursor
                    try:
                        if self.find_among_b(TamilStemmer.a_21) == 0:
                            raise lab7()
                        raise lab6()
                    except lab7: pass
                    self.cursor = self.limit - v_8
                    raise lab5()
                except lab6: pass
                self.cursor = self.limit - v_7
                v_9 = self.limit - self.cursor
                if self.find_among_b(TamilStemmer.a_22) == 0:
                    return False
                if not self.eq_s_b(u"\u0BCD"):
                    return False
                self.cursor = self.limit - v_9
            except lab5: pass
            self.bra = self.cursor
            if not self.slice_from(u"\u0BCD"):
                return False
            self.cursor = self.limit - v_6
        except lab0: pass
        self.B_found_vetrumai_urupu = True
        v_10 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if not self.eq_s_b(u"\u0BBF\u0BA9\u0BCD"):
                raise lab8()
            self.bra = self.cursor
            if not self.slice_from(u"\u0BCD"):
                return False
        except lab8: pass
        self.cursor = self.limit - v_10
        self.cursor = self.limit_backward
        self.__r_fix_endings()
        return True

    def __r_remove_tense_suffixes(self):
        while True:
            v_1 = self.cursor
            try:
                if not self.__r_remove_tense_suffix():
                    raise lab0()
                continue
            except lab0: pass
            self.cursor = v_1
            break
        return True

    def __r_remove_tense_suffix(self):
        self.B_found_a_match = False
        if not self.__r_has_min_length():
            return False
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_1 = self.limit - self.cursor
        try:
            v_2 = self.limit - self.cursor
            self.ket = self.cursor
            among_var = self.find_among_b(TamilStemmer.a_25)
            if among_var == 0:
                raise lab0()
            self.bra = self.cursor
            if among_var == 1:
                if not self.slice_del():
                    return False

            elif among_var == 2:
                v_3 = self.limit - self.cursor
                try:
                    if self.find_among_b(TamilStemmer.a_23) == 0:
                        raise lab1()
                    raise lab0()
                except lab1: pass
                self.cursor = self.limit - v_3
                if not self.slice_del():
                    return False

            elif among_var == 3:
                v_4 = self.limit - self.cursor
                try:
                    if self.find_among_b(TamilStemmer.a_24) == 0:
                        raise lab2()
                    raise lab0()
                except lab2: pass
                self.cursor = self.limit - v_4
                if not self.slice_del():
                    return False

            elif among_var == 4:
                v_5 = self.limit - self.cursor
                try:
                    if not self.eq_s_b(u"\u0B9A"):
                        raise lab3()
                    raise lab0()
                except lab3: pass
                self.cursor = self.limit - v_5
                if not self.slice_from(u"\u0BCD"):
                    return False
            elif among_var == 5:
                if not self.slice_from(u"\u0BCD"):
                    return False
            else:
                v_6 = self.limit - self.cursor
                if not self.eq_s_b(u"\u0BCD"):
                    raise lab0()
                self.cursor = self.limit - v_6
                if not self.slice_del():
                    return False

            self.B_found_a_match = True
            self.cursor = self.limit - v_2
        except lab0: pass
        self.cursor = self.limit - v_1
        v_7 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if self.find_among_b(TamilStemmer.a_26) == 0:
                raise lab4()
            self.bra = self.cursor
            if not self.slice_del():
                return False

            self.B_found_a_match = True
        except lab4: pass
        self.cursor = self.limit - v_7
        self.cursor = self.limit_backward
        self.__r_fix_endings()
        if not self.B_found_a_match:
            return False
        return True

    def _stem(self):
        self.B_found_vetrumai_urupu = False
        v_1 = self.cursor
        self.__r_fix_ending()
        self.cursor = v_1
        if not self.__r_has_min_length():
            return False
        v_2 = self.cursor
        self.__r_remove_question_prefixes()
        self.cursor = v_2
        v_3 = self.cursor
        self.__r_remove_pronoun_prefixes()
        self.cursor = v_3
        v_4 = self.cursor
        self.__r_remove_question_suffixes()
        self.cursor = v_4
        v_5 = self.cursor
        self.__r_remove_um()
        self.cursor = v_5
        v_6 = self.cursor
        self.__r_remove_common_word_endings()
        self.cursor = v_6
        v_7 = self.cursor
        self.__r_remove_vetrumai_urupukal()
        self.cursor = v_7
        v_8 = self.cursor
        self.__r_remove_plural_suffix()
        self.cursor = v_8
        v_9 = self.cursor
        self.__r_remove_command_suffixes()
        self.cursor = v_9
        v_10 = self.cursor
        self.__r_remove_tense_suffixes()
        self.cursor = v_10
        return True

    a_0 = [
        Among(u"\u0BB5\u0BC1", -1, 3),
        Among(u"\u0BB5\u0BC2", -1, 4),
        Among(u"\u0BB5\u0BCA", -1, 2),
        Among(u"\u0BB5\u0BCB", -1, 1)
    ]

    a_1 = [
        Among(u"\u0B95", -1, -1),
        Among(u"\u0B99", -1, -1),
        Among(u"\u0B9A", -1, -1),
        Among(u"\u0B9E", -1, -1),
        Among(u"\u0BA4", -1, -1),
        Among(u"\u0BA8", -1, -1),
        Among(u"\u0BAA", -1, -1),
        Among(u"\u0BAE", -1, -1),
        Among(u"\u0BAF", -1, -1),
        Among(u"\u0BB5", -1, -1)
    ]

    a_2 = [
        Among(u"\u0BBF", -1, -1),
        Among(u"\u0BC0", -1, -1),
        Among(u"\u0BC8", -1, -1)
    ]

    a_3 = [
        Among(u"\u0BBE", -1, -1),
        Among(u"\u0BBF", -1, -1),
        Among(u"\u0BC0", -1, -1),
        Among(u"\u0BC1", -1, -1),
        Among(u"\u0BC2", -1, -1),
        Among(u"\u0BC6", -1, -1),
        Among(u"\u0BC7", -1, -1),
        Among(u"\u0BC8", -1, -1)
    ]

    a_4 = [
        Among(u"", -1, 2),
        Among(u"\u0BC8", 0, 1),
        Among(u"\u0BCD", 0, 1)
    ]

    a_5 = [
        Among(u"\u0BA8\u0BCD\u0BA4", -1, 1),
        Among(u"\u0BAF", -1, 1),
        Among(u"\u0BB5", -1, 1),
        Among(u"\u0BA9\u0BC1", -1, 8),
        Among(u"\u0BC1\u0B95\u0BCD", -1, 7),
        Among(u"\u0BC1\u0B95\u0BCD\u0B95\u0BCD", -1, 7),
        Among(u"\u0B9F\u0BCD\u0B95\u0BCD", -1, 3),
        Among(u"\u0BB1\u0BCD\u0B95\u0BCD", -1, 4),
        Among(u"\u0B99\u0BCD", -1, 9),
        Among(u"\u0B9F\u0BCD\u0B9F\u0BCD", -1, 5),
        Among(u"\u0BA4\u0BCD\u0BA4\u0BCD", -1, 6),
        Among(u"\u0BA8\u0BCD\u0BA4\u0BCD", -1, 1),
        Among(u"\u0BA8\u0BCD", -1, 1),
        Among(u"\u0B9F\u0BCD\u0BAA\u0BCD", -1, 3),
        Among(u"\u0BAF\u0BCD", -1, 2),
        Among(u"\u0BA9\u0BCD\u0BB1\u0BCD", -1, 4),
        Among(u"\u0BB5\u0BCD", -1, 1)
    ]

    a_6 = [
        Among(u"\u0B95", -1, -1),
        Among(u"\u0B9A", -1, -1),
        Among(u"\u0B9F", -1, -1),
        Among(u"\u0BA4", -1, -1),
        Among(u"\u0BAA", -1, -1),
        Among(u"\u0BB1", -1, -1)
    ]

    a_7 = [
        Among(u"\u0B95", -1, -1),
        Among(u"\u0B9A", -1, -1),
        Among(u"\u0B9F", -1, -1),
        Among(u"\u0BA4", -1, -1),
        Among(u"\u0BAA", -1, -1),
        Among(u"\u0BB1", -1, -1)
    ]

    a_8 = [
        Among(u"\u0B9E", -1, -1),
        Among(u"\u0BA3", -1, -1),
        Among(u"\u0BA8", -1, -1),
        Among(u"\u0BA9", -1, -1),
        Among(u"\u0BAE", -1, -1),
        Among(u"\u0BAF", -1, -1),
        Among(u"\u0BB0", -1, -1),
        Among(u"\u0BB2", -1, -1),
        Among(u"\u0BB3", -1, -1),
        Among(u"\u0BB4", -1, -1),
        Among(u"\u0BB5", -1, -1)
    ]

    a_9 = [
        Among(u"\u0BBE", -1, -1),
        Among(u"\u0BBF", -1, -1),
        Among(u"\u0BC0", -1, -1),
        Among(u"\u0BC1", -1, -1),
        Among(u"\u0BC2", -1, -1),
        Among(u"\u0BC6", -1, -1),
        Among(u"\u0BC7", -1, -1),
        Among(u"\u0BC8", -1, -1),
        Among(u"\u0BCD", -1, -1)
    ]

    a_10 = [
        Among(u"\u0B85", -1, -1),
        Among(u"\u0B87", -1, -1),
        Among(u"\u0B89", -1, -1)
    ]

    a_11 = [
        Among(u"\u0B95", -1, -1),
        Among(u"\u0B99", -1, -1),
        Among(u"\u0B9A", -1, -1),
        Among(u"\u0B9E", -1, -1),
        Among(u"\u0BA4", -1, -1),
        Among(u"\u0BA8", -1, -1),
        Among(u"\u0BAA", -1, -1),
        Among(u"\u0BAE", -1, -1),
        Among(u"\u0BAF", -1, -1),
        Among(u"\u0BB5", -1, -1)
    ]

    a_12 = [
        Among(u"\u0B95", -1, -1),
        Among(u"\u0B9A", -1, -1),
        Among(u"\u0B9F", -1, -1),
        Among(u"\u0BA4", -1, -1),
        Among(u"\u0BAA", -1, -1),
        Among(u"\u0BB1", -1, -1)
    ]

    a_13 = [
        Among(u"\u0B95\u0BB3\u0BCD", -1, 4),
        Among(u"\u0BC1\u0B99\u0BCD\u0B95\u0BB3\u0BCD", 0, 1),
        Among(u"\u0B9F\u0BCD\u0B95\u0BB3\u0BCD", 0, 3),
        Among(u"\u0BB1\u0BCD\u0B95\u0BB3\u0BCD", 0, 2)
    ]

    a_14 = [
        Among(u"\u0BBE", -1, -1),
        Among(u"\u0BC7", -1, -1),
        Among(u"\u0BCB", -1, -1)
    ]

    a_15 = [
        Among(u"\u0BAA\u0BBF", -1, -1),
        Among(u"\u0BB5\u0BBF", -1, -1)
    ]

    a_16 = [
        Among(u"\u0BBE", -1, -1),
        Among(u"\u0BBF", -1, -1),
        Among(u"\u0BC0", -1, -1),
        Among(u"\u0BC1", -1, -1),
        Among(u"\u0BC2", -1, -1),
        Among(u"\u0BC6", -1, -1),
        Among(u"\u0BC7", -1, -1),
        Among(u"\u0BC8", -1, -1)
    ]

    a_17 = [
        Among(u"\u0BAA\u0B9F\u0BCD\u0B9F", -1, 3),
        Among(u"\u0BAA\u0B9F\u0BCD\u0B9F\u0BA3", -1, 3),
        Among(u"\u0BA4\u0BBE\u0BA9", -1, 3),
        Among(u"\u0BAA\u0B9F\u0BBF\u0BA4\u0BBE\u0BA9", 2, 3),
        Among(u"\u0BC6\u0BA9", -1, 1),
        Among(u"\u0BBE\u0B95\u0BBF\u0BAF", -1, 1),
        Among(u"\u0B95\u0BC1\u0BB0\u0BBF\u0BAF", -1, 3),
        Among(u"\u0BC1\u0B9F\u0BC8\u0BAF", -1, 1),
        Among(u"\u0BB2\u0BCD\u0BB2", -1, 2),
        Among(u"\u0BC1\u0BB3\u0BCD\u0BB3", -1, 1),
        Among(u"\u0BBE\u0B95\u0BBF", -1, 1),
        Among(u"\u0BAA\u0B9F\u0BBF", -1, 3),
        Among(u"\u0BBF\u0BA9\u0BCD\u0BB1\u0BBF", -1, 1),
        Among(u"\u0BAA\u0BB1\u0BCD\u0BB1\u0BBF", -1, 3),
        Among(u"\u0BAA\u0B9F\u0BC1", -1, 3),
        Among(u"\u0BB5\u0BBF\u0B9F\u0BC1", -1, 3),
        Among(u"\u0BAA\u0B9F\u0BCD\u0B9F\u0BC1", -1, 3),
        Among(u"\u0BB5\u0BBF\u0B9F\u0BCD\u0B9F\u0BC1", -1, 3),
        Among(u"\u0BAA\u0B9F\u0BCD\u0B9F\u0BA4\u0BC1", -1, 3),
        Among(u"\u0BC6\u0BA9\u0BCD\u0BB1\u0BC1", -1, 1),
        Among(u"\u0BC1\u0B9F\u0BC8", -1, 1),
        Among(u"\u0BBF\u0BB2\u0BCD\u0BB2\u0BC8", -1, 1),
        Among(u"\u0BC1\u0B9F\u0BA9\u0BCD", -1, 1),
        Among(u"\u0BBF\u0B9F\u0BAE\u0BCD", -1, 1),
        Among(u"\u0BC6\u0BB2\u0BCD\u0BB2\u0BBE\u0BAE\u0BCD", -1, 3),
        Among(u"\u0BC6\u0BA9\u0BC1\u0BAE\u0BCD", -1, 1)
    ]

    a_18 = [
        Among(u"\u0BBE", -1, -1),
        Among(u"\u0BBF", -1, -1),
        Among(u"\u0BC0", -1, -1),
        Among(u"\u0BC1", -1, -1),
        Among(u"\u0BC2", -1, -1),
        Among(u"\u0BC6", -1, -1),
        Among(u"\u0BC7", -1, -1),
        Among(u"\u0BC8", -1, -1)
    ]

    a_19 = [
        Among(u"\u0BBE", -1, -1),
        Among(u"\u0BBF", -1, -1),
        Among(u"\u0BC0", -1, -1),
        Among(u"\u0BC1", -1, -1),
        Among(u"\u0BC2", -1, -1),
        Among(u"\u0BC6", -1, -1),
        Among(u"\u0BC7", -1, -1),
        Among(u"\u0BC8", -1, -1)
    ]

    a_20 = [
        Among(u"\u0BB5\u0BBF\u0B9F", -1, 2),
        Among(u"\u0BC0", -1, 7),
        Among(u"\u0BCA\u0B9F\u0BC1", -1, 2),
        Among(u"\u0BCB\u0B9F\u0BC1", -1, 2),
        Among(u"\u0BA4\u0BC1", -1, 6),
        Among(u"\u0BBF\u0BB0\u0BC1\u0BA8\u0BCD\u0BA4\u0BC1", 4, 2),
        Among(u"\u0BBF\u0BA9\u0BCD\u0BB1\u0BC1", -1, 2),
        Among(u"\u0BC1\u0B9F\u0BC8", -1, 2),
        Among(u"\u0BA9\u0BC8", -1, 1),
        Among(u"\u0B95\u0BA3\u0BCD", -1, 1),
        Among(u"\u0BBF\u0BA9\u0BCD", -1, 3),
        Among(u"\u0BAE\u0BC1\u0BA9\u0BCD", -1, 1),
        Among(u"\u0BBF\u0B9F\u0BAE\u0BCD", -1, 4),
        Among(u"\u0BBF\u0BB1\u0BCD", -1, 2),
        Among(u"\u0BAE\u0BC7\u0BB1\u0BCD", -1, 1),
        Among(u"\u0BB2\u0BCD", -1, 5),
        Among(u"\u0BBE\u0BAE\u0BB2\u0BCD", 15, 2),
        Among(u"\u0BBE\u0BB2\u0BCD", 15, 2),
        Among(u"\u0BBF\u0BB2\u0BCD", 15, 2),
        Among(u"\u0BAE\u0BC7\u0BB2\u0BCD", 15, 1),
        Among(u"\u0BC1\u0BB3\u0BCD", -1, 2),
        Among(u"\u0B95\u0BC0\u0BB4\u0BCD", -1, 1)
    ]

    a_21 = [
        Among(u"\u0B95", -1, -1),
        Among(u"\u0B9A", -1, -1),
        Among(u"\u0B9F", -1, -1),
        Among(u"\u0BA4", -1, -1),
        Among(u"\u0BAA", -1, -1),
        Among(u"\u0BB1", -1, -1)
    ]

    a_22 = [
        Among(u"\u0B95", -1, -1),
        Among(u"\u0B9A", -1, -1),
        Among(u"\u0B9F", -1, -1),
        Among(u"\u0BA4", -1, -1),
        Among(u"\u0BAA", -1, -1),
        Among(u"\u0BB1", -1, -1)
    ]

    a_23 = [
        Among(u"\u0B85", -1, -1),
        Among(u"\u0B86", -1, -1),
        Among(u"\u0B87", -1, -1),
        Among(u"\u0B88", -1, -1),
        Among(u"\u0B89", -1, -1),
        Among(u"\u0B8A", -1, -1),
        Among(u"\u0B8E", -1, -1),
        Among(u"\u0B8F", -1, -1),
        Among(u"\u0B90", -1, -1),
        Among(u"\u0B92", -1, -1),
        Among(u"\u0B93", -1, -1),
        Among(u"\u0B94", -1, -1)
    ]

    a_24 = [
        Among(u"\u0BBE", -1, -1),
        Among(u"\u0BBF", -1, -1),
        Among(u"\u0BC0", -1, -1),
        Among(u"\u0BC1", -1, -1),
        Among(u"\u0BC2", -1, -1),
        Among(u"\u0BC6", -1, -1),
        Among(u"\u0BC7", -1, -1),
        Among(u"\u0BC8", -1, -1)
    ]

    a_25 = [
        Among(u"\u0B95", -1, 1),
        Among(u"\u0BA4", -1, 1),
        Among(u"\u0BA9", -1, 1),
        Among(u"\u0BAA", -1, 1),
        Among(u"\u0BAF", -1, 1),
        Among(u"\u0BBE", -1, 5),
        Among(u"\u0B95\u0BC1", -1, 6),
        Among(u"\u0BAA\u0B9F\u0BC1", -1, 1),
        Among(u"\u0BA4\u0BC1", -1, 3),
        Among(u"\u0BBF\u0BB1\u0BCD\u0BB1\u0BC1", -1, 1),
        Among(u"\u0BA9\u0BC8", -1, 1),
        Among(u"\u0BB5\u0BC8", -1, 1),
        Among(u"\u0BA9\u0BA9\u0BCD", -1, 1),
        Among(u"\u0BAA\u0BA9\u0BCD", -1, 1),
        Among(u"\u0BB5\u0BA9\u0BCD", -1, 2),
        Among(u"\u0BBE\u0BA9\u0BCD", -1, 4),
        Among(u"\u0BA9\u0BBE\u0BA9\u0BCD", 15, 1),
        Among(u"\u0BAE\u0BBF\u0BA9\u0BCD", -1, 1),
        Among(u"\u0BA9\u0BC6\u0BA9\u0BCD", -1, 1),
        Among(u"\u0BC7\u0BA9\u0BCD", -1, 5),
        Among(u"\u0BA9\u0BAE\u0BCD", -1, 1),
        Among(u"\u0BAA\u0BAE\u0BCD", -1, 1),
        Among(u"\u0BBE\u0BAE\u0BCD", -1, 5),
        Among(u"\u0B95\u0BC1\u0BAE\u0BCD", -1, 1),
        Among(u"\u0B9F\u0BC1\u0BAE\u0BCD", -1, 5),
        Among(u"\u0BA4\u0BC1\u0BAE\u0BCD", -1, 1),
        Among(u"\u0BB1\u0BC1\u0BAE\u0BCD", -1, 1),
        Among(u"\u0BC6\u0BAE\u0BCD", -1, 5),
        Among(u"\u0BC7\u0BAE\u0BCD", -1, 5),
        Among(u"\u0BCB\u0BAE\u0BCD", -1, 5),
        Among(u"\u0BBE\u0BAF\u0BCD", -1, 5),
        Among(u"\u0BA9\u0BB0\u0BCD", -1, 1),
        Among(u"\u0BAA\u0BB0\u0BCD", -1, 1),
        Among(u"\u0BC0\u0BAF\u0BB0\u0BCD", -1, 5),
        Among(u"\u0BB5\u0BB0\u0BCD", -1, 1),
        Among(u"\u0BBE\u0BB0\u0BCD", -1, 5),
        Among(u"\u0BA9\u0BBE\u0BB0\u0BCD", 35, 1),
        Among(u"\u0BAE\u0BBE\u0BB0\u0BCD", 35, 1),
        Among(u"\u0B95\u0BCA\u0BA3\u0BCD\u0B9F\u0BBF\u0BB0\u0BCD", -1, 1),
        Among(u"\u0BA9\u0BBF\u0BB0\u0BCD", -1, 5),
        Among(u"\u0BC0\u0BB0\u0BCD", -1, 5),
        Among(u"\u0BA9\u0BB3\u0BCD", -1, 1),
        Among(u"\u0BAA\u0BB3\u0BCD", -1, 1),
        Among(u"\u0BB5\u0BB3\u0BCD", -1, 1),
        Among(u"\u0BBE\u0BB3\u0BCD", -1, 5),
        Among(u"\u0BA9\u0BBE\u0BB3\u0BCD", 44, 1)
    ]

    a_26 = [
        Among(u"\u0B95\u0BBF\u0BB1", -1, -1),
        Among(u"\u0B95\u0BBF\u0BA9\u0BCD\u0BB1", -1, -1),
        Among(u"\u0BBE\u0BA8\u0BBF\u0BA9\u0BCD\u0BB1", -1, -1),
        Among(u"\u0B95\u0BBF\u0BB1\u0BCD", -1, -1),
        Among(u"\u0B95\u0BBF\u0BA9\u0BCD\u0BB1\u0BCD", -1, -1),
        Among(u"\u0BBE\u0BA8\u0BBF\u0BA9\u0BCD\u0BB1\u0BCD", -1, -1)
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

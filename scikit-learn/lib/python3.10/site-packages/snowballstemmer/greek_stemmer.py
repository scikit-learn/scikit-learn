#-*- coding: utf-8 -*-
# Generated from greek.sbl by Snowball 3.0.1 - https://snowballstem.org/

from .basestemmer import BaseStemmer
from .among import Among


class GreekStemmer(BaseStemmer):
    '''
    This class implements the stemming algorithm defined by a snowball script.
    Generated from greek.sbl by Snowball 3.0.1 - https://snowballstem.org/
    '''

    g_v = {u"α", u"ε", u"η", u"ι", u"ο", u"υ", u"ω"}

    g_v2 = {u"α", u"ε", u"η", u"ι", u"ο", u"ω"}

    B_test1 = False

    def __r_has_min_length(self):
        return len(self.current) >= 3

    def __r_tolower(self):
        while True:
            v_1 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                among_var = self.find_among_b(GreekStemmer.a_0)
                self.bra = self.cursor
                if among_var == 1:
                    if not self.slice_from(u"α"):
                        return False
                elif among_var == 2:
                    if not self.slice_from(u"β"):
                        return False
                elif among_var == 3:
                    if not self.slice_from(u"γ"):
                        return False
                elif among_var == 4:
                    if not self.slice_from(u"δ"):
                        return False
                elif among_var == 5:
                    if not self.slice_from(u"ε"):
                        return False
                elif among_var == 6:
                    if not self.slice_from(u"ζ"):
                        return False
                elif among_var == 7:
                    if not self.slice_from(u"η"):
                        return False
                elif among_var == 8:
                    if not self.slice_from(u"θ"):
                        return False
                elif among_var == 9:
                    if not self.slice_from(u"ι"):
                        return False
                elif among_var == 10:
                    if not self.slice_from(u"κ"):
                        return False
                elif among_var == 11:
                    if not self.slice_from(u"λ"):
                        return False
                elif among_var == 12:
                    if not self.slice_from(u"μ"):
                        return False
                elif among_var == 13:
                    if not self.slice_from(u"ν"):
                        return False
                elif among_var == 14:
                    if not self.slice_from(u"ξ"):
                        return False
                elif among_var == 15:
                    if not self.slice_from(u"ο"):
                        return False
                elif among_var == 16:
                    if not self.slice_from(u"π"):
                        return False
                elif among_var == 17:
                    if not self.slice_from(u"ρ"):
                        return False
                elif among_var == 18:
                    if not self.slice_from(u"σ"):
                        return False
                elif among_var == 19:
                    if not self.slice_from(u"τ"):
                        return False
                elif among_var == 20:
                    if not self.slice_from(u"υ"):
                        return False
                elif among_var == 21:
                    if not self.slice_from(u"φ"):
                        return False
                elif among_var == 22:
                    if not self.slice_from(u"χ"):
                        return False
                elif among_var == 23:
                    if not self.slice_from(u"ψ"):
                        return False
                elif among_var == 24:
                    if not self.slice_from(u"ω"):
                        return False
                else:
                    if self.cursor <= self.limit_backward:
                        raise lab0()
                    self.cursor -= 1
                continue
            except lab0: pass
            self.cursor = self.limit - v_1
            break
        return True

    def __r_step_1(self):
        self.ket = self.cursor
        among_var = self.find_among_b(GreekStemmer.a_1)
        if among_var == 0:
            return False
        self.bra = self.cursor
        if among_var == 1:
            if not self.slice_from(u"φα"):
                return False
        elif among_var == 2:
            if not self.slice_from(u"σκα"):
                return False
        elif among_var == 3:
            if not self.slice_from(u"ολο"):
                return False
        elif among_var == 4:
            if not self.slice_from(u"σο"):
                return False
        elif among_var == 5:
            if not self.slice_from(u"τατο"):
                return False
        elif among_var == 6:
            if not self.slice_from(u"κρε"):
                return False
        elif among_var == 7:
            if not self.slice_from(u"περ"):
                return False
        elif among_var == 8:
            if not self.slice_from(u"τερ"):
                return False
        elif among_var == 9:
            if not self.slice_from(u"φω"):
                return False
        elif among_var == 10:
            if not self.slice_from(u"καθεστ"):
                return False
        else:
            if not self.slice_from(u"γεγον"):
                return False
        self.B_test1 = False
        return True

    def __r_step_s1(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_3) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        self.ket = self.cursor
        self.bra = self.cursor
        among_var = self.find_among_b(GreekStemmer.a_2)
        if among_var == 0:
            return False
        if self.cursor > self.limit_backward:
            return False
        if among_var == 1:
            if not self.slice_from(u"ι"):
                return False
        else:
            if not self.slice_from(u"ιζ"):
                return False
        return True

    def __r_step_s2(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_5) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        self.ket = self.cursor
        self.bra = self.cursor
        if self.find_among_b(GreekStemmer.a_4) == 0:
            return False
        if self.cursor > self.limit_backward:
            return False
        if not self.slice_from(u"ων"):
            return False
        return True

    def __r_step_s3(self):
        try:
            v_1 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                if not self.eq_s_b(u"ισα"):
                    raise lab1()
                self.bra = self.cursor
                if self.cursor > self.limit_backward:
                    raise lab1()
                if not self.slice_from(u"ισ"):
                    return False
                raise lab0()
            except lab1: pass
            self.cursor = self.limit - v_1
            self.ket = self.cursor
        except lab0: pass
        if self.find_among_b(GreekStemmer.a_7) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        self.ket = self.cursor
        self.bra = self.cursor
        among_var = self.find_among_b(GreekStemmer.a_6)
        if among_var == 0:
            return False
        if self.cursor > self.limit_backward:
            return False
        if among_var == 1:
            if not self.slice_from(u"ι"):
                return False
        else:
            if not self.slice_from(u"ισ"):
                return False
        return True

    def __r_step_s4(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_9) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        self.ket = self.cursor
        self.bra = self.cursor
        if self.find_among_b(GreekStemmer.a_8) == 0:
            return False
        if self.cursor > self.limit_backward:
            return False
        if not self.slice_from(u"ι"):
            return False
        return True

    def __r_step_s5(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_11) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        self.ket = self.cursor
        self.bra = self.cursor
        among_var = self.find_among_b(GreekStemmer.a_10)
        if among_var == 0:
            return False
        if self.cursor > self.limit_backward:
            return False
        if among_var == 1:
            if not self.slice_from(u"ι"):
                return False
        else:
            if not self.slice_from(u"ιστ"):
                return False
        return True

    def __r_step_s6(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_14) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        try:
            v_1 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                self.bra = self.cursor
                among_var = self.find_among_b(GreekStemmer.a_12)
                if among_var == 0:
                    raise lab1()
                if self.cursor > self.limit_backward:
                    raise lab1()
                if among_var == 1:
                    if not self.slice_from(u"ισμ"):
                        return False
                else:
                    if not self.slice_from(u"ι"):
                        return False
                raise lab0()
            except lab1: pass
            self.cursor = self.limit - v_1
            self.ket = self.cursor
            among_var = self.find_among_b(GreekStemmer.a_13)
            if among_var == 0:
                return False
            self.bra = self.cursor
            if among_var == 1:
                if not self.slice_from(u"αγνωστ"):
                    return False
            elif among_var == 2:
                if not self.slice_from(u"ατομ"):
                    return False
            elif among_var == 3:
                if not self.slice_from(u"γνωστ"):
                    return False
            elif among_var == 4:
                if not self.slice_from(u"εθν"):
                    return False
            elif among_var == 5:
                if not self.slice_from(u"εκλεκτ"):
                    return False
            elif among_var == 6:
                if not self.slice_from(u"σκεπτ"):
                    return False
            elif among_var == 7:
                if not self.slice_from(u"τοπ"):
                    return False
            elif among_var == 8:
                if not self.slice_from(u"αλεξανδρ"):
                    return False
            elif among_var == 9:
                if not self.slice_from(u"βυζαντ"):
                    return False
            else:
                if not self.slice_from(u"θεατρ"):
                    return False
        except lab0: pass
        return True

    def __r_step_s7(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_16) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        self.ket = self.cursor
        self.bra = self.cursor
        if self.find_among_b(GreekStemmer.a_15) == 0:
            return False
        if self.cursor > self.limit_backward:
            return False
        if not self.slice_from(u"αρακ"):
            return False
        return True

    def __r_step_s8(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_18) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        try:
            v_1 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                self.bra = self.cursor
                among_var = self.find_among_b(GreekStemmer.a_17)
                if among_var == 0:
                    raise lab1()
                if self.cursor > self.limit_backward:
                    raise lab1()
                if among_var == 1:
                    if not self.slice_from(u"ακ"):
                        return False
                else:
                    if not self.slice_from(u"ιτσ"):
                        return False
                raise lab0()
            except lab1: pass
            self.cursor = self.limit - v_1
            self.ket = self.cursor
            self.bra = self.cursor
            if not self.eq_s_b(u"κορ"):
                return False
            if not self.slice_from(u"ιτσ"):
                return False
        except lab0: pass
        return True

    def __r_step_s9(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_21) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        try:
            v_1 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                self.bra = self.cursor
                if self.find_among_b(GreekStemmer.a_19) == 0:
                    raise lab1()
                if self.cursor > self.limit_backward:
                    raise lab1()
                if not self.slice_from(u"ιδ"):
                    return False
                raise lab0()
            except lab1: pass
            self.cursor = self.limit - v_1
            self.ket = self.cursor
            self.bra = self.cursor
            if self.find_among_b(GreekStemmer.a_20) == 0:
                return False
            if not self.slice_from(u"ιδ"):
                return False
        except lab0: pass
        return True

    def __r_step_s10(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_23) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        self.ket = self.cursor
        self.bra = self.cursor
        if self.find_among_b(GreekStemmer.a_22) == 0:
            return False
        if self.cursor > self.limit_backward:
            return False
        if not self.slice_from(u"ισκ"):
            return False
        return True

    def __r_step_2a(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_24) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        v_1 = self.limit - self.cursor
        try:
            if self.find_among_b(GreekStemmer.a_25) == 0:
                raise lab0()
            return False
        except lab0: pass
        self.cursor = self.limit - v_1
        c = self.cursor
        self.insert(self.cursor, self.cursor, u"αδ")
        self.cursor = c
        return True

    def __r_step_2b(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_26) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.ket = self.cursor
        self.bra = self.cursor
        if self.find_among_b(GreekStemmer.a_27) == 0:
            return False
        if not self.slice_from(u"εδ"):
            return False
        return True

    def __r_step_2c(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_28) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.ket = self.cursor
        self.bra = self.cursor
        if self.find_among_b(GreekStemmer.a_29) == 0:
            return False
        if not self.slice_from(u"ουδ"):
            return False
        return True

    def __r_step_2d(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_30) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        self.ket = self.cursor
        self.bra = self.cursor
        if self.find_among_b(GreekStemmer.a_31) == 0:
            return False
        if self.cursor > self.limit_backward:
            return False
        if not self.slice_from(u"ε"):
            return False
        return True

    def __r_step_3(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_32) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        self.ket = self.cursor
        self.bra = self.cursor
        if not self.in_grouping_b(GreekStemmer.g_v):
            return False
        if not self.slice_from(u"ι"):
            return False
        return True

    def __r_step_4(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_33) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        try:
            v_1 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                self.bra = self.cursor
                if not self.in_grouping_b(GreekStemmer.g_v):
                    raise lab1()
                if not self.slice_from(u"ικ"):
                    return False
                raise lab0()
            except lab1: pass
            self.cursor = self.limit - v_1
            self.ket = self.cursor
        except lab0: pass
        self.bra = self.cursor
        if self.find_among_b(GreekStemmer.a_34) == 0:
            return False
        if self.cursor > self.limit_backward:
            return False
        if not self.slice_from(u"ικ"):
            return False
        return True

    def __r_step_5a(self):
        v_1 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if not self.eq_s_b(u"αγαμε"):
                raise lab0()
            self.bra = self.cursor
            if self.cursor > self.limit_backward:
                raise lab0()
            if not self.slice_from(u"αγαμ"):
                return False
        except lab0: pass
        self.cursor = self.limit - v_1
        v_2 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if self.find_among_b(GreekStemmer.a_35) == 0:
                raise lab1()
            self.bra = self.cursor
            if not self.slice_del():
                return False

            self.B_test1 = False
        except lab1: pass
        self.cursor = self.limit - v_2
        self.ket = self.cursor
        if not self.eq_s_b(u"αμε"):
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        self.ket = self.cursor
        self.bra = self.cursor
        if self.find_among_b(GreekStemmer.a_36) == 0:
            return False
        if self.cursor > self.limit_backward:
            return False
        if not self.slice_from(u"αμ"):
            return False
        return True

    def __r_step_5b(self):
        v_1 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if self.find_among_b(GreekStemmer.a_38) == 0:
                raise lab0()
            self.bra = self.cursor
            if not self.slice_del():
                return False

            self.B_test1 = False
            self.ket = self.cursor
            self.bra = self.cursor
            if self.find_among_b(GreekStemmer.a_37) == 0:
                raise lab0()
            if self.cursor > self.limit_backward:
                raise lab0()
            if not self.slice_from(u"αγαν"):
                return False
        except lab0: pass
        self.cursor = self.limit - v_1
        self.ket = self.cursor
        if not self.eq_s_b(u"ανε"):
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        try:
            v_2 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                self.bra = self.cursor
                if not self.in_grouping_b(GreekStemmer.g_v2):
                    raise lab2()
                if not self.slice_from(u"αν"):
                    return False
                raise lab1()
            except lab2: pass
            self.cursor = self.limit - v_2
            self.ket = self.cursor
        except lab1: pass
        self.bra = self.cursor
        if self.find_among_b(GreekStemmer.a_39) == 0:
            return False
        if self.cursor > self.limit_backward:
            return False
        if not self.slice_from(u"αν"):
            return False
        return True

    def __r_step_5c(self):
        v_1 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if self.find_among_b(GreekStemmer.a_40) == 0:
                raise lab0()
            self.bra = self.cursor
            if not self.slice_del():
                return False

            self.B_test1 = False
        except lab0: pass
        self.cursor = self.limit - v_1
        self.ket = self.cursor
        if not self.eq_s_b(u"ετε"):
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        try:
            v_2 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                self.bra = self.cursor
                if not self.in_grouping_b(GreekStemmer.g_v2):
                    raise lab2()
                if not self.slice_from(u"ετ"):
                    return False
                raise lab1()
            except lab2: pass
            self.cursor = self.limit - v_2
            try:
                self.ket = self.cursor
                self.bra = self.cursor
                if self.find_among_b(GreekStemmer.a_41) == 0:
                    raise lab3()
                if not self.slice_from(u"ετ"):
                    return False
                raise lab1()
            except lab3: pass
            self.cursor = self.limit - v_2
            self.ket = self.cursor
        except lab1: pass
        self.bra = self.cursor
        if self.find_among_b(GreekStemmer.a_42) == 0:
            return False
        if self.cursor > self.limit_backward:
            return False
        if not self.slice_from(u"ετ"):
            return False
        return True

    def __r_step_5d(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_43) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        try:
            v_1 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                self.bra = self.cursor
                if not self.eq_s_b(u"αρχ"):
                    raise lab1()
                if self.cursor > self.limit_backward:
                    raise lab1()
                if not self.slice_from(u"οντ"):
                    return False
                raise lab0()
            except lab1: pass
            self.cursor = self.limit - v_1
            self.ket = self.cursor
            self.bra = self.cursor
            if not self.eq_s_b(u"κρε"):
                return False
            if not self.slice_from(u"ωντ"):
                return False
        except lab0: pass
        return True

    def __r_step_5e(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_44) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        self.ket = self.cursor
        self.bra = self.cursor
        if not self.eq_s_b(u"ον"):
            return False
        if self.cursor > self.limit_backward:
            return False
        if not self.slice_from(u"ομαστ"):
            return False
        return True

    def __r_step_5f(self):
        v_1 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if not self.eq_s_b(u"ιεστε"):
                raise lab0()
            self.bra = self.cursor
            if not self.slice_del():
                return False

            self.B_test1 = False
            self.ket = self.cursor
            self.bra = self.cursor
            if self.find_among_b(GreekStemmer.a_45) == 0:
                raise lab0()
            if self.cursor > self.limit_backward:
                raise lab0()
            if not self.slice_from(u"ιεστ"):
                return False
        except lab0: pass
        self.cursor = self.limit - v_1
        self.ket = self.cursor
        if not self.eq_s_b(u"εστε"):
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        self.ket = self.cursor
        self.bra = self.cursor
        if self.find_among_b(GreekStemmer.a_46) == 0:
            return False
        if self.cursor > self.limit_backward:
            return False
        if not self.slice_from(u"ιεστ"):
            return False
        return True

    def __r_step_5g(self):
        v_1 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if self.find_among_b(GreekStemmer.a_47) == 0:
                raise lab0()
            self.bra = self.cursor
            if not self.slice_del():
                return False

            self.B_test1 = False
        except lab0: pass
        self.cursor = self.limit - v_1
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_50) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        try:
            v_2 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                self.bra = self.cursor
                if self.find_among_b(GreekStemmer.a_48) == 0:
                    raise lab2()
                if not self.slice_from(u"ηκ"):
                    return False
                raise lab1()
            except lab2: pass
            self.cursor = self.limit - v_2
            self.ket = self.cursor
            self.bra = self.cursor
            if self.find_among_b(GreekStemmer.a_49) == 0:
                return False
            if self.cursor > self.limit_backward:
                return False
            if not self.slice_from(u"ηκ"):
                return False
        except lab1: pass
        return True

    def __r_step_5h(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_53) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        try:
            v_1 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                self.bra = self.cursor
                if self.find_among_b(GreekStemmer.a_51) == 0:
                    raise lab1()
                if not self.slice_from(u"ουσ"):
                    return False
                raise lab0()
            except lab1: pass
            self.cursor = self.limit - v_1
            self.ket = self.cursor
            self.bra = self.cursor
            if self.find_among_b(GreekStemmer.a_52) == 0:
                return False
            if self.cursor > self.limit_backward:
                return False
            if not self.slice_from(u"ουσ"):
                return False
        except lab0: pass
        return True

    def __r_step_5i(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_56) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        try:
            v_1 = self.limit - self.cursor
            try:
                self.ket = self.cursor
                self.bra = self.cursor
                if not self.eq_s_b(u"κολλ"):
                    raise lab1()
                if not self.slice_from(u"αγ"):
                    return False
                raise lab0()
            except lab1: pass
            self.cursor = self.limit - v_1
            try:
                v_2 = self.limit - self.cursor
                try:
                    self.ket = self.cursor
                    self.bra = self.cursor
                    among_var = self.find_among_b(GreekStemmer.a_54)
                    if among_var == 0:
                        raise lab3()
                    if among_var == 1:
                        if not self.slice_from(u"αγ"):
                            return False
                    raise lab2()
                except lab3: pass
                self.cursor = self.limit - v_2
                self.ket = self.cursor
                self.bra = self.cursor
                if self.find_among_b(GreekStemmer.a_55) == 0:
                    return False
                if self.cursor > self.limit_backward:
                    return False
                if not self.slice_from(u"αγ"):
                    return False
            except lab2: pass
        except lab0: pass
        return True

    def __r_step_5j(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_57) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        self.ket = self.cursor
        self.bra = self.cursor
        if self.find_among_b(GreekStemmer.a_58) == 0:
            return False
        if self.cursor > self.limit_backward:
            return False
        if not self.slice_from(u"ησ"):
            return False
        return True

    def __r_step_5k(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_59) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        self.ket = self.cursor
        self.bra = self.cursor
        if self.find_among_b(GreekStemmer.a_60) == 0:
            return False
        if self.cursor > self.limit_backward:
            return False
        if not self.slice_from(u"ηστ"):
            return False
        return True

    def __r_step_5l(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_61) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        self.ket = self.cursor
        self.bra = self.cursor
        if self.find_among_b(GreekStemmer.a_62) == 0:
            return False
        if self.cursor > self.limit_backward:
            return False
        if not self.slice_from(u"ουν"):
            return False
        return True

    def __r_step_5m(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_63) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        self.B_test1 = False
        self.ket = self.cursor
        self.bra = self.cursor
        if self.find_among_b(GreekStemmer.a_64) == 0:
            return False
        if self.cursor > self.limit_backward:
            return False
        if not self.slice_from(u"ουμ"):
            return False
        return True

    def __r_step_6(self):
        v_1 = self.limit - self.cursor
        try:
            self.ket = self.cursor
            if self.find_among_b(GreekStemmer.a_65) == 0:
                raise lab0()
            self.bra = self.cursor
            if not self.slice_from(u"μα"):
                return False
        except lab0: pass
        self.cursor = self.limit - v_1
        if not self.B_test1:
            return False
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_66) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        return True

    def __r_step_7(self):
        self.ket = self.cursor
        if self.find_among_b(GreekStemmer.a_67) == 0:
            return False
        self.bra = self.cursor
        if not self.slice_del():
            return False

        return True

    def _stem(self):
        self.limit_backward = self.cursor
        self.cursor = self.limit
        v_1 = self.limit - self.cursor
        self.__r_tolower()
        self.cursor = self.limit - v_1
        if not self.__r_has_min_length():
            return False
        self.B_test1 = True
        v_2 = self.limit - self.cursor
        self.__r_step_1()
        self.cursor = self.limit - v_2
        v_3 = self.limit - self.cursor
        self.__r_step_s1()
        self.cursor = self.limit - v_3
        v_4 = self.limit - self.cursor
        self.__r_step_s2()
        self.cursor = self.limit - v_4
        v_5 = self.limit - self.cursor
        self.__r_step_s3()
        self.cursor = self.limit - v_5
        v_6 = self.limit - self.cursor
        self.__r_step_s4()
        self.cursor = self.limit - v_6
        v_7 = self.limit - self.cursor
        self.__r_step_s5()
        self.cursor = self.limit - v_7
        v_8 = self.limit - self.cursor
        self.__r_step_s6()
        self.cursor = self.limit - v_8
        v_9 = self.limit - self.cursor
        self.__r_step_s7()
        self.cursor = self.limit - v_9
        v_10 = self.limit - self.cursor
        self.__r_step_s8()
        self.cursor = self.limit - v_10
        v_11 = self.limit - self.cursor
        self.__r_step_s9()
        self.cursor = self.limit - v_11
        v_12 = self.limit - self.cursor
        self.__r_step_s10()
        self.cursor = self.limit - v_12
        v_13 = self.limit - self.cursor
        self.__r_step_2a()
        self.cursor = self.limit - v_13
        v_14 = self.limit - self.cursor
        self.__r_step_2b()
        self.cursor = self.limit - v_14
        v_15 = self.limit - self.cursor
        self.__r_step_2c()
        self.cursor = self.limit - v_15
        v_16 = self.limit - self.cursor
        self.__r_step_2d()
        self.cursor = self.limit - v_16
        v_17 = self.limit - self.cursor
        self.__r_step_3()
        self.cursor = self.limit - v_17
        v_18 = self.limit - self.cursor
        self.__r_step_4()
        self.cursor = self.limit - v_18
        v_19 = self.limit - self.cursor
        self.__r_step_5a()
        self.cursor = self.limit - v_19
        v_20 = self.limit - self.cursor
        self.__r_step_5b()
        self.cursor = self.limit - v_20
        v_21 = self.limit - self.cursor
        self.__r_step_5c()
        self.cursor = self.limit - v_21
        v_22 = self.limit - self.cursor
        self.__r_step_5d()
        self.cursor = self.limit - v_22
        v_23 = self.limit - self.cursor
        self.__r_step_5e()
        self.cursor = self.limit - v_23
        v_24 = self.limit - self.cursor
        self.__r_step_5f()
        self.cursor = self.limit - v_24
        v_25 = self.limit - self.cursor
        self.__r_step_5g()
        self.cursor = self.limit - v_25
        v_26 = self.limit - self.cursor
        self.__r_step_5h()
        self.cursor = self.limit - v_26
        v_27 = self.limit - self.cursor
        self.__r_step_5j()
        self.cursor = self.limit - v_27
        v_28 = self.limit - self.cursor
        self.__r_step_5i()
        self.cursor = self.limit - v_28
        v_29 = self.limit - self.cursor
        self.__r_step_5k()
        self.cursor = self.limit - v_29
        v_30 = self.limit - self.cursor
        self.__r_step_5l()
        self.cursor = self.limit - v_30
        v_31 = self.limit - self.cursor
        self.__r_step_5m()
        self.cursor = self.limit - v_31
        v_32 = self.limit - self.cursor
        self.__r_step_6()
        self.cursor = self.limit - v_32
        v_33 = self.limit - self.cursor
        self.__r_step_7()
        self.cursor = self.limit - v_33
        self.cursor = self.limit_backward
        return True

    a_0 = [
        Among(u"", -1, 25),
        Among(u"Ά", 0, 1),
        Among(u"Έ", 0, 5),
        Among(u"Ή", 0, 7),
        Among(u"Ί", 0, 9),
        Among(u"Ό", 0, 15),
        Among(u"Ύ", 0, 20),
        Among(u"Ώ", 0, 24),
        Among(u"ΐ", 0, 7),
        Among(u"Α", 0, 1),
        Among(u"Β", 0, 2),
        Among(u"Γ", 0, 3),
        Among(u"Δ", 0, 4),
        Among(u"Ε", 0, 5),
        Among(u"Ζ", 0, 6),
        Among(u"Η", 0, 7),
        Among(u"Θ", 0, 8),
        Among(u"Ι", 0, 9),
        Among(u"Κ", 0, 10),
        Among(u"Λ", 0, 11),
        Among(u"Μ", 0, 12),
        Among(u"Ν", 0, 13),
        Among(u"Ξ", 0, 14),
        Among(u"Ο", 0, 15),
        Among(u"Π", 0, 16),
        Among(u"Ρ", 0, 17),
        Among(u"Σ", 0, 18),
        Among(u"Τ", 0, 19),
        Among(u"Υ", 0, 20),
        Among(u"Φ", 0, 21),
        Among(u"Χ", 0, 22),
        Among(u"Ψ", 0, 23),
        Among(u"Ω", 0, 24),
        Among(u"Ϊ", 0, 9),
        Among(u"Ϋ", 0, 20),
        Among(u"ά", 0, 1),
        Among(u"έ", 0, 5),
        Among(u"ή", 0, 7),
        Among(u"ί", 0, 9),
        Among(u"ΰ", 0, 20),
        Among(u"ς", 0, 18),
        Among(u"ϊ", 0, 7),
        Among(u"ϋ", 0, 20),
        Among(u"ό", 0, 15),
        Among(u"ύ", 0, 20),
        Among(u"ώ", 0, 24)
    ]

    a_1 = [
        Among(u"σκαγια", -1, 2),
        Among(u"φαγια", -1, 1),
        Among(u"ολογια", -1, 3),
        Among(u"σογια", -1, 4),
        Among(u"τατογια", -1, 5),
        Among(u"κρεατα", -1, 6),
        Among(u"περατα", -1, 7),
        Among(u"τερατα", -1, 8),
        Among(u"γεγονοτα", -1, 11),
        Among(u"καθεστωτα", -1, 10),
        Among(u"φωτα", -1, 9),
        Among(u"περατη", -1, 7),
        Among(u"σκαγιων", -1, 2),
        Among(u"φαγιων", -1, 1),
        Among(u"ολογιων", -1, 3),
        Among(u"σογιων", -1, 4),
        Among(u"τατογιων", -1, 5),
        Among(u"κρεατων", -1, 6),
        Among(u"περατων", -1, 7),
        Among(u"τερατων", -1, 8),
        Among(u"γεγονοτων", -1, 11),
        Among(u"καθεστωτων", -1, 10),
        Among(u"φωτων", -1, 9),
        Among(u"κρεασ", -1, 6),
        Among(u"περασ", -1, 7),
        Among(u"τερασ", -1, 8),
        Among(u"γεγονοσ", -1, 11),
        Among(u"κρεατοσ", -1, 6),
        Among(u"περατοσ", -1, 7),
        Among(u"τερατοσ", -1, 8),
        Among(u"γεγονοτοσ", -1, 11),
        Among(u"καθεστωτοσ", -1, 10),
        Among(u"φωτοσ", -1, 9),
        Among(u"καθεστωσ", -1, 10),
        Among(u"φωσ", -1, 9),
        Among(u"σκαγιου", -1, 2),
        Among(u"φαγιου", -1, 1),
        Among(u"ολογιου", -1, 3),
        Among(u"σογιου", -1, 4),
        Among(u"τατογιου", -1, 5)
    ]

    a_2 = [
        Among(u"πα", -1, 1),
        Among(u"ξαναπα", 0, 1),
        Among(u"επα", 0, 1),
        Among(u"περιπα", 0, 1),
        Among(u"αναμπα", 0, 1),
        Among(u"εμπα", 0, 1),
        Among(u"β", -1, 2),
        Among(u"δανε", -1, 1),
        Among(u"βαθυρι", -1, 2),
        Among(u"βαρκ", -1, 2),
        Among(u"μαρκ", -1, 2),
        Among(u"λ", -1, 2),
        Among(u"μ", -1, 2),
        Among(u"κορν", -1, 2),
        Among(u"αθρο", -1, 1),
        Among(u"συναθρο", 14, 1),
        Among(u"π", -1, 2),
        Among(u"ιμπ", 16, 2),
        Among(u"ρ", -1, 2),
        Among(u"μαρ", 18, 2),
        Among(u"αμπαρ", 18, 2),
        Among(u"γκρ", 18, 2),
        Among(u"βολβορ", 18, 2),
        Among(u"γλυκορ", 18, 2),
        Among(u"πιπερορ", 18, 2),
        Among(u"πρ", 18, 2),
        Among(u"μπρ", 25, 2),
        Among(u"αρρ", 18, 2),
        Among(u"γλυκυρ", 18, 2),
        Among(u"πολυρ", 18, 2),
        Among(u"λου", -1, 2)
    ]

    a_3 = [
        Among(u"ιζα", -1, 1),
        Among(u"ιζε", -1, 1),
        Among(u"ιζαμε", -1, 1),
        Among(u"ιζουμε", -1, 1),
        Among(u"ιζανε", -1, 1),
        Among(u"ιζουνε", -1, 1),
        Among(u"ιζατε", -1, 1),
        Among(u"ιζετε", -1, 1),
        Among(u"ιζει", -1, 1),
        Among(u"ιζαν", -1, 1),
        Among(u"ιζουν", -1, 1),
        Among(u"ιζεσ", -1, 1),
        Among(u"ιζεισ", -1, 1),
        Among(u"ιζω", -1, 1)
    ]

    a_4 = [
        Among(u"βι", -1, 1),
        Among(u"λι", -1, 1),
        Among(u"αλ", -1, 1),
        Among(u"εν", -1, 1),
        Among(u"σ", -1, 1),
        Among(u"χ", -1, 1),
        Among(u"υψ", -1, 1),
        Among(u"ζω", -1, 1)
    ]

    a_5 = [
        Among(u"ωθηκα", -1, 1),
        Among(u"ωθηκε", -1, 1),
        Among(u"ωθηκαμε", -1, 1),
        Among(u"ωθηκανε", -1, 1),
        Among(u"ωθηκατε", -1, 1),
        Among(u"ωθηκαν", -1, 1),
        Among(u"ωθηκεσ", -1, 1)
    ]

    a_6 = [
        Among(u"ξαναπα", -1, 1),
        Among(u"επα", -1, 1),
        Among(u"περιπα", -1, 1),
        Among(u"αναμπα", -1, 1),
        Among(u"εμπα", -1, 1),
        Among(u"χαρτοπα", -1, 1),
        Among(u"εξαρχα", -1, 1),
        Among(u"γε", -1, 2),
        Among(u"γκε", -1, 2),
        Among(u"κλε", -1, 1),
        Among(u"εκλε", 9, 1),
        Among(u"απεκλε", 10, 1),
        Among(u"αποκλε", 9, 1),
        Among(u"εσωκλε", 9, 1),
        Among(u"δανε", -1, 1),
        Among(u"πε", -1, 1),
        Among(u"επε", 15, 1),
        Among(u"μετεπε", 16, 1),
        Among(u"εσε", -1, 1),
        Among(u"γκ", -1, 2),
        Among(u"μ", -1, 2),
        Among(u"πουκαμ", 20, 2),
        Among(u"κομ", 20, 2),
        Among(u"αν", -1, 2),
        Among(u"ολο", -1, 2),
        Among(u"αθρο", -1, 1),
        Among(u"συναθρο", 25, 1),
        Among(u"π", -1, 2),
        Among(u"λαρ", -1, 2),
        Among(u"δημοκρατ", -1, 2),
        Among(u"αφ", -1, 2),
        Among(u"γιγαντοαφ", 30, 2)
    ]

    a_7 = [
        Among(u"ισα", -1, 1),
        Among(u"ισαμε", -1, 1),
        Among(u"ισανε", -1, 1),
        Among(u"ισε", -1, 1),
        Among(u"ισατε", -1, 1),
        Among(u"ισαν", -1, 1),
        Among(u"ισεσ", -1, 1)
    ]

    a_8 = [
        Among(u"ξαναπα", -1, 1),
        Among(u"επα", -1, 1),
        Among(u"περιπα", -1, 1),
        Among(u"αναμπα", -1, 1),
        Among(u"εμπα", -1, 1),
        Among(u"χαρτοπα", -1, 1),
        Among(u"εξαρχα", -1, 1),
        Among(u"κλε", -1, 1),
        Among(u"εκλε", 7, 1),
        Among(u"απεκλε", 8, 1),
        Among(u"αποκλε", 7, 1),
        Among(u"εσωκλε", 7, 1),
        Among(u"δανε", -1, 1),
        Among(u"πε", -1, 1),
        Among(u"επε", 13, 1),
        Among(u"μετεπε", 14, 1),
        Among(u"εσε", -1, 1),
        Among(u"αθρο", -1, 1),
        Among(u"συναθρο", 17, 1)
    ]

    a_9 = [
        Among(u"ισουμε", -1, 1),
        Among(u"ισουνε", -1, 1),
        Among(u"ισετε", -1, 1),
        Among(u"ισει", -1, 1),
        Among(u"ισουν", -1, 1),
        Among(u"ισεισ", -1, 1),
        Among(u"ισω", -1, 1)
    ]

    a_10 = [
        Among(u"ατα", -1, 2),
        Among(u"φα", -1, 2),
        Among(u"ηφα", 1, 2),
        Among(u"μεγ", -1, 2),
        Among(u"λυγ", -1, 2),
        Among(u"ηδ", -1, 2),
        Among(u"κλε", -1, 1),
        Among(u"εσωκλε", 6, 1),
        Among(u"πλε", -1, 1),
        Among(u"δανε", -1, 1),
        Among(u"σε", -1, 1),
        Among(u"ασε", 10, 1),
        Among(u"καθ", -1, 2),
        Among(u"εχθ", -1, 2),
        Among(u"κακ", -1, 2),
        Among(u"μακ", -1, 2),
        Among(u"σκ", -1, 2),
        Among(u"φιλ", -1, 2),
        Among(u"κυλ", -1, 2),
        Among(u"μ", -1, 2),
        Among(u"γεμ", 19, 2),
        Among(u"αχν", -1, 2),
        Among(u"συναθρο", -1, 1),
        Among(u"π", -1, 2),
        Among(u"απ", 23, 2),
        Among(u"εμπ", 23, 2),
        Among(u"ευπ", 23, 2),
        Among(u"αρ", -1, 2),
        Among(u"αορ", -1, 2),
        Among(u"γυρ", -1, 2),
        Among(u"χρ", -1, 2),
        Among(u"χωρ", -1, 2),
        Among(u"κτ", -1, 2),
        Among(u"ακτ", 32, 2),
        Among(u"χτ", -1, 2),
        Among(u"αχτ", 34, 2),
        Among(u"ταχ", -1, 2),
        Among(u"σχ", -1, 2),
        Among(u"ασχ", 37, 2),
        Among(u"υψ", -1, 2)
    ]

    a_11 = [
        Among(u"ιστα", -1, 1),
        Among(u"ιστε", -1, 1),
        Among(u"ιστη", -1, 1),
        Among(u"ιστοι", -1, 1),
        Among(u"ιστων", -1, 1),
        Among(u"ιστο", -1, 1),
        Among(u"ιστεσ", -1, 1),
        Among(u"ιστησ", -1, 1),
        Among(u"ιστοσ", -1, 1),
        Among(u"ιστουσ", -1, 1),
        Among(u"ιστου", -1, 1)
    ]

    a_12 = [
        Among(u"εγκλε", -1, 1),
        Among(u"αποκλε", -1, 1),
        Among(u"δανε", -1, 2),
        Among(u"αντιδανε", 2, 2),
        Among(u"σε", -1, 1),
        Among(u"μετασε", 4, 1),
        Among(u"μικροσε", 4, 1)
    ]

    a_13 = [
        Among(u"ατομικ", -1, 2),
        Among(u"εθνικ", -1, 4),
        Among(u"τοπικ", -1, 7),
        Among(u"εκλεκτικ", -1, 5),
        Among(u"σκεπτικ", -1, 6),
        Among(u"γνωστικ", -1, 3),
        Among(u"αγνωστικ", 5, 1),
        Among(u"αλεξανδριν", -1, 8),
        Among(u"θεατριν", -1, 10),
        Among(u"βυζαντιν", -1, 9)
    ]

    a_14 = [
        Among(u"ισμοι", -1, 1),
        Among(u"ισμων", -1, 1),
        Among(u"ισμο", -1, 1),
        Among(u"ισμοσ", -1, 1),
        Among(u"ισμουσ", -1, 1),
        Among(u"ισμου", -1, 1)
    ]

    a_15 = [
        Among(u"σ", -1, 1),
        Among(u"χ", -1, 1)
    ]

    a_16 = [
        Among(u"ουδακια", -1, 1),
        Among(u"αρακια", -1, 1),
        Among(u"ουδακι", -1, 1),
        Among(u"αρακι", -1, 1)
    ]

    a_17 = [
        Among(u"β", -1, 2),
        Among(u"βαμβ", 0, 1),
        Among(u"σλοβ", 0, 1),
        Among(u"τσεχοσλοβ", 2, 1),
        Among(u"καρδ", -1, 2),
        Among(u"ζ", -1, 2),
        Among(u"τζ", 5, 1),
        Among(u"κ", -1, 1),
        Among(u"καπακ", 7, 1),
        Among(u"σοκ", 7, 1),
        Among(u"σκ", 7, 1),
        Among(u"βαλ", -1, 2),
        Among(u"μαλ", -1, 1),
        Among(u"γλ", -1, 2),
        Among(u"τριπολ", -1, 2),
        Among(u"πλ", -1, 1),
        Among(u"λουλ", -1, 1),
        Among(u"φυλ", -1, 1),
        Among(u"καιμ", -1, 1),
        Among(u"κλιμ", -1, 1),
        Among(u"φαρμ", -1, 1),
        Among(u"γιαν", -1, 2),
        Among(u"σπαν", -1, 1),
        Among(u"ηγουμεν", -1, 2),
        Among(u"κον", -1, 1),
        Among(u"μακρυν", -1, 2),
        Among(u"π", -1, 2),
        Among(u"κατραπ", 26, 1),
        Among(u"ρ", -1, 1),
        Among(u"βρ", 28, 1),
        Among(u"λαβρ", 29, 1),
        Among(u"αμβρ", 29, 1),
        Among(u"μερ", 28, 1),
        Among(u"πατερ", 28, 2),
        Among(u"ανθρ", 28, 1),
        Among(u"κορ", 28, 1),
        Among(u"σ", -1, 1),
        Among(u"ναγκασ", 36, 1),
        Among(u"τοσ", 36, 2),
        Among(u"μουστ", -1, 1),
        Among(u"ρυ", -1, 1),
        Among(u"φ", -1, 1),
        Among(u"σφ", 41, 1),
        Among(u"αλισφ", 42, 1),
        Among(u"νυφ", 41, 2),
        Among(u"χ", -1, 1)
    ]

    a_18 = [
        Among(u"ακια", -1, 1),
        Among(u"αρακια", 0, 1),
        Among(u"ιτσα", -1, 1),
        Among(u"ακι", -1, 1),
        Among(u"αρακι", 3, 1),
        Among(u"ιτσων", -1, 1),
        Among(u"ιτσασ", -1, 1),
        Among(u"ιτσεσ", -1, 1)
    ]

    a_19 = [
        Among(u"ψαλ", -1, 1),
        Among(u"αιφν", -1, 1),
        Among(u"ολο", -1, 1),
        Among(u"ιρ", -1, 1)
    ]

    a_20 = [
        Among(u"ε", -1, 1),
        Among(u"παιχν", -1, 1)
    ]

    a_21 = [
        Among(u"ιδια", -1, 1),
        Among(u"ιδιων", -1, 1),
        Among(u"ιδιο", -1, 1)
    ]

    a_22 = [
        Among(u"ιβ", -1, 1),
        Among(u"δ", -1, 1),
        Among(u"φραγκ", -1, 1),
        Among(u"λυκ", -1, 1),
        Among(u"οβελ", -1, 1),
        Among(u"μην", -1, 1),
        Among(u"ρ", -1, 1)
    ]

    a_23 = [
        Among(u"ισκε", -1, 1),
        Among(u"ισκο", -1, 1),
        Among(u"ισκοσ", -1, 1),
        Among(u"ισκου", -1, 1)
    ]

    a_24 = [
        Among(u"αδων", -1, 1),
        Among(u"αδεσ", -1, 1)
    ]

    a_25 = [
        Among(u"γιαγι", -1, -1),
        Among(u"θει", -1, -1),
        Among(u"οκ", -1, -1),
        Among(u"μαμ", -1, -1),
        Among(u"μαν", -1, -1),
        Among(u"μπαμπ", -1, -1),
        Among(u"πεθερ", -1, -1),
        Among(u"πατερ", -1, -1),
        Among(u"κυρ", -1, -1),
        Among(u"νταντ", -1, -1)
    ]

    a_26 = [
        Among(u"εδων", -1, 1),
        Among(u"εδεσ", -1, 1)
    ]

    a_27 = [
        Among(u"μιλ", -1, 1),
        Among(u"δαπ", -1, 1),
        Among(u"γηπ", -1, 1),
        Among(u"ιπ", -1, 1),
        Among(u"εμπ", -1, 1),
        Among(u"οπ", -1, 1),
        Among(u"κρασπ", -1, 1),
        Among(u"υπ", -1, 1)
    ]

    a_28 = [
        Among(u"ουδων", -1, 1),
        Among(u"ουδεσ", -1, 1)
    ]

    a_29 = [
        Among(u"τραγ", -1, 1),
        Among(u"φε", -1, 1),
        Among(u"καλιακ", -1, 1),
        Among(u"αρκ", -1, 1),
        Among(u"σκ", -1, 1),
        Among(u"πεταλ", -1, 1),
        Among(u"βελ", -1, 1),
        Among(u"λουλ", -1, 1),
        Among(u"φλ", -1, 1),
        Among(u"χν", -1, 1),
        Among(u"πλεξ", -1, 1),
        Among(u"σπ", -1, 1),
        Among(u"φρ", -1, 1),
        Among(u"σ", -1, 1),
        Among(u"λιχ", -1, 1)
    ]

    a_30 = [
        Among(u"εων", -1, 1),
        Among(u"εωσ", -1, 1)
    ]

    a_31 = [
        Among(u"δ", -1, 1),
        Among(u"ιδ", 0, 1),
        Among(u"θ", -1, 1),
        Among(u"γαλ", -1, 1),
        Among(u"ελ", -1, 1),
        Among(u"ν", -1, 1),
        Among(u"π", -1, 1),
        Among(u"παρ", -1, 1)
    ]

    a_32 = [
        Among(u"ια", -1, 1),
        Among(u"ιων", -1, 1),
        Among(u"ιου", -1, 1)
    ]

    a_33 = [
        Among(u"ικα", -1, 1),
        Among(u"ικων", -1, 1),
        Among(u"ικο", -1, 1),
        Among(u"ικου", -1, 1)
    ]

    a_34 = [
        Among(u"αδ", -1, 1),
        Among(u"συναδ", 0, 1),
        Among(u"καταδ", 0, 1),
        Among(u"αντιδ", -1, 1),
        Among(u"ενδ", -1, 1),
        Among(u"φυλοδ", -1, 1),
        Among(u"υποδ", -1, 1),
        Among(u"πρωτοδ", -1, 1),
        Among(u"εξωδ", -1, 1),
        Among(u"ηθ", -1, 1),
        Among(u"ανηθ", 9, 1),
        Among(u"ξικ", -1, 1),
        Among(u"αλ", -1, 1),
        Among(u"αμμοχαλ", 12, 1),
        Among(u"συνομηλ", -1, 1),
        Among(u"μπολ", -1, 1),
        Among(u"μουλ", -1, 1),
        Among(u"τσαμ", -1, 1),
        Among(u"βρωμ", -1, 1),
        Among(u"αμαν", -1, 1),
        Among(u"μπαν", -1, 1),
        Among(u"καλλιν", -1, 1),
        Among(u"ποστελν", -1, 1),
        Among(u"φιλον", -1, 1),
        Among(u"καλπ", -1, 1),
        Among(u"γερ", -1, 1),
        Among(u"χασ", -1, 1),
        Among(u"μποσ", -1, 1),
        Among(u"πλιατσ", -1, 1),
        Among(u"πετσ", -1, 1),
        Among(u"πιτσ", -1, 1),
        Among(u"φυσ", -1, 1),
        Among(u"μπαγιατ", -1, 1),
        Among(u"νιτ", -1, 1),
        Among(u"πικαντ", -1, 1),
        Among(u"σερτ", -1, 1)
    ]

    a_35 = [
        Among(u"αγαμε", -1, 1),
        Among(u"ηκαμε", -1, 1),
        Among(u"ηθηκαμε", 1, 1),
        Among(u"ησαμε", -1, 1),
        Among(u"ουσαμε", -1, 1)
    ]

    a_36 = [
        Among(u"βουβ", -1, 1),
        Among(u"ξεθ", -1, 1),
        Among(u"πεθ", -1, 1),
        Among(u"αποθ", -1, 1),
        Among(u"αποκ", -1, 1),
        Among(u"ουλ", -1, 1),
        Among(u"αναπ", -1, 1),
        Among(u"πικρ", -1, 1),
        Among(u"ποτ", -1, 1),
        Among(u"αποστ", -1, 1),
        Among(u"χ", -1, 1),
        Among(u"σιχ", 10, 1)
    ]

    a_37 = [
        Among(u"τρ", -1, 1),
        Among(u"τσ", -1, 1)
    ]

    a_38 = [
        Among(u"αγανε", -1, 1),
        Among(u"ηκανε", -1, 1),
        Among(u"ηθηκανε", 1, 1),
        Among(u"ησανε", -1, 1),
        Among(u"ουσανε", -1, 1),
        Among(u"οντανε", -1, 1),
        Among(u"ιοντανε", 5, 1),
        Among(u"ουντανε", -1, 1),
        Among(u"ιουντανε", 7, 1),
        Among(u"οτανε", -1, 1),
        Among(u"ιοτανε", 9, 1)
    ]

    a_39 = [
        Among(u"ταβ", -1, 1),
        Among(u"νταβ", 0, 1),
        Among(u"ψηλοταβ", 0, 1),
        Among(u"λιβ", -1, 1),
        Among(u"κλιβ", 3, 1),
        Among(u"ξηροκλιβ", 4, 1),
        Among(u"γ", -1, 1),
        Among(u"αγ", 6, 1),
        Among(u"τραγ", 7, 1),
        Among(u"τσαγ", 7, 1),
        Among(u"αθιγγ", 6, 1),
        Among(u"τσιγγ", 6, 1),
        Among(u"ατσιγγ", 11, 1),
        Among(u"στεγ", 6, 1),
        Among(u"απηγ", 6, 1),
        Among(u"σιγ", 6, 1),
        Among(u"ανοργ", 6, 1),
        Among(u"ενοργ", 6, 1),
        Among(u"καλπουζ", -1, 1),
        Among(u"θ", -1, 1),
        Among(u"μωαμεθ", 19, 1),
        Among(u"πιθ", 19, 1),
        Among(u"απιθ", 21, 1),
        Among(u"δεκ", -1, 1),
        Among(u"πελεκ", -1, 1),
        Among(u"ικ", -1, 1),
        Among(u"ανικ", 25, 1),
        Among(u"βουλκ", -1, 1),
        Among(u"βασκ", -1, 1),
        Among(u"βραχυκ", -1, 1),
        Among(u"γαλ", -1, 1),
        Among(u"καταγαλ", 30, 1),
        Among(u"ολογαλ", 30, 1),
        Among(u"βαθυγαλ", 30, 1),
        Among(u"μελ", -1, 1),
        Among(u"καστελ", -1, 1),
        Among(u"πορτολ", -1, 1),
        Among(u"πλ", -1, 1),
        Among(u"διπλ", 37, 1),
        Among(u"λαοπλ", 37, 1),
        Among(u"ψυχοπλ", 37, 1),
        Among(u"ουλ", -1, 1),
        Among(u"μ", -1, 1),
        Among(u"ολιγοδαμ", 42, 1),
        Among(u"μουσουλμ", 42, 1),
        Among(u"δραδουμ", 42, 1),
        Among(u"βραχμ", 42, 1),
        Among(u"ν", -1, 1),
        Among(u"αμερικαν", 47, 1),
        Among(u"π", -1, 1),
        Among(u"αδαπ", 49, 1),
        Among(u"χαμηλοδαπ", 49, 1),
        Among(u"πολυδαπ", 49, 1),
        Among(u"κοπ", 49, 1),
        Among(u"υποκοπ", 53, 1),
        Among(u"τσοπ", 49, 1),
        Among(u"σπ", 49, 1),
        Among(u"ερ", -1, 1),
        Among(u"γερ", 57, 1),
        Among(u"βετερ", 57, 1),
        Among(u"λουθηρ", -1, 1),
        Among(u"κορμορ", -1, 1),
        Among(u"περιτρ", -1, 1),
        Among(u"ουρ", -1, 1),
        Among(u"σ", -1, 1),
        Among(u"βασ", 64, 1),
        Among(u"πολισ", 64, 1),
        Among(u"σαρακατσ", 64, 1),
        Among(u"θυσ", 64, 1),
        Among(u"διατ", -1, 1),
        Among(u"πλατ", -1, 1),
        Among(u"τσαρλατ", -1, 1),
        Among(u"τετ", -1, 1),
        Among(u"πουριτ", -1, 1),
        Among(u"σουλτ", -1, 1),
        Among(u"μαιντ", -1, 1),
        Among(u"ζωντ", -1, 1),
        Among(u"καστ", -1, 1),
        Among(u"φ", -1, 1),
        Among(u"διαφ", 78, 1),
        Among(u"στεφ", 78, 1),
        Among(u"φωτοστεφ", 80, 1),
        Among(u"περηφ", 78, 1),
        Among(u"υπερηφ", 82, 1),
        Among(u"κοιλαρφ", 78, 1),
        Among(u"πενταρφ", 78, 1),
        Among(u"ορφ", 78, 1),
        Among(u"χ", -1, 1),
        Among(u"αμηχ", 87, 1),
        Among(u"βιομηχ", 87, 1),
        Among(u"μεγλοβιομηχ", 89, 1),
        Among(u"καπνοβιομηχ", 89, 1),
        Among(u"μικροβιομηχ", 89, 1),
        Among(u"πολυμηχ", 87, 1),
        Among(u"λιχ", 87, 1)
    ]

    a_40 = [
        Among(u"ησετε", -1, 1)
    ]

    a_41 = [
        Among(u"ενδ", -1, 1),
        Among(u"συνδ", -1, 1),
        Among(u"οδ", -1, 1),
        Among(u"διαθ", -1, 1),
        Among(u"καθ", -1, 1),
        Among(u"ραθ", -1, 1),
        Among(u"ταθ", -1, 1),
        Among(u"τιθ", -1, 1),
        Among(u"εκθ", -1, 1),
        Among(u"ενθ", -1, 1),
        Among(u"συνθ", -1, 1),
        Among(u"ροθ", -1, 1),
        Among(u"υπερθ", -1, 1),
        Among(u"σθ", -1, 1),
        Among(u"ευθ", -1, 1),
        Among(u"αρκ", -1, 1),
        Among(u"ωφελ", -1, 1),
        Among(u"βολ", -1, 1),
        Among(u"αιν", -1, 1),
        Among(u"πον", -1, 1),
        Among(u"ρον", -1, 1),
        Among(u"συν", -1, 1),
        Among(u"βαρ", -1, 1),
        Among(u"βρ", -1, 1),
        Among(u"αιρ", -1, 1),
        Among(u"φορ", -1, 1),
        Among(u"ευρ", -1, 1),
        Among(u"πυρ", -1, 1),
        Among(u"χωρ", -1, 1),
        Among(u"νετ", -1, 1),
        Among(u"σχ", -1, 1)
    ]

    a_42 = [
        Among(u"παγ", -1, 1),
        Among(u"δ", -1, 1),
        Among(u"αδ", 1, 1),
        Among(u"θ", -1, 1),
        Among(u"αθ", 3, 1),
        Among(u"τοκ", -1, 1),
        Among(u"σκ", -1, 1),
        Among(u"παρακαλ", -1, 1),
        Among(u"σκελ", -1, 1),
        Among(u"απλ", -1, 1),
        Among(u"εμ", -1, 1),
        Among(u"αν", -1, 1),
        Among(u"βεν", -1, 1),
        Among(u"βαρον", -1, 1),
        Among(u"κοπ", -1, 1),
        Among(u"σερπ", -1, 1),
        Among(u"αβαρ", -1, 1),
        Among(u"εναρ", -1, 1),
        Among(u"αβρ", -1, 1),
        Among(u"μπορ", -1, 1),
        Among(u"θαρρ", -1, 1),
        Among(u"ντρ", -1, 1),
        Among(u"υ", -1, 1),
        Among(u"νιφ", -1, 1),
        Among(u"συρφ", -1, 1)
    ]

    a_43 = [
        Among(u"οντασ", -1, 1),
        Among(u"ωντασ", -1, 1)
    ]

    a_44 = [
        Among(u"ομαστε", -1, 1),
        Among(u"ιομαστε", 0, 1)
    ]

    a_45 = [
        Among(u"π", -1, 1),
        Among(u"απ", 0, 1),
        Among(u"ακαταπ", 1, 1),
        Among(u"συμπ", 0, 1),
        Among(u"ασυμπ", 3, 1),
        Among(u"αμεταμφ", -1, 1)
    ]

    a_46 = [
        Among(u"ζ", -1, 1),
        Among(u"αλ", -1, 1),
        Among(u"παρακαλ", 1, 1),
        Among(u"εκτελ", -1, 1),
        Among(u"μ", -1, 1),
        Among(u"ξ", -1, 1),
        Among(u"προ", -1, 1),
        Among(u"αρ", -1, 1),
        Among(u"νισ", -1, 1)
    ]

    a_47 = [
        Among(u"ηθηκα", -1, 1),
        Among(u"ηθηκε", -1, 1),
        Among(u"ηθηκεσ", -1, 1)
    ]

    a_48 = [
        Among(u"πιθ", -1, 1),
        Among(u"οθ", -1, 1),
        Among(u"ναρθ", -1, 1),
        Among(u"σκουλ", -1, 1),
        Among(u"σκωλ", -1, 1),
        Among(u"σφ", -1, 1)
    ]

    a_49 = [
        Among(u"θ", -1, 1),
        Among(u"διαθ", 0, 1),
        Among(u"παρακαταθ", 0, 1),
        Among(u"συνθ", 0, 1),
        Among(u"προσθ", 0, 1)
    ]

    a_50 = [
        Among(u"ηκα", -1, 1),
        Among(u"ηκε", -1, 1),
        Among(u"ηκεσ", -1, 1)
    ]

    a_51 = [
        Among(u"φαγ", -1, 1),
        Among(u"ληγ", -1, 1),
        Among(u"φρυδ", -1, 1),
        Among(u"μαντιλ", -1, 1),
        Among(u"μαλλ", -1, 1),
        Among(u"ομ", -1, 1),
        Among(u"βλεπ", -1, 1),
        Among(u"ποδαρ", -1, 1),
        Among(u"κυματ", -1, 1),
        Among(u"πρωτ", -1, 1),
        Among(u"λαχ", -1, 1),
        Among(u"πανταχ", -1, 1)
    ]

    a_52 = [
        Among(u"τσα", -1, 1),
        Among(u"χαδ", -1, 1),
        Among(u"μεδ", -1, 1),
        Among(u"λαμπιδ", -1, 1),
        Among(u"δε", -1, 1),
        Among(u"πλε", -1, 1),
        Among(u"μεσαζ", -1, 1),
        Among(u"δεσποζ", -1, 1),
        Among(u"αιθ", -1, 1),
        Among(u"φαρμακ", -1, 1),
        Among(u"αγκ", -1, 1),
        Among(u"ανηκ", -1, 1),
        Among(u"λ", -1, 1),
        Among(u"μ", -1, 1),
        Among(u"αμ", 13, 1),
        Among(u"βρομ", 13, 1),
        Among(u"υποτειν", -1, 1),
        Among(u"εκλιπ", -1, 1),
        Among(u"ρ", -1, 1),
        Among(u"ενδιαφερ", 18, 1),
        Among(u"αναρρ", 18, 1),
        Among(u"πατ", -1, 1),
        Among(u"καθαρευ", -1, 1),
        Among(u"δευτερευ", -1, 1),
        Among(u"λεχ", -1, 1)
    ]

    a_53 = [
        Among(u"ουσα", -1, 1),
        Among(u"ουσε", -1, 1),
        Among(u"ουσεσ", -1, 1)
    ]

    a_54 = [
        Among(u"πελ", -1, 1),
        Among(u"λλ", -1, 1),
        Among(u"σμην", -1, 1),
        Among(u"ρπ", -1, 1),
        Among(u"πρ", -1, 1),
        Among(u"φρ", -1, 1),
        Among(u"χορτ", -1, 1),
        Among(u"οφ", -1, 1),
        Among(u"ψοφ", 7, -1),
        Among(u"σφ", -1, 1),
        Among(u"λοχ", -1, 1),
        Among(u"ναυλοχ", 10, -1)
    ]

    a_55 = [
        Among(u"αμαλλι", -1, 1),
        Among(u"λ", -1, 1),
        Among(u"αμαλ", 1, 1),
        Among(u"μ", -1, 1),
        Among(u"ουλαμ", 3, 1),
        Among(u"εν", -1, 1),
        Among(u"δερβεν", 5, 1),
        Among(u"π", -1, 1),
        Among(u"αειπ", 7, 1),
        Among(u"αρτιπ", 7, 1),
        Among(u"συμπ", 7, 1),
        Among(u"νεοπ", 7, 1),
        Among(u"κροκαλοπ", 7, 1),
        Among(u"ολοπ", 7, 1),
        Among(u"προσωποπ", 7, 1),
        Among(u"σιδηροπ", 7, 1),
        Among(u"δροσοπ", 7, 1),
        Among(u"ασπ", 7, 1),
        Among(u"ανυπ", 7, 1),
        Among(u"ρ", -1, 1),
        Among(u"ασπαρ", 19, 1),
        Among(u"χαρ", 19, 1),
        Among(u"αχαρ", 21, 1),
        Among(u"απερ", 19, 1),
        Among(u"τρ", 19, 1),
        Among(u"ουρ", 19, 1),
        Among(u"τ", -1, 1),
        Among(u"διατ", 26, 1),
        Among(u"επιτ", 26, 1),
        Among(u"συντ", 26, 1),
        Among(u"ομοτ", 26, 1),
        Among(u"νομοτ", 30, 1),
        Among(u"αποτ", 26, 1),
        Among(u"υποτ", 26, 1),
        Among(u"αβαστ", 26, 1),
        Among(u"αιμοστ", 26, 1),
        Among(u"προστ", 26, 1),
        Among(u"ανυστ", 26, 1),
        Among(u"ναυ", -1, 1),
        Among(u"αφ", -1, 1),
        Among(u"ξεφ", -1, 1),
        Among(u"αδηφ", -1, 1),
        Among(u"παμφ", -1, 1),
        Among(u"πολυφ", -1, 1)
    ]

    a_56 = [
        Among(u"αγα", -1, 1),
        Among(u"αγε", -1, 1),
        Among(u"αγεσ", -1, 1)
    ]

    a_57 = [
        Among(u"ησα", -1, 1),
        Among(u"ησε", -1, 1),
        Among(u"ησου", -1, 1)
    ]

    a_58 = [
        Among(u"ν", -1, 1),
        Among(u"δωδεκαν", 0, 1),
        Among(u"επταν", 0, 1),
        Among(u"μεγαλον", 0, 1),
        Among(u"ερημον", 0, 1),
        Among(u"χερσον", 0, 1)
    ]

    a_59 = [
        Among(u"ηστε", -1, 1)
    ]

    a_60 = [
        Among(u"σβ", -1, 1),
        Among(u"ασβ", 0, 1),
        Among(u"απλ", -1, 1),
        Among(u"αειμν", -1, 1),
        Among(u"χρ", -1, 1),
        Among(u"αχρ", 4, 1),
        Among(u"κοινοχρ", 4, 1),
        Among(u"δυσχρ", 4, 1),
        Among(u"ευχρ", 4, 1),
        Among(u"παλιμψ", -1, 1)
    ]

    a_61 = [
        Among(u"ουνε", -1, 1),
        Among(u"ηθουνε", 0, 1),
        Among(u"ησουνε", 0, 1)
    ]

    a_62 = [
        Among(u"σπι", -1, 1),
        Among(u"ν", -1, 1),
        Among(u"εξων", 1, 1),
        Among(u"ρ", -1, 1),
        Among(u"στραβομουτσ", -1, 1),
        Among(u"κακομουτσ", -1, 1)
    ]

    a_63 = [
        Among(u"ουμε", -1, 1),
        Among(u"ηθουμε", 0, 1),
        Among(u"ησουμε", 0, 1)
    ]

    a_64 = [
        Among(u"αζ", -1, 1),
        Among(u"ωριοπλ", -1, 1),
        Among(u"ασουσ", -1, 1),
        Among(u"παρασουσ", 2, 1),
        Among(u"αλλοσουσ", -1, 1),
        Among(u"φ", -1, 1),
        Among(u"χ", -1, 1)
    ]

    a_65 = [
        Among(u"ματα", -1, 1),
        Among(u"ματων", -1, 1),
        Among(u"ματοσ", -1, 1)
    ]

    a_66 = [
        Among(u"α", -1, 1),
        Among(u"ιουμα", 0, 1),
        Among(u"ομουνα", 0, 1),
        Among(u"ιομουνα", 2, 1),
        Among(u"οσουνα", 0, 1),
        Among(u"ιοσουνα", 4, 1),
        Among(u"ε", -1, 1),
        Among(u"αγατε", 6, 1),
        Among(u"ηκατε", 6, 1),
        Among(u"ηθηκατε", 8, 1),
        Among(u"ησατε", 6, 1),
        Among(u"ουσατε", 6, 1),
        Among(u"ειτε", 6, 1),
        Among(u"ηθειτε", 12, 1),
        Among(u"ιεμαστε", 6, 1),
        Among(u"ουμαστε", 6, 1),
        Among(u"ιουμαστε", 15, 1),
        Among(u"ιεσαστε", 6, 1),
        Among(u"οσαστε", 6, 1),
        Among(u"ιοσαστε", 18, 1),
        Among(u"η", -1, 1),
        Among(u"ι", -1, 1),
        Among(u"αμαι", 21, 1),
        Among(u"ιεμαι", 21, 1),
        Among(u"ομαι", 21, 1),
        Among(u"ουμαι", 21, 1),
        Among(u"ασαι", 21, 1),
        Among(u"εσαι", 21, 1),
        Among(u"ιεσαι", 27, 1),
        Among(u"αται", 21, 1),
        Among(u"εται", 21, 1),
        Among(u"ιεται", 30, 1),
        Among(u"ονται", 21, 1),
        Among(u"ουνται", 21, 1),
        Among(u"ιουνται", 33, 1),
        Among(u"ει", 21, 1),
        Among(u"αει", 35, 1),
        Among(u"ηθει", 35, 1),
        Among(u"ησει", 35, 1),
        Among(u"οι", 21, 1),
        Among(u"αν", -1, 1),
        Among(u"αγαν", 40, 1),
        Among(u"ηκαν", 40, 1),
        Among(u"ηθηκαν", 42, 1),
        Among(u"ησαν", 40, 1),
        Among(u"ουσαν", 40, 1),
        Among(u"οντουσαν", 45, 1),
        Among(u"ιοντουσαν", 46, 1),
        Among(u"ονταν", 40, 1),
        Among(u"ιονταν", 48, 1),
        Among(u"ουνταν", 40, 1),
        Among(u"ιουνταν", 50, 1),
        Among(u"οταν", 40, 1),
        Among(u"ιοταν", 52, 1),
        Among(u"ομασταν", 40, 1),
        Among(u"ιομασταν", 54, 1),
        Among(u"οσασταν", 40, 1),
        Among(u"ιοσασταν", 56, 1),
        Among(u"ουν", -1, 1),
        Among(u"ηθουν", 58, 1),
        Among(u"ομουν", 58, 1),
        Among(u"ιομουν", 60, 1),
        Among(u"ησουν", 58, 1),
        Among(u"οσουν", 58, 1),
        Among(u"ιοσουν", 63, 1),
        Among(u"ων", -1, 1),
        Among(u"ηδων", 65, 1),
        Among(u"ο", -1, 1),
        Among(u"ασ", -1, 1),
        Among(u"εσ", -1, 1),
        Among(u"ηδεσ", 69, 1),
        Among(u"ησεσ", 69, 1),
        Among(u"ησ", -1, 1),
        Among(u"εισ", -1, 1),
        Among(u"ηθεισ", 73, 1),
        Among(u"οσ", -1, 1),
        Among(u"υσ", -1, 1),
        Among(u"ουσ", 76, 1),
        Among(u"υ", -1, 1),
        Among(u"ου", 78, 1),
        Among(u"ω", -1, 1),
        Among(u"αω", 80, 1),
        Among(u"ηθω", 80, 1),
        Among(u"ησω", 80, 1)
    ]

    a_67 = [
        Among(u"οτερ", -1, 1),
        Among(u"εστερ", -1, 1),
        Among(u"υτερ", -1, 1),
        Among(u"ωτερ", -1, 1),
        Among(u"οτατ", -1, 1),
        Among(u"εστατ", -1, 1),
        Among(u"υτατ", -1, 1),
        Among(u"ωτατ", -1, 1)
    ]


class lab0(BaseException): pass


class lab1(BaseException): pass


class lab2(BaseException): pass


class lab3(BaseException): pass

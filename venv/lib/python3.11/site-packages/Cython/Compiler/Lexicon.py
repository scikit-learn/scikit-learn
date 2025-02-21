# -*- coding: utf-8 -*-
# cython: language_level=3, py2_import=True
#
#   Cython Scanner - Lexical Definitions
#

from __future__ import absolute_import, unicode_literals

raw_prefixes = "rR"
bytes_prefixes = "bB"
string_prefixes = "fFuU" + bytes_prefixes
char_prefixes = "cC"
any_string_prefix = raw_prefixes + string_prefixes + char_prefixes
IDENT = 'IDENT'


def make_lexicon():
    from ..Plex import \
        Str, Any, AnyBut, AnyChar, Rep, Rep1, Opt, Bol, Eol, Eof, \
        TEXT, IGNORE, Method, State, Lexicon, Range

    nonzero_digit = Any("123456789")
    digit = Any("0123456789")
    bindigit = Any("01")
    octdigit = Any("01234567")
    hexdigit = Any("0123456789ABCDEFabcdef")
    indentation = Bol + Rep(Any(" \t"))

    # The list of valid unicode identifier characters are pretty slow to generate at runtime,
    # and require Python3, so are just included directly here
    # (via the generated code block at the bottom of the file)
    unicode_start_character = (Any(unicode_start_ch_any) | Range(unicode_start_ch_range))
    unicode_continuation_character = (
        unicode_start_character |
        Any(unicode_continuation_ch_any) | Range(unicode_continuation_ch_range))

    def underscore_digits(d):
        return Rep1(d) + Rep(Str("_") + Rep1(d))

    def prefixed_digits(prefix, digits):
        return prefix + Opt(Str("_")) + underscore_digits(digits)

    decimal = underscore_digits(digit)
    dot = Str(".")
    exponent = Any("Ee") + Opt(Any("+-")) + decimal
    decimal_fract = (decimal + dot + Opt(decimal)) | (dot + decimal)

    #name = letter + Rep(letter | digit)
    name = unicode_start_character + Rep(unicode_continuation_character)
    intconst = (prefixed_digits(nonzero_digit, digit) |  # decimal literals with underscores must not start with '0'
                (Str("0") + (prefixed_digits(Any("Xx"), hexdigit) |
                             prefixed_digits(Any("Oo"), octdigit) |
                             prefixed_digits(Any("Bb"), bindigit) )) |
                underscore_digits(Str('0'))  # 0_0_0_0... is allowed as a decimal literal
                | Rep1(digit)  # FIXME: remove these Py2 style decimal/octal literals (PY_VERSION_HEX < 3)
                )
    intsuffix = (Opt(Any("Uu")) + Opt(Any("Ll")) + Opt(Any("Ll"))) | (Opt(Any("Ll")) + Opt(Any("Ll")) + Opt(Any("Uu")))
    intliteral = intconst + intsuffix
    fltconst = (decimal_fract + Opt(exponent)) | (decimal + exponent)
    imagconst = (intconst | fltconst) + Any("jJ")

    # invalid combinations of prefixes are caught in p_string_literal
    beginstring = Opt(Rep(Any(string_prefixes + raw_prefixes)) |
                      Any(char_prefixes)
                      ) + (Str("'") | Str('"') | Str("'''") | Str('"""'))
    two_oct = octdigit + octdigit
    three_oct = octdigit + octdigit + octdigit
    two_hex = hexdigit + hexdigit
    four_hex = two_hex + two_hex
    escapeseq = Str("\\") + (two_oct | three_oct |
                             Str('N{') + Rep(AnyBut('}')) + Str('}') |
                             Str('u') + four_hex | Str('x') + two_hex |
                             Str('U') + four_hex + four_hex | AnyChar)

    bra = Any("([{")
    ket = Any(")]}")
    ellipsis = Str("...")
    punct = Any(":,;+-*/|&<>=.%`~^?!@")
    diphthong = Str("==", "<>", "!=", "<=", ">=", "<<", ">>", "**", "//",
                    "+=", "-=", "*=", "/=", "%=", "|=", "^=", "&=",
                    "<<=", ">>=", "**=", "//=", "->", "@=", "&&", "||", ':=')
    spaces = Rep1(Any(" \t\f"))
    escaped_newline = Str("\\\n")
    lineterm = Eol + Opt(Str("\n"))

    comment = Str("#") + Rep(AnyBut("\n"))

    return Lexicon([
        (name, Method('normalize_ident')),
        (intliteral, Method('strip_underscores', symbol='INT')),
        (fltconst, Method('strip_underscores', symbol='FLOAT')),
        (imagconst, Method('strip_underscores', symbol='IMAG')),
        (ellipsis | punct | diphthong, TEXT),

        (bra, Method('open_bracket_action')),
        (ket, Method('close_bracket_action')),
        (lineterm, Method('newline_action')),

        (beginstring, Method('begin_string_action')),

        (comment, IGNORE),
        (spaces, IGNORE),
        (escaped_newline, IGNORE),

        State('INDENT', [
            (comment + lineterm, Method('commentline')),
            (Opt(spaces) + Opt(comment) + lineterm, IGNORE),
            (indentation, Method('indentation_action')),
            (Eof, Method('eof_action'))
        ]),

        State('SQ_STRING', [
            (escapeseq, 'ESCAPE'),
            (Rep1(AnyBut("'\"\n\\")), 'CHARS'),
            (Str('"'), 'CHARS'),
            (Str("\n"), Method('unclosed_string_action')),
            (Str("'"), Method('end_string_action')),
            (Eof, 'EOF')
        ]),

        State('DQ_STRING', [
            (escapeseq, 'ESCAPE'),
            (Rep1(AnyBut('"\n\\')), 'CHARS'),
            (Str("'"), 'CHARS'),
            (Str("\n"), Method('unclosed_string_action')),
            (Str('"'), Method('end_string_action')),
            (Eof, 'EOF')
        ]),

        State('TSQ_STRING', [
            (escapeseq, 'ESCAPE'),
            (Rep1(AnyBut("'\"\n\\")), 'CHARS'),
            (Any("'\""), 'CHARS'),
            (Str("\n"), 'NEWLINE'),
            (Str("'''"), Method('end_string_action')),
            (Eof, 'EOF')
        ]),

        State('TDQ_STRING', [
            (escapeseq, 'ESCAPE'),
            (Rep1(AnyBut('"\'\n\\')), 'CHARS'),
            (Any("'\""), 'CHARS'),
            (Str("\n"), 'NEWLINE'),
            (Str('"""'), Method('end_string_action')),
            (Eof, 'EOF')
        ]),

        (Eof, Method('eof_action'))
        ],

        # FIXME: Plex 1.9 needs different args here from Plex 1.1.4
        #debug_flags = scanner_debug_flags,
        #debug_file = scanner_dump_file
        )


# BEGIN GENERATED CODE
# Generated with 'cython-generate-lexicon.py' from:
# cpython 3.12.0a7+ (heads/master:4cd1cc843a, Apr 11 2023, 10:32:26) [GCC 11.3.0]

unicode_start_ch_any = (
    u"\u005f\u00aa\u00b5\u00ba\u02ec\u02ee\u037f\u0386\u038c\u0559\u06d5"
    u"\u06ff\u0710\u07b1\u07fa\u081a\u0824\u0828\u093d\u0950\u09b2\u09bd"
    u"\u09ce\u09fc\u0a5e\u0abd\u0ad0\u0af9\u0b3d\u0b71\u0b83\u0b9c\u0bd0"
    u"\u0c3d\u0c5d\u0c80\u0cbd\u0d3d\u0d4e\u0dbd\u0e32\u0e84\u0ea5\u0eb2"
    u"\u0ebd\u0ec6\u0f00\u103f\u1061\u108e\u10c7\u10cd\u1258\u12c0\u17d7"
    u"\u17dc\u18aa\u1aa7\u1cfa\u1f59\u1f5b\u1f5d\u1fbe\u2071\u207f\u2102"
    u"\u2107\u2115\u2124\u2126\u2128\u214e\u2d27\u2d2d\u2d6f\ua7d3\ua8fb"
    u"\ua9cf\uaa7a\uaab1\uaac0\uaac2\ufb1d\ufb3e\ufe71\ufe73\ufe77\ufe79"
    u"\ufe7b\ufe7d\U00010808\U0001083c\U00010a00\U00010f27\U00011075\U00011144\U00011147\U00011176\U000111da"
    u"\U000111dc\U00011288\U0001133d\U00011350\U000114c7\U00011644\U000116b8\U00011909\U0001193f\U00011941\U000119e1"
    u"\U000119e3\U00011a00\U00011a3a\U00011a50\U00011a9d\U00011c40\U00011d46\U00011d98\U00011f02\U00011fb0\U00016f50"
    u"\U00016fe3\U0001b132\U0001b155\U0001d4a2\U0001d4bb\U0001d546\U0001e14e\U0001e94b\U0001ee24\U0001ee27\U0001ee39"
    u"\U0001ee3b\U0001ee42\U0001ee47\U0001ee49\U0001ee4b\U0001ee54\U0001ee57\U0001ee59\U0001ee5b\U0001ee5d\U0001ee5f"
    u"\U0001ee64\U0001ee7e"
)
unicode_start_ch_range = (
    u"\u0041\u005a\u0061\u007a\u00c0\u00d6\u00d8\u00f6\u00f8\u02c1\u02c6"
    u"\u02d1\u02e0\u02e4\u0370\u0374\u0376\u0377\u037b\u037d\u0388\u038a"
    u"\u038e\u03a1\u03a3\u03f5\u03f7\u0481\u048a\u052f\u0531\u0556\u0560"
    u"\u0588\u05d0\u05ea\u05ef\u05f2\u0620\u064a\u066e\u066f\u0671\u06d3"
    u"\u06e5\u06e6\u06ee\u06ef\u06fa\u06fc\u0712\u072f\u074d\u07a5\u07ca"
    u"\u07ea\u07f4\u07f5\u0800\u0815\u0840\u0858\u0860\u086a\u0870\u0887"
    u"\u0889\u088e\u08a0\u08c9\u0904\u0939\u0958\u0961\u0971\u0980\u0985"
    u"\u098c\u098f\u0990\u0993\u09a8\u09aa\u09b0\u09b6\u09b9\u09dc\u09dd"
    u"\u09df\u09e1\u09f0\u09f1\u0a05\u0a0a\u0a0f\u0a10\u0a13\u0a28\u0a2a"
    u"\u0a30\u0a32\u0a33\u0a35\u0a36\u0a38\u0a39\u0a59\u0a5c\u0a72\u0a74"
    u"\u0a85\u0a8d\u0a8f\u0a91\u0a93\u0aa8\u0aaa\u0ab0\u0ab2\u0ab3\u0ab5"
    u"\u0ab9\u0ae0\u0ae1\u0b05\u0b0c\u0b0f\u0b10\u0b13\u0b28\u0b2a\u0b30"
    u"\u0b32\u0b33\u0b35\u0b39\u0b5c\u0b5d\u0b5f\u0b61\u0b85\u0b8a\u0b8e"
    u"\u0b90\u0b92\u0b95\u0b99\u0b9a\u0b9e\u0b9f\u0ba3\u0ba4\u0ba8\u0baa"
    u"\u0bae\u0bb9\u0c05\u0c0c\u0c0e\u0c10\u0c12\u0c28\u0c2a\u0c39\u0c58"
    u"\u0c5a\u0c60\u0c61\u0c85\u0c8c\u0c8e\u0c90\u0c92\u0ca8\u0caa\u0cb3"
    u"\u0cb5\u0cb9\u0cdd\u0cde\u0ce0\u0ce1\u0cf1\u0cf2\u0d04\u0d0c\u0d0e"
    u"\u0d10\u0d12\u0d3a\u0d54\u0d56\u0d5f\u0d61\u0d7a\u0d7f\u0d85\u0d96"
    u"\u0d9a\u0db1\u0db3\u0dbb\u0dc0\u0dc6\u0e01\u0e30\u0e40\u0e46\u0e81"
    u"\u0e82\u0e86\u0e8a\u0e8c\u0ea3\u0ea7\u0eb0\u0ec0\u0ec4\u0edc\u0edf"
    u"\u0f40\u0f47\u0f49\u0f6c\u0f88\u0f8c\u1000\u102a\u1050\u1055\u105a"
    u"\u105d\u1065\u1066\u106e\u1070\u1075\u1081\u10a0\u10c5\u10d0\u10fa"
    u"\u10fc\u1248\u124a\u124d\u1250\u1256\u125a\u125d\u1260\u1288\u128a"
    u"\u128d\u1290\u12b0\u12b2\u12b5\u12b8\u12be\u12c2\u12c5\u12c8\u12d6"
    u"\u12d8\u1310\u1312\u1315\u1318\u135a\u1380\u138f\u13a0\u13f5\u13f8"
    u"\u13fd\u1401\u166c\u166f\u167f\u1681\u169a\u16a0\u16ea\u16ee\u16f8"
    u"\u1700\u1711\u171f\u1731\u1740\u1751\u1760\u176c\u176e\u1770\u1780"
    u"\u17b3\u1820\u1878\u1880\u18a8\u18b0\u18f5\u1900\u191e\u1950\u196d"
    u"\u1970\u1974\u1980\u19ab\u19b0\u19c9\u1a00\u1a16\u1a20\u1a54\u1b05"
    u"\u1b33\u1b45\u1b4c\u1b83\u1ba0\u1bae\u1baf\u1bba\u1be5\u1c00\u1c23"
    u"\u1c4d\u1c4f\u1c5a\u1c7d\u1c80\u1c88\u1c90\u1cba\u1cbd\u1cbf\u1ce9"
    u"\u1cec\u1cee\u1cf3\u1cf5\u1cf6\u1d00\u1dbf\u1e00\u1f15\u1f18\u1f1d"
    u"\u1f20\u1f45\u1f48\u1f4d\u1f50\u1f57\u1f5f\u1f7d\u1f80\u1fb4\u1fb6"
    u"\u1fbc\u1fc2\u1fc4\u1fc6\u1fcc\u1fd0\u1fd3\u1fd6\u1fdb\u1fe0\u1fec"
    u"\u1ff2\u1ff4\u1ff6\u1ffc\u2090\u209c\u210a\u2113\u2118\u211d\u212a"
    u"\u2139\u213c\u213f\u2145\u2149\u2160\u2188\u2c00\u2ce4\u2ceb\u2cee"
    u"\u2cf2\u2cf3\u2d00\u2d25\u2d30\u2d67\u2d80\u2d96\u2da0\u2da6\u2da8"
    u"\u2dae\u2db0\u2db6\u2db8\u2dbe\u2dc0\u2dc6\u2dc8\u2dce\u2dd0\u2dd6"
    u"\u2dd8\u2dde\u3005\u3007\u3021\u3029\u3031\u3035\u3038\u303c\u3041"
    u"\u3096\u309d\u309f\u30a1\u30fa\u30fc\u30ff\u3105\u312f\u3131\u318e"
    u"\u31a0\u31bf\u31f0\u31ff\u3400\u4dbf\u4e00\ua48c\ua4d0\ua4fd\ua500"
    u"\ua60c\ua610\ua61f\ua62a\ua62b\ua640\ua66e\ua67f\ua69d\ua6a0\ua6ef"
    u"\ua717\ua71f\ua722\ua788\ua78b\ua7ca\ua7d0\ua7d1\ua7d5\ua7d9\ua7f2"
    u"\ua801\ua803\ua805\ua807\ua80a\ua80c\ua822\ua840\ua873\ua882\ua8b3"
    u"\ua8f2\ua8f7\ua8fd\ua8fe\ua90a\ua925\ua930\ua946\ua960\ua97c\ua984"
    u"\ua9b2\ua9e0\ua9e4\ua9e6\ua9ef\ua9fa\ua9fe\uaa00\uaa28\uaa40\uaa42"
    u"\uaa44\uaa4b\uaa60\uaa76\uaa7e\uaaaf\uaab5\uaab6\uaab9\uaabd\uaadb"
    u"\uaadd\uaae0\uaaea\uaaf2\uaaf4\uab01\uab06\uab09\uab0e\uab11\uab16"
    u"\uab20\uab26\uab28\uab2e\uab30\uab5a\uab5c\uab69\uab70\uabe2\uac00"
    u"\ud7a3\ud7b0\ud7c6\ud7cb\ud7fb\uf900\ufa6d\ufa70\ufad9\ufb00\ufb06"
    u"\ufb13\ufb17\ufb1f\ufb28\ufb2a\ufb36\ufb38\ufb3c\ufb40\ufb41\ufb43"
    u"\ufb44\ufb46\ufbb1\ufbd3\ufc5d\ufc64\ufd3d\ufd50\ufd8f\ufd92\ufdc7"
    u"\ufdf0\ufdf9\ufe7f\ufefc\uff21\uff3a\uff41\uff5a\uff66\uff9d\uffa0"
    u"\uffbe\uffc2\uffc7\uffca\uffcf\uffd2\uffd7\uffda\uffdc\U00010000\U0001000b"
    u"\U0001000d\U00010026\U00010028\U0001003a\U0001003c\U0001003d\U0001003f\U0001004d\U00010050\U0001005d\U00010080"
    u"\U000100fa\U00010140\U00010174\U00010280\U0001029c\U000102a0\U000102d0\U00010300\U0001031f\U0001032d\U0001034a"
    u"\U00010350\U00010375\U00010380\U0001039d\U000103a0\U000103c3\U000103c8\U000103cf\U000103d1\U000103d5\U00010400"
    u"\U0001049d\U000104b0\U000104d3\U000104d8\U000104fb\U00010500\U00010527\U00010530\U00010563\U00010570\U0001057a"
    u"\U0001057c\U0001058a\U0001058c\U00010592\U00010594\U00010595\U00010597\U000105a1\U000105a3\U000105b1\U000105b3"
    u"\U000105b9\U000105bb\U000105bc\U00010600\U00010736\U00010740\U00010755\U00010760\U00010767\U00010780\U00010785"
    u"\U00010787\U000107b0\U000107b2\U000107ba\U00010800\U00010805\U0001080a\U00010835\U00010837\U00010838\U0001083f"
    u"\U00010855\U00010860\U00010876\U00010880\U0001089e\U000108e0\U000108f2\U000108f4\U000108f5\U00010900\U00010915"
    u"\U00010920\U00010939\U00010980\U000109b7\U000109be\U000109bf\U00010a10\U00010a13\U00010a15\U00010a17\U00010a19"
    u"\U00010a35\U00010a60\U00010a7c\U00010a80\U00010a9c\U00010ac0\U00010ac7\U00010ac9\U00010ae4\U00010b00\U00010b35"
    u"\U00010b40\U00010b55\U00010b60\U00010b72\U00010b80\U00010b91\U00010c00\U00010c48\U00010c80\U00010cb2\U00010cc0"
    u"\U00010cf2\U00010d00\U00010d23\U00010e80\U00010ea9\U00010eb0\U00010eb1\U00010f00\U00010f1c\U00010f30\U00010f45"
    u"\U00010f70\U00010f81\U00010fb0\U00010fc4\U00010fe0\U00010ff6\U00011003\U00011037\U00011071\U00011072\U00011083"
    u"\U000110af\U000110d0\U000110e8\U00011103\U00011126\U00011150\U00011172\U00011183\U000111b2\U000111c1\U000111c4"
    u"\U00011200\U00011211\U00011213\U0001122b\U0001123f\U00011240\U00011280\U00011286\U0001128a\U0001128d\U0001128f"
    u"\U0001129d\U0001129f\U000112a8\U000112b0\U000112de\U00011305\U0001130c\U0001130f\U00011310\U00011313\U00011328"
    u"\U0001132a\U00011330\U00011332\U00011333\U00011335\U00011339\U0001135d\U00011361\U00011400\U00011434\U00011447"
    u"\U0001144a\U0001145f\U00011461\U00011480\U000114af\U000114c4\U000114c5\U00011580\U000115ae\U000115d8\U000115db"
    u"\U00011600\U0001162f\U00011680\U000116aa\U00011700\U0001171a\U00011740\U00011746\U00011800\U0001182b\U000118a0"
    u"\U000118df\U000118ff\U00011906\U0001190c\U00011913\U00011915\U00011916\U00011918\U0001192f\U000119a0\U000119a7"
    u"\U000119aa\U000119d0\U00011a0b\U00011a32\U00011a5c\U00011a89\U00011ab0\U00011af8\U00011c00\U00011c08\U00011c0a"
    u"\U00011c2e\U00011c72\U00011c8f\U00011d00\U00011d06\U00011d08\U00011d09\U00011d0b\U00011d30\U00011d60\U00011d65"
    u"\U00011d67\U00011d68\U00011d6a\U00011d89\U00011ee0\U00011ef2\U00011f04\U00011f10\U00011f12\U00011f33\U00012000"
    u"\U00012399\U00012400\U0001246e\U00012480\U00012543\U00012f90\U00012ff0\U00013000\U0001342f\U00013441\U00013446"
    u"\U00014400\U00014646\U00016800\U00016a38\U00016a40\U00016a5e\U00016a70\U00016abe\U00016ad0\U00016aed\U00016b00"
    u"\U00016b2f\U00016b40\U00016b43\U00016b63\U00016b77\U00016b7d\U00016b8f\U00016e40\U00016e7f\U00016f00\U00016f4a"
    u"\U00016f93\U00016f9f\U00016fe0\U00016fe1\U00017000\U000187f7\U00018800\U00018cd5\U00018d00\U00018d08\U0001aff0"
    u"\U0001aff3\U0001aff5\U0001affb\U0001affd\U0001affe\U0001b000\U0001b122\U0001b150\U0001b152\U0001b164\U0001b167"
    u"\U0001b170\U0001b2fb\U0001bc00\U0001bc6a\U0001bc70\U0001bc7c\U0001bc80\U0001bc88\U0001bc90\U0001bc99\U0001d400"
    u"\U0001d454\U0001d456\U0001d49c\U0001d49e\U0001d49f\U0001d4a5\U0001d4a6\U0001d4a9\U0001d4ac\U0001d4ae\U0001d4b9"
    u"\U0001d4bd\U0001d4c3\U0001d4c5\U0001d505\U0001d507\U0001d50a\U0001d50d\U0001d514\U0001d516\U0001d51c\U0001d51e"
    u"\U0001d539\U0001d53b\U0001d53e\U0001d540\U0001d544\U0001d54a\U0001d550\U0001d552\U0001d6a5\U0001d6a8\U0001d6c0"
    u"\U0001d6c2\U0001d6da\U0001d6dc\U0001d6fa\U0001d6fc\U0001d714\U0001d716\U0001d734\U0001d736\U0001d74e\U0001d750"
    u"\U0001d76e\U0001d770\U0001d788\U0001d78a\U0001d7a8\U0001d7aa\U0001d7c2\U0001d7c4\U0001d7cb\U0001df00\U0001df1e"
    u"\U0001df25\U0001df2a\U0001e030\U0001e06d\U0001e100\U0001e12c\U0001e137\U0001e13d\U0001e290\U0001e2ad\U0001e2c0"
    u"\U0001e2eb\U0001e4d0\U0001e4eb\U0001e7e0\U0001e7e6\U0001e7e8\U0001e7eb\U0001e7ed\U0001e7ee\U0001e7f0\U0001e7fe"
    u"\U0001e800\U0001e8c4\U0001e900\U0001e943\U0001ee00\U0001ee03\U0001ee05\U0001ee1f\U0001ee21\U0001ee22\U0001ee29"
    u"\U0001ee32\U0001ee34\U0001ee37\U0001ee4d\U0001ee4f\U0001ee51\U0001ee52\U0001ee61\U0001ee62\U0001ee67\U0001ee6a"
    u"\U0001ee6c\U0001ee72\U0001ee74\U0001ee77\U0001ee79\U0001ee7c\U0001ee80\U0001ee89\U0001ee8b\U0001ee9b\U0001eea1"
    u"\U0001eea3\U0001eea5\U0001eea9\U0001eeab\U0001eebb\U00020000\U0002a6df\U0002a700\U0002b739\U0002b740\U0002b81d"
    u"\U0002b820\U0002cea1\U0002ceb0\U0002ebe0\U0002f800\U0002fa1d\U00030000\U0003134a"
)
unicode_continuation_ch_any = (
    u"\u00b7\u0387\u05bf\u05c7\u0670\u0711\u07fd\u09bc\u09d7\u09fe\u0a3c"
    u"\u0a51\u0a75\u0abc\u0b3c\u0b82\u0bd7\u0c3c\u0cbc\u0cf3\u0d57\u0dca"
    u"\u0dd6\u0e31\u0eb1\u0f35\u0f37\u0f39\u0fc6\u17dd\u18a9\u1ced\u1cf4"
    u"\u2054\u20e1\u2d7f\ua66f\ua802\ua806\ua80b\ua82c\ua9e5\uaa43\uaab0"
    u"\uaac1\ufb1e\uff3f\U000101fd\U000102e0\U00010a3f\U000110c2\U00011173\U0001123e\U00011241\U00011357"
    u"\U0001145e\U00011940\U000119e4\U00011a47\U00011d3a\U00011d47\U00011f03\U00013440\U00016f4f\U00016fe4\U0001da75"
    u"\U0001da84\U0001e08f\U0001e2ae"
)
unicode_continuation_ch_range = (
    u"\u0030\u0039\u0300\u036f\u0483\u0487\u0591\u05bd\u05c1\u05c2\u05c4"
    u"\u05c5\u0610\u061a\u064b\u0669\u06d6\u06dc\u06df\u06e4\u06e7\u06e8"
    u"\u06ea\u06ed\u06f0\u06f9\u0730\u074a\u07a6\u07b0\u07c0\u07c9\u07eb"
    u"\u07f3\u0816\u0819\u081b\u0823\u0825\u0827\u0829\u082d\u0859\u085b"
    u"\u0898\u089f\u08ca\u08e1\u08e3\u0903\u093a\u093c\u093e\u094f\u0951"
    u"\u0957\u0962\u0963\u0966\u096f\u0981\u0983\u09be\u09c4\u09c7\u09c8"
    u"\u09cb\u09cd\u09e2\u09e3\u09e6\u09ef\u0a01\u0a03\u0a3e\u0a42\u0a47"
    u"\u0a48\u0a4b\u0a4d\u0a66\u0a71\u0a81\u0a83\u0abe\u0ac5\u0ac7\u0ac9"
    u"\u0acb\u0acd\u0ae2\u0ae3\u0ae6\u0aef\u0afa\u0aff\u0b01\u0b03\u0b3e"
    u"\u0b44\u0b47\u0b48\u0b4b\u0b4d\u0b55\u0b57\u0b62\u0b63\u0b66\u0b6f"
    u"\u0bbe\u0bc2\u0bc6\u0bc8\u0bca\u0bcd\u0be6\u0bef\u0c00\u0c04\u0c3e"
    u"\u0c44\u0c46\u0c48\u0c4a\u0c4d\u0c55\u0c56\u0c62\u0c63\u0c66\u0c6f"
    u"\u0c81\u0c83\u0cbe\u0cc4\u0cc6\u0cc8\u0cca\u0ccd\u0cd5\u0cd6\u0ce2"
    u"\u0ce3\u0ce6\u0cef\u0d00\u0d03\u0d3b\u0d3c\u0d3e\u0d44\u0d46\u0d48"
    u"\u0d4a\u0d4d\u0d62\u0d63\u0d66\u0d6f\u0d81\u0d83\u0dcf\u0dd4\u0dd8"
    u"\u0ddf\u0de6\u0def\u0df2\u0df3\u0e33\u0e3a\u0e47\u0e4e\u0e50\u0e59"
    u"\u0eb3\u0ebc\u0ec8\u0ece\u0ed0\u0ed9\u0f18\u0f19\u0f20\u0f29\u0f3e"
    u"\u0f3f\u0f71\u0f84\u0f86\u0f87\u0f8d\u0f97\u0f99\u0fbc\u102b\u103e"
    u"\u1040\u1049\u1056\u1059\u105e\u1060\u1062\u1064\u1067\u106d\u1071"
    u"\u1074\u1082\u108d\u108f\u109d\u135d\u135f\u1369\u1371\u1712\u1715"
    u"\u1732\u1734\u1752\u1753\u1772\u1773\u17b4\u17d3\u17e0\u17e9\u180b"
    u"\u180d\u180f\u1819\u1920\u192b\u1930\u193b\u1946\u194f\u19d0\u19da"
    u"\u1a17\u1a1b\u1a55\u1a5e\u1a60\u1a7c\u1a7f\u1a89\u1a90\u1a99\u1ab0"
    u"\u1abd\u1abf\u1ace\u1b00\u1b04\u1b34\u1b44\u1b50\u1b59\u1b6b\u1b73"
    u"\u1b80\u1b82\u1ba1\u1bad\u1bb0\u1bb9\u1be6\u1bf3\u1c24\u1c37\u1c40"
    u"\u1c49\u1c50\u1c59\u1cd0\u1cd2\u1cd4\u1ce8\u1cf7\u1cf9\u1dc0\u1dff"
    u"\u203f\u2040\u20d0\u20dc\u20e5\u20f0\u2cef\u2cf1\u2de0\u2dff\u302a"
    u"\u302f\u3099\u309a\ua620\ua629\ua674\ua67d\ua69e\ua69f\ua6f0\ua6f1"
    u"\ua823\ua827\ua880\ua881\ua8b4\ua8c5\ua8d0\ua8d9\ua8e0\ua8f1\ua8ff"
    u"\ua909\ua926\ua92d\ua947\ua953\ua980\ua983\ua9b3\ua9c0\ua9d0\ua9d9"
    u"\ua9f0\ua9f9\uaa29\uaa36\uaa4c\uaa4d\uaa50\uaa59\uaa7b\uaa7d\uaab2"
    u"\uaab4\uaab7\uaab8\uaabe\uaabf\uaaeb\uaaef\uaaf5\uaaf6\uabe3\uabea"
    u"\uabec\uabed\uabf0\uabf9\ufe00\ufe0f\ufe20\ufe2f\ufe33\ufe34\ufe4d"
    u"\ufe4f\uff10\uff19\uff9e\uff9f\U00010376\U0001037a\U000104a0\U000104a9\U00010a01\U00010a03"
    u"\U00010a05\U00010a06\U00010a0c\U00010a0f\U00010a38\U00010a3a\U00010ae5\U00010ae6\U00010d24\U00010d27\U00010d30"
    u"\U00010d39\U00010eab\U00010eac\U00010efd\U00010eff\U00010f46\U00010f50\U00010f82\U00010f85\U00011000\U00011002"
    u"\U00011038\U00011046\U00011066\U00011070\U00011073\U00011074\U0001107f\U00011082\U000110b0\U000110ba\U000110f0"
    u"\U000110f9\U00011100\U00011102\U00011127\U00011134\U00011136\U0001113f\U00011145\U00011146\U00011180\U00011182"
    u"\U000111b3\U000111c0\U000111c9\U000111cc\U000111ce\U000111d9\U0001122c\U00011237\U000112df\U000112ea\U000112f0"
    u"\U000112f9\U00011300\U00011303\U0001133b\U0001133c\U0001133e\U00011344\U00011347\U00011348\U0001134b\U0001134d"
    u"\U00011362\U00011363\U00011366\U0001136c\U00011370\U00011374\U00011435\U00011446\U00011450\U00011459\U000114b0"
    u"\U000114c3\U000114d0\U000114d9\U000115af\U000115b5\U000115b8\U000115c0\U000115dc\U000115dd\U00011630\U00011640"
    u"\U00011650\U00011659\U000116ab\U000116b7\U000116c0\U000116c9\U0001171d\U0001172b\U00011730\U00011739\U0001182c"
    u"\U0001183a\U000118e0\U000118e9\U00011930\U00011935\U00011937\U00011938\U0001193b\U0001193e\U00011942\U00011943"
    u"\U00011950\U00011959\U000119d1\U000119d7\U000119da\U000119e0\U00011a01\U00011a0a\U00011a33\U00011a39\U00011a3b"
    u"\U00011a3e\U00011a51\U00011a5b\U00011a8a\U00011a99\U00011c2f\U00011c36\U00011c38\U00011c3f\U00011c50\U00011c59"
    u"\U00011c92\U00011ca7\U00011ca9\U00011cb6\U00011d31\U00011d36\U00011d3c\U00011d3d\U00011d3f\U00011d45\U00011d50"
    u"\U00011d59\U00011d8a\U00011d8e\U00011d90\U00011d91\U00011d93\U00011d97\U00011da0\U00011da9\U00011ef3\U00011ef6"
    u"\U00011f00\U00011f01\U00011f34\U00011f3a\U00011f3e\U00011f42\U00011f50\U00011f59\U00013447\U00013455\U00016a60"
    u"\U00016a69\U00016ac0\U00016ac9\U00016af0\U00016af4\U00016b30\U00016b36\U00016b50\U00016b59\U00016f51\U00016f87"
    u"\U00016f8f\U00016f92\U00016ff0\U00016ff1\U0001bc9d\U0001bc9e\U0001cf00\U0001cf2d\U0001cf30\U0001cf46\U0001d165"
    u"\U0001d169\U0001d16d\U0001d172\U0001d17b\U0001d182\U0001d185\U0001d18b\U0001d1aa\U0001d1ad\U0001d242\U0001d244"
    u"\U0001d7ce\U0001d7ff\U0001da00\U0001da36\U0001da3b\U0001da6c\U0001da9b\U0001da9f\U0001daa1\U0001daaf\U0001e000"
    u"\U0001e006\U0001e008\U0001e018\U0001e01b\U0001e021\U0001e023\U0001e024\U0001e026\U0001e02a\U0001e130\U0001e136"
    u"\U0001e140\U0001e149\U0001e2ec\U0001e2f9\U0001e4ec\U0001e4f9\U0001e8d0\U0001e8d6\U0001e944\U0001e94a\U0001e950"
    u"\U0001e959\U0001fbf0\U0001fbf9"
)

# END GENERATED CODE

#!/usr/local/bin/python

"""
Based on: http://wxpsvg.googlecode.com/svn/trunk/svg/pathdata.py
According to that project, this file is licensed under the LGPL
"""

try:
    from pyparsing import (ParserElement, Literal, Word, CaselessLiteral, 
        Optional, Combine, Forward, ZeroOrMore, nums, oneOf, Group, ParseException, OneOrMore)
except ImportError:
    import sys
    sys.exit("pyparsing is required")
    
    
#ParserElement.enablePackrat()

def Command(char):
    """ Case insensitive but case preserving"""
    return CaselessPreservingLiteral(char)
    
def Arguments(token):
    return Group(token)
    
    
class CaselessPreservingLiteral(CaselessLiteral):
    """ Like CaselessLiteral, but returns the match as found
        instead of as defined.
    """
    def __init__( self, matchString ):
        super().__init__(matchString.upper())
        self.name = "'%s'" % matchString
        self.errmsg = "Expected " + self.name
        self.myException.msg = self.errmsg

    def parseImpl( self, instring, loc, doActions=True ):
        test = instring[ loc:loc+self.matchLen ]
        if test.upper() == self.match:
            return loc+self.matchLen, test
        #~ raise ParseException( instring, loc, self.errmsg )
        exc = self.myException
        exc.loc = loc
        exc.pstr = instring
        raise exc   
    
def Sequence(token):
    """ A sequence of the token"""
    return OneOrMore(token+maybeComma)

digit_sequence = Word(nums)

sign = oneOf("+ -")

def convertToFloat(s, loc, toks):
    try:
        return float(toks[0])
    except:
        raise ParseException(loc, "invalid float format %s"%toks[0])

exponent = CaselessLiteral("e")+Optional(sign)+Word(nums)

#note that almost all these fields are optional, 
#and this can match almost anything. We rely on Pythons built-in
#float() function to clear out invalid values - loosely matching like this
#speeds up parsing quite a lot
floatingPointConstant = Combine(
    Optional(sign) + 
    Optional(Word(nums)) + 
    Optional(Literal(".") + Optional(Word(nums)))+
    Optional(exponent)
)

floatingPointConstant.setParseAction(convertToFloat)

number = floatingPointConstant

#same as FP constant but don't allow a - sign
nonnegativeNumber = Combine(
    Optional(Word(nums)) + 
    Optional(Literal(".") + Optional(Word(nums)))+
    Optional(exponent)
)
nonnegativeNumber.setParseAction(convertToFloat)

coordinate = number

#comma or whitespace can separate values all over the place in SVG
maybeComma = Optional(Literal(',')).suppress()

coordinateSequence = Sequence(coordinate)

coordinatePair = (coordinate + maybeComma + coordinate).setParseAction(lambda t: tuple(t))
coordinatePairSequence = Sequence(coordinatePair)

coordinatePairPair = coordinatePair + maybeComma + coordinatePair
coordinatePairPairSequence = Sequence(Group(coordinatePairPair))

coordinatePairTriple = coordinatePair + maybeComma + coordinatePair + maybeComma + coordinatePair
coordinatePairTripleSequence = Sequence(Group(coordinatePairTriple))

#commands
lineTo = Group(Command("L") + Arguments(coordinatePairSequence))
curve = Group(Command("C") + Arguments(coordinatePairSequence))

moveTo = Group(Command("M") + Arguments(coordinatePairSequence))

closePath = Group(Command("Z")).setParseAction(lambda t: ('Z', (None,)))

flag = oneOf("1 0").setParseAction(lambda t: bool(int((t[0]))))

arcRadius = (
    nonnegativeNumber + maybeComma + #rx
    nonnegativeNumber #ry
).setParseAction(lambda t: tuple(t))

arcFlags = (flag + maybeComma + flag).setParseAction(lambda t: tuple(t))

ellipticalArcArgument = Group(
    arcRadius + maybeComma + #rx, ry
    number + maybeComma +#rotation
    arcFlags + #large-arc-flag, sweep-flag
    coordinatePair #(x,y)
)

ellipticalArc = Group(Command("A") + Arguments(Sequence(ellipticalArcArgument)))

smoothQuadraticBezierCurveto = Group(Command("T") + Arguments(coordinatePairSequence))

quadraticBezierCurveto = Group(Command("Q") + Arguments(coordinatePairPairSequence))

smoothCurve = Group(Command("S") + Arguments(coordinatePairPairSequence))

#curve = Group(Command("C") + Arguments(coordinatePairTripleSequence))

horizontalLine = Group(Command("H") + Arguments(coordinateSequence))
verticalLine = Group(Command("V") + Arguments(coordinateSequence))

drawToCommand = (
    lineTo | moveTo | closePath | ellipticalArc | smoothQuadraticBezierCurveto |
    quadraticBezierCurveto | smoothCurve | curve | horizontalLine | verticalLine
    )

#~ number.debug = True
moveToDrawToCommands = moveTo + ZeroOrMore(drawToCommand)

path = ZeroOrMore(moveToDrawToCommands)
path.keepTabs = True

def get_points(d):
    commands = path.parseString(d)
    points = []
    currentset = None
    for command in commands:
        if command[0] == 'M' or command[0] == 'm':
            currentset = []
            points.append(currentset)
            currentset.append(command[1][-1])
        elif command[0] == 'L' or command[0] == 'l':
            currentset.extend(command[1])
        elif command[0] == 'C' or command[0] == 'c':
            currentset.extend(command[1])
    return points

if __name__ == "__main__":
    s = ("M 242.96145,653.59282 L 244.83646,650.1553 L 247.02397,649.8428 "
         "L 247.33647,650.62405 L 245.30521,653.59282 L 242.96145,653.59282 z "
         "M 252.80525,649.99905 L 258.74278,652.49906 L 260.77404,652.18656 "
         "L 262.33654,648.43654 L 261.71154,645.15528 L 257.64902,644.68653 "
         "L 253.74275,646.40528 L 252.80525,649.99905 z M 282.49289,659.6866 "
         "L 286.08665,664.99912 L 288.43041,664.68662 L 289.52417,664.21787 "
         "L 290.93042,665.46787 L 294.52419,665.31162 L 295.4617,663.90537 "
         "L 292.64918,662.18661 L 290.77417,658.59284 L 288.74291,655.15533 "
         "L 283.11789,657.96784 L 282.49289,659.6866 z M 302.02423,668.28039 "
         "L 303.27423,666.40538 L 307.8055,667.34288 L 308.43051,666.87413 "
         "L 314.36803,667.49913 L 314.05553,668.74914 L 311.55552,670.15539 "
         "L 307.33675,669.84289 L 302.02423,668.28039 z M 307.1805,673.28041 "
         "L 309.05551,677.03043 L 312.02427,675.93667 L 312.33677,674.37416 "
         "L 310.77427,672.3429 L 307.1805,672.0304 L 307.1805,673.28041 z "
         "M 313.89928,672.18665 L 316.08679,669.37414 L 320.61806,671.7179 "
         "L 324.83683,672.81166 L 329.0556,675.46792 L 329.0556,677.34293 "
         "L 325.61809,679.06169 L 320.93056,679.99919 L 318.5868,678.59293 "
         "L 313.89928,672.18665 z M 329.99311,687.18672 L 331.55561,685.93672 "
         "L 334.83688,687.49923 L 342.18066,690.93674 L 345.46193,692.968 "
         "L 347.02443,695.31176 L 348.89944,699.53053 L 352.80571,702.03054 "
         "L 352.49321,703.28055 L 348.74319,706.40556 L 344.68067,707.81182 "
         "L 343.27442,707.18682 L 340.30565,708.90557 L 337.96189,712.03059 "
         "L 335.77438,714.8431 L 334.05562,714.68685 L 330.61811,712.18684 "
         "L 330.30561,707.81182 L 330.93061,705.46806 L 329.3681,699.99928 "
         "L 327.33684,698.28052 L 327.18059,695.78051 L 329.3681,694.84301 "
         "L 331.39936,691.87425 L 331.86811,690.93674 L 330.30561,689.21798 "
         "L 329.99311,687.18672 z ")
    print(path.parseString(s))

from __future__ import unicode_literals

import unittest
from io import StringIO
import string

from .. import Scanning
from ..Symtab import ModuleScope
from ..TreeFragment import StringParseContext
from ..Errors import init_thread

# generate some fake code - just a bunch of lines of the form "a0 a1 ..."
code = []
for ch in string.ascii_lowercase:
    line = " ".join(["%s%s" % (ch, n) for n in range(10)])
    code.append(line)
code = "\n".join(code)

init_thread()


class TestScanning(unittest.TestCase):
    def make_scanner(self):
        source = Scanning.StringSourceDescriptor("fake code", code)
        buf = StringIO(code)
        context = StringParseContext("fake context")
        scope = ModuleScope("fake_module", None, None)

        return Scanning.PyrexScanner(buf, source, scope=scope, context=context)

    def test_put_back_positions(self):
        scanner = self.make_scanner()

        self.assertEqual(scanner.sy, "IDENT")
        self.assertEqual(scanner.systring, "a0")
        scanner.next()
        self.assertEqual(scanner.sy, "IDENT")
        self.assertEqual(scanner.systring, "a1")
        a1pos = scanner.position()
        self.assertEqual(a1pos[1:], (1, 3))
        a2peek = scanner.peek()  # shouldn't mess up the position
        self.assertEqual(a1pos, scanner.position())
        scanner.next()
        self.assertEqual(a2peek, (scanner.sy, scanner.systring))

        # find next line
        while scanner.sy != "NEWLINE":
            scanner.next()

        line_sy = []
        line_systring = []
        line_pos = []

        scanner.next()
        while scanner.sy != "NEWLINE":
            line_sy.append(scanner.sy)
            line_systring.append(scanner.systring)
            line_pos.append(scanner.position())
            scanner.next()

        for sy, systring, pos in zip(
            line_sy[::-1], line_systring[::-1], line_pos[::-1]
        ):
            scanner.put_back(sy, systring, pos)

        n = 0
        while scanner.sy != "NEWLINE":
            self.assertEqual(scanner.sy, line_sy[n])
            self.assertEqual(scanner.systring, line_systring[n])
            self.assertEqual(scanner.position(), line_pos[n])
            scanner.next()
            n += 1

        self.assertEqual(n, len(line_pos))

    def test_tentatively_scan(self):
        scanner = self.make_scanner()
        with Scanning.tentatively_scan(scanner) as errors:
            while scanner.sy != "NEWLINE":
                scanner.next()
        self.assertFalse(errors)

        scanner.next()
        self.assertEqual(scanner.systring, "b0")
        pos = scanner.position()
        with Scanning.tentatively_scan(scanner) as errors:
            while scanner.sy != "NEWLINE":
                scanner.next()
                if scanner.systring == "b7":
                    scanner.error("Oh no not b7!")
                    break
        self.assertTrue(errors)
        self.assertEqual(scanner.systring, "b0")  # state has been restored
        self.assertEqual(scanner.position(), pos)
        scanner.next()
        self.assertEqual(scanner.systring, "b1")  # and we can keep going again
        scanner.next()
        self.assertEqual(scanner.systring, "b2")  # and we can keep going again

        with Scanning.tentatively_scan(scanner) as error:
            scanner.error("Something has gone wrong with the current symbol")
        self.assertEqual(scanner.systring, "b2")
        scanner.next()
        self.assertEqual(scanner.systring, "b3")

        # test a few combinations of nested scanning
        sy1, systring1 = scanner.sy, scanner.systring
        pos1 = scanner.position()
        with Scanning.tentatively_scan(scanner):
            scanner.next()
            sy2, systring2 = scanner.sy, scanner.systring
            pos2 = scanner.position()
            with Scanning.tentatively_scan(scanner):
                with Scanning.tentatively_scan(scanner):
                    scanner.next()
                    scanner.next()
                    scanner.error("Ooops")
                self.assertEqual((scanner.sy, scanner.systring), (sy2, systring2))
            self.assertEqual((scanner.sy, scanner.systring), (sy2, systring2))
            scanner.error("eee")
        self.assertEqual((scanner.sy, scanner.systring), (sy1, systring1))
        with Scanning.tentatively_scan(scanner):
            scanner.next()
            scanner.next()
            with Scanning.tentatively_scan(scanner):
                scanner.next()
                # no error - but this block should be unwound by the outer block too
            scanner.next()
            scanner.error("Oooops")
        self.assertEqual((scanner.sy, scanner.systring), (sy1, systring1))




if __name__ == "__main__":
    unittest.main()

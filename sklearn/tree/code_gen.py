import contextlib
import tempfile
import os
import subprocess

from distutils import sysconfig

CXX_COMPILER = sysconfig.get_config_var('CXX')

EVALUATE_FN_NAME = "evaluate"
ALWAYS_INLINE = "__attribute__((__always_inline__))"


class CodeGenerator(object):
    def __init__(self):
        self._lines = []
        self._indent = 0

    @property
    def lines(self):
        return self._lines

    def write(self, line):
        self._lines.append("  " * self._indent + line)

    @contextlib.contextmanager
    def bracketed(self, preamble, postamble):
        assert self._indent >= 0
        self.write(preamble)
        self._indent += 1
        yield
        self._indent -= 1
        self.write(postamble)


def code_gen_tree(tree, evaluate_fn=EVALUATE_FN_NAME, gen=None):
    """
    Generates C code representing the evaluation of a tree.

    Writes code similar to:
    ```
        extern "C" {
          __attribute__((__always_inline__)) float evaluate(float* f) {
            if (f[9] <= 0.175931170583) {
              return 0.0;
            }
            else {
              return 1.0;
            }
          }
        }
    ```

    to the given CodeGenerator object.
    """
    if gen is None:
        gen = CodeGenerator()

    def recur(node):
        if tree.children_left[node] == -1:
            assert tree.value[node].size == 1
            gen.write("return {0};".format(tree.value[node][0][0]))
            return

        branch = "if (f[{feature}] <= {threshold}) {{".format(
            feature=tree.feature[node],
            threshold=tree.threshold[node])
        with gen.bracketed(branch, "}"):
            recur(tree.children_left[node])

        with gen.bracketed("else {", "}"):
            recur(tree.children_right[node])

    with gen.bracketed('extern "C" {', "}"):
        fn_decl = "{inline} float {name}(float* f) {{".format(
            inline=ALWAYS_INLINE,
            name=evaluate_fn)
        with gen.bracketed(fn_decl, "}"):
            recur(0)
    return gen.lines


def code_gen_ensemble(trees, individual_learner_weight, initial_value,
                      gen=None):
    """
    Writes code similar to:
    ```
        extern "C" {
          __attribute__((__always_inline__)) float evaluate_partial_0(float* f) {
            if (f[4] <= 0.662200987339) {
              return 1.0;
            }
            else {
              if (f[8] <= 0.804652512074) {
                return 0.0;
              }
              else {
                return 1.0;
              }
            }
          }
        }
        extern "C" {
          __attribute__((__always_inline__)) float evaluate_partial_1(float* f) {
            if (f[4] <= 0.694428026676) {
              return 1.0;
            }
            else {
              if (f[7] <= 0.4402526021) {
                return 1.0;
              }
              else {
                return 0.0;
              }
            }
          }
        }

        extern "C" {
          float evaluate(float* f) {
            float result = 0.0;
            result += evaluate_partial_0(f) * 0.1;
            result += evaluate_partial_1(f) * 0.1;
            return result;
          }
        }
    ```

    to the given CodeGenerator object.
    """

    if gen is None:
        gen = CodeGenerator()

    for i, tree in enumerate(trees):
        name = "{name}_{index}".format(name=EVALUATE_FN_NAME, index=i)
        code_gen_tree(tree, name, gen)

    with gen.bracketed('extern "C" {', "}"):
        fn_decl = "float {name}(float* f) {{".format(name=EVALUATE_FN_NAME)
        with gen.bracketed(fn_decl, "}"):
            gen.write("float result = {0};".format(initial_value))
            for i, _ in enumerate(trees):
                increment = "result += {name}_{index}(f) * {weight};".format(
                    name=EVALUATE_FN_NAME,
                    index=i,
                    weight=individual_learner_weight)
                gen.write(increment)
            gen.write("return result;")
    return gen.lines


def compile_code_to_object(code):
    # XXX - should we clean up the temporary directory left from mkdtemp()?
    tmpdir = tempfile.mkdtemp()
    cpp_f = os.path.join(tmpdir, "tree.cpp")
    so_f = os.path.join(tmpdir, "tree.so")
    o_f = os.path.join(tmpdir, "tree.o")

    with open(cpp_f, 'w') as f:
        f.write(code)

    subprocess.check_call([CXX_COMPILER, cpp_f,  "-c", "-o", o_f, "-O3"])
    subprocess.check_call([CXX_COMPILER, "-shared", o_f, "-dynamiclib",
                           "-fpic", "-flto", "-o", so_f, "-O3"])
    return so_f

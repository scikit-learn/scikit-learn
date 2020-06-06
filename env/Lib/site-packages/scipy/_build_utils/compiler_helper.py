"""
Helpers for detection of compiler features
"""
import tempfile
import os

def try_compile(compiler, code=None, flags=[], ext=None):
    """Returns True if the compiler is able to compile the given code"""
    from distutils.errors import CompileError

    code = code or 'int main (int argc, char **argv) { return 0; }'
    ext = ext or compiler.src_extensions[0]

    with tempfile.TemporaryDirectory() as temp_dir:
        fname = os.path.join(temp_dir, 'main'+ext)
        with open(fname, 'w') as f:
            f.write(code)

        try:
            compiler.compile([fname], output_dir=temp_dir, extra_postargs=flags)
        except CompileError:
            return False
    return True


def has_flag(compiler, flag, ext=None):
    """Returns True if the compiler supports the given flag"""
    return try_compile(compiler, flags=[flag], ext=ext)


def get_cxx_std_flag(compiler):
    """Detects compiler flag for c++14, c++11 or None if not detected"""
    # GNU C compiler documentation use single dash:
    #    https://gcc.gnu.org/onlinedocs/gcc/Standards.html
    # but silently understands two dahes like --std=c++11 too.
    # Other GCC compatible compilers, like Intel C Compiler on Linux do not.
    gnu_flags = ['-std=c++14', '-std=c++11']
    flags_by_cc = {
        'msvc': ['/std:c++14', None],
        'intelw': ['/Qstd=c++14', '/Qstd=c++11'],
        'intelem': ['-std=c++14', '-std=c++11']
    }
    flags = flags_by_cc.get(compiler.compiler_type, gnu_flags)

    for flag in flags:
        if flag is None:
            return None

        if has_flag(compiler, flag):
            return flag

    from numpy.distutils import log
    log.warn('Could not detect c++ standard flag')
    return None


def try_add_flag(args, compiler, flag, ext=None):
    """Appends flag to the list of arguments if supported by the compiler"""
    if try_compile(compiler, flags=args+[flag], ext=ext):
        args.append(flag)

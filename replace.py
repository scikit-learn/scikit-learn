import re
from pathlib import Path
from itertools import chain

sklearn_root = Path("./sklearn")

pyx_paths = chain(
    sklearn_root.glob("**/*.pyx"),
    sklearn_root.glob("**/*.pyx.tp"),
    sklearn_root.glob("**/*.pxi"),
)


nogil_colon = re.compile(r"\) nogil:")
nogil_except_neg1 = re.compile(r"\) nogil except -1:")
nogil_except_star = re.compile(r"\) nogil except [*]:")

for pyx_path in pyx_paths:
    orig_contents = pyx_path.read_text()
    new_contents = nogil_colon.sub(") noexcept nogil:", orig_contents)
    new_contents = nogil_except_neg1.sub(") except -1 nogil:", new_contents)
    new_contents = nogil_except_star.sub(") except * nogil:", new_contents)

    if new_contents != orig_contents:
        pyx_path.write_text(new_contents)


nogil_no_colon = re.compile(r"\) nogil\n")
nogil_except_neg1_no_colon = re.compile(r"\) nogil except -1")
nogil_except_star_no_colon = re.compile(r"\) nogil except \*")
pxd_paths = chain(sklearn_root.glob("**/*.pxd"), sklearn_root.glob("**/*.pxd.tp"))

for pxd_path in pxd_paths:
    orig_contents = pxd_path.read_text()
    new_contents = nogil_no_colon.sub(") noexcept nogil\n", orig_contents)
    new_contents = nogil_except_neg1_no_colon.sub(") except -1 nogil", new_contents)
    new_contents = nogil_except_star_no_colon.sub(") except * nogil", new_contents)

    if new_contents != orig_contents:
        pxd_path.write_text(new_contents)
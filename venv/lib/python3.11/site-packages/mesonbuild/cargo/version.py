# SPDX-License-Identifier: Apache-2.0
# Copyright © 2022-2023 Intel Corporation

"""Convert Cargo versions into Meson compatible ones."""

from __future__ import annotations
import typing as T


def convert(cargo_ver: str) -> T.List[str]:
    """Convert a Cargo compatible version into a Meson compatible one.

    :param cargo_ver: The version, as Cargo specifies
    :return: A list of version constraints, as Meson understands them
    """
    # Cleanup, just for safety
    cargo_ver = cargo_ver.strip()
    cargo_vers = [c.strip() for c in cargo_ver.split(',')]

    out: T.List[str] = []

    for ver in cargo_vers:
        # This covers >= and =< as well
        # https://doc.rust-lang.org/cargo/reference/specifying-dependencies.html#comparison-requirements
        if ver.startswith(('>', '<', '=')):
            out.append(ver)

        elif ver.startswith('~'):
            # Rust has these tilde requirements, which means that it is >= to
            # the version, but less than the next version
            # https://doc.rust-lang.org/cargo/reference/specifying-dependencies.html#tilde-requirements
            # we convert those into a pair of constraints
            v = ver[1:].split('.')
            out.append(f'>= {".".join(v)}')
            if len(v) == 3:
                out.append(f'< {v[0]}.{int(v[1]) + 1}.0')
            elif len(v) == 2:
                out.append(f'< {v[0]}.{int(v[1]) + 1}')
            else:
                out.append(f'< {int(v[0]) + 1}')

        elif '*' in ver:
            # Rust has astrisk requirements,, which are like 1.* == ~1
            # https://doc.rust-lang.org/cargo/reference/specifying-dependencies.html#wildcard-requirements
            v = ver.split('.')[:-1]
            if v:
                out.append(f'>= {".".join(v)}')
            if len(v) == 2:
                out.append(f'< {v[0]}.{int(v[1]) + 1}')
            elif len(v) == 1:
                out.append(f'< {int(v[0]) + 1}')

        else:
            # a Caret version is equivalent to the default strategy
            # https://doc.rust-lang.org/cargo/reference/specifying-dependencies.html#caret-requirements
            if ver.startswith('^'):
                ver = ver[1:]

            # If there is no qualifier, then it means this or the next non-zero version
            # That means that if this is `1.1.0``, then we need `>= 1.1.0` && `< 2.0.0`
            # Or if we have `0.1.0`, then we need `>= 0.1.0` && `< 0.2.0`
            # Or if we have `0.1`, then we need `>= 0.1.0` && `< 0.2.0`
            # Or if we have `0.0.0`, then we need `< 1.0.0`
            # Or if we have `0.0`, then we need `< 1.0.0`
            # Or if we have `0`, then we need `< 1.0.0`
            # Or if we have `0.0.3`, then we need `>= 0.0.3` && `< 0.0.4`
            # https://doc.rust-lang.org/cargo/reference/specifying-dependencies.html#specifying-dependencies-from-cratesio
            #
            # this works much like the ~ versions, but in reverse. Tilde starts
            # at the patch version and works up, to the major version, while
            # bare numbers start at the major version and work down to the patch
            # version
            vers = ver.split('.')
            min_: T.List[str] = []
            max_: T.List[str] = []
            bumped = False
            for v_ in vers:
                if v_ != '0' and not bumped:
                    min_.append(v_)
                    max_.append(str(int(v_) + 1))
                    bumped = True
                else:
                    min_.append(v_)
                    if not bumped:
                        max_.append('0')

            # If there is no minimum, don't emit one
            if set(min_) != {'0'}:
                out.append('>= {}'.format('.'.join(min_)))
            if set(max_) != {'0'}:
                out.append('< {}'.format('.'.join(max_)))
            else:
                out.append('< 1')

    return out

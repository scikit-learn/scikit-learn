"""This module provides contains the code for intersphinx command-line utilities."""

from __future__ import annotations

import sys

from sphinx.ext.intersphinx._load import _fetch_inventory


def inspect_main(argv: list[str], /) -> int:
    """Debug functionality to print out an inventory"""
    if len(argv) < 1:
        print(
            'Print out an inventory file.\n'
            'Error: must specify local path or URL to an inventory file.',
            file=sys.stderr,
        )
        return 1

    class MockConfig:
        intersphinx_timeout: int | None = None
        tls_verify = False
        tls_cacerts: str | dict[str, str] | None = None
        user_agent: str = ''

    try:
        filename = argv[0]
        inv_data = _fetch_inventory(
            target_uri='',
            inv_location=filename,
            config=MockConfig(),  # type: ignore[arg-type]
            srcdir='',  # type: ignore[arg-type]
        )
        for key in sorted(inv_data or {}):
            print(key)
            inv_entries = sorted(inv_data[key].items())
            for entry, (_proj, _ver, url_path, display_name) in inv_entries:
                display_name = display_name * (display_name != '-')
                print(f'    {entry:<40} {display_name:<40}: {url_path}')
    except ValueError as exc:
        print(exc.args[0] % exc.args[1:], file=sys.stderr)
        return 1
    except Exception as exc:
        print(f'Unknown error: {exc!r}', file=sys.stderr)
        return 1
    else:
        return 0

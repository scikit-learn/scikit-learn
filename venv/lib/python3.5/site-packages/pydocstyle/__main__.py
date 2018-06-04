#! /usr/bin/env python
"""Static analysis tool for checking docstring conventions and style.

The repository is located at:
http://github.com/PyCQA/pydocstyle

"""


__all__ = ()


def main():
    from pydocstyle import cli
    cli.main()


if __name__ == '__main__':
    main()

"""Base 'sphinx' command.

Subcommands are loaded lazily from the ``_COMMANDS`` table for performance.

All subcommand modules must define three attributes:

- ``parser_description``, a description of the subcommand. The first paragraph
  is taken as the short description for the command.
- ``set_up_parser``, a callable taking and returning an ``ArgumentParser``. This
  function is responsible for adding options and arguments to the subcommand's
  parser.
- ``run``, a callable taking parsed arguments and returning an exit code. This
  function is responsible for running the main body of the subcommand and
  returning the exit status.

The entire ``sphinx._cli`` namespace is private, only the command line interface
has backwards-compatibility guarantees.
"""

from __future__ import annotations

import argparse
import importlib
import locale
import sys
from typing import TYPE_CHECKING

from sphinx._cli.util.colour import (
    bold,
    disable_colour,
    enable_colour,
    terminal_supports_colour,
    underline,
)
from sphinx.locale import __, init_console

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Iterator, Sequence
    from typing import NoReturn, TypeAlias

    _PARSER_SETUP: TypeAlias = Callable[
        [argparse.ArgumentParser], argparse.ArgumentParser
    ]
    _RUNNER: TypeAlias = Callable[[argparse.Namespace], int]

    from typing import Protocol

    class _SubcommandModule(Protocol):
        parser_description: str
        set_up_parser: _PARSER_SETUP  # takes and returns argument parser
        run: _RUNNER  # takes parsed args, returns exit code


# Map of command name to import path.
_COMMANDS: dict[str, str] = {}


def _load_subcommand_descriptions() -> Iterator[tuple[str, str]]:
    for command, module_name in _COMMANDS.items():
        module: _SubcommandModule = importlib.import_module(module_name)
        try:
            description = module.parser_description
        except AttributeError:
            # log an error here, but don't fail the full enumeration
            print(f'Failed to load the description for {command}', file=sys.stderr)
        else:
            yield command, description.split('\n\n', 1)[0]


class _RootArgumentParser(argparse.ArgumentParser):
    def format_help(self) -> str:
        help_fragments: list[str] = [
            bold(underline(__('Usage:'))),
            ' ',
            __('{0} [OPTIONS] <COMMAND> [<ARGS>]').format(bold(self.prog)),
            '\n',
            '\n',
            __('  The Sphinx documentation generator.'),
            '\n',
        ]

        if commands := list(_load_subcommand_descriptions()):
            command_lengths = map(len, next(zip(*commands, strict=True), ()))
            command_max_length = min(max(command_lengths), 22)
            help_fragments += [
                '\n',
                bold(underline(__('Commands:'))),
                '\n',
            ]
            help_fragments += [
                f'  {command_name: <{command_max_length}}  {command_desc}'
                for command_name, command_desc in commands
            ]
            help_fragments.append('\n')

        # self._action_groups[1] is self._optionals
        # Uppercase the title of the Optionals group
        self._optionals.title = __('Options')
        for argument_group in self._action_groups[1:]:
            if arguments := [
                action
                for action in argument_group._group_actions
                if action.help != argparse.SUPPRESS
            ]:
                help_fragments += self._format_optional_arguments(
                    arguments,
                    argument_group.title or '',
                )

        help_fragments += [
            '\n',
            __(
                'For more information, visit https://www.sphinx-doc.org/en/master/man/.'
            ),
            '\n',
        ]
        return ''.join(help_fragments)

    def _format_optional_arguments(
        self,
        actions: Iterable[argparse.Action],
        title: str,
    ) -> Iterator[str]:
        yield '\n'
        yield bold(underline(title + ':'))
        yield '\n'

        for action in actions:
            prefix = '    ' * all(o[1] == '-' for o in action.option_strings)
            opt = prefix + '  ' + ', '.join(map(bold, action.option_strings))
            if action.nargs != 0:
                opt += ' ' + self._format_metavar(
                    action.nargs, action.metavar, action.choices, action.dest
                )
            yield opt
            yield '\n'
            if action_help := (action.help or '').strip():
                yield from (f'        {line}\n' for line in action_help.splitlines())

    @staticmethod
    def _format_metavar(
        nargs: int | str | None,
        metavar: str | tuple[str, ...] | None,
        choices: Iterable[str] | None,
        dest: str,
    ) -> str:
        if metavar is None:
            if choices is not None:
                metavar = '{' + ', '.join(sorted(choices)) + '}'
            else:
                metavar = dest.upper()
        if nargs is None:
            return f'{metavar}'
        elif nargs == argparse.OPTIONAL:
            return f'[{metavar}]'
        elif nargs == argparse.ZERO_OR_MORE:
            if len(metavar) == 2:
                return f'[{metavar[0]} [{metavar[1]} ...]]'
            else:
                return f'[{metavar} ...]'
        elif nargs == argparse.ONE_OR_MORE:
            return f'{metavar} [{metavar} ...]'
        elif nargs == argparse.REMAINDER:
            return '...'
        elif nargs == argparse.PARSER:
            return f'{metavar} ...'
        msg = 'invalid nargs value'
        raise ValueError(msg)

    def error(self, message: str) -> NoReturn:
        sys.stderr.write(
            __(
                '{0}: error: {1}\n' "Run '{0} --help' for information"  # NoQA: COM812
            ).format(self.prog, message)
        )
        raise SystemExit(2)


def _create_parser() -> _RootArgumentParser:
    parser = _RootArgumentParser(
        prog='sphinx',
        description=__('   Manage documentation with Sphinx.'),
        epilog=__(
            'For more information, visit https://www.sphinx-doc.org/en/master/man/.'
        ),
        add_help=False,
        allow_abbrev=False,
    )
    parser.add_argument(
        '-V',
        '--version',
        action='store_true',
        default=argparse.SUPPRESS,
        help=__('Show the version and exit.'),
    )
    parser.add_argument(
        '-h',
        '-?',
        '--help',
        action='store_true',
        default=argparse.SUPPRESS,
        help=__('Show this message and exit.'),
    )

    # logging control
    log_control = parser.add_argument_group(__('Logging'))
    log_control.add_argument(
        '-v',
        '--verbose',
        action='count',
        dest='verbosity',
        default=0,
        help=__('Increase verbosity (can be repeated)'),
    )
    log_control.add_argument(
        '-q',
        '--quiet',
        action='store_const',
        dest='verbosity',
        const=-1,
        help=__('Only print errors and warnings.'),
    )
    log_control.add_argument(
        '--silent',
        action='store_const',
        dest='verbosity',
        const=-2,
        help=__('No output at all'),
    )

    parser.add_argument(
        'COMMAND',
        nargs=argparse.REMAINDER,
        metavar=__('<command>'),
    )
    return parser


def _parse_command(argv: Sequence[str] = ()) -> tuple[str, Sequence[str]]:
    parser = _create_parser()
    args = parser.parse_args(argv)
    command_name, *command_argv = args.COMMAND or ('help',)
    command_name = command_name.lower()

    if terminal_supports_colour():
        enable_colour()
    else:
        disable_colour()

    # Handle '--version' or '-V' passed to the main command or any subcommand
    if 'version' in args or {'-V', '--version'}.intersection(command_argv):
        from sphinx import __display_version__

        sys.stderr.write(f'sphinx {__display_version__}\n')
        raise SystemExit(0)

    # Handle '--help' or '-h' passed to the main command (subcommands may have
    # their own help text)
    if 'help' in args or command_name == 'help':
        sys.stderr.write(parser.format_help())
        raise SystemExit(0)

    if command_name not in _COMMANDS:
        sys.stderr.write(
            __(
                f'sphinx: {command_name!r} is not a sphinx command. '
                "See 'sphinx --help'.\n"
            )
        )
        raise SystemExit(2)

    return command_name, command_argv


def _load_subcommand(command_name: str) -> tuple[str, _PARSER_SETUP, _RUNNER]:
    try:
        module: _SubcommandModule = importlib.import_module(_COMMANDS[command_name])
    except KeyError:
        msg = f'invalid command name {command_name!r}.'
        raise ValueError(msg) from None
    return module.parser_description, module.set_up_parser, module.run


def _create_sub_parser(
    command_name: str,
    description: str,
    parser_setup: _PARSER_SETUP,
) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=f'sphinx {command_name}',
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        allow_abbrev=False,
    )
    return parser_setup(parser)


def run(argv: Sequence[str] = (), /) -> int:
    locale.setlocale(locale.LC_ALL, '')
    init_console()

    argv = argv or sys.argv[1:]
    try:
        cmd_name, cmd_argv = _parse_command(argv)
        cmd_description, set_up_parser, runner = _load_subcommand(cmd_name)
        cmd_parser = _create_sub_parser(cmd_name, cmd_description, set_up_parser)
        cmd_args = cmd_parser.parse_args(cmd_argv)
        return runner(cmd_args)
    except SystemExit as exc:
        return exc.code  # type: ignore[return-value]
    except (Exception, KeyboardInterrupt):
        return 2


if __name__ == '__main__':
    raise SystemExit(run())

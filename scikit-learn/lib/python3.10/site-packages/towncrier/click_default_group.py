# SPDX-FileCopyrightText: 2015 Heungsub Lee <19982+sublee@users.noreply.github.com>
#
# SPDX-License-Identifier: BSD-3-Clause

# Vendored from
# https://github.com/click-contrib/click-default-group/tree/b671ae5325d186fe5ea7abb584f15852a1e931aa
# Because the PyPI package could not be installed on modern Pips anymore and
# the project looks unmaintaintained.

"""
click_default_group
~~~~~~~~~~~~~~~~~~~

Define a default subcommand by `default=True`:

.. sourcecode:: python

   import click
   from click_default_group import DefaultGroup

   @click.group(cls=DefaultGroup, default_if_no_args=True)
   def cli():
       pass

   @cli.command(default=True)
   def foo():
       click.echo('foo')

   @cli.command()
   def bar():
       click.echo('bar')

Then you can invoke that without explicit subcommand name:

.. sourcecode:: console

   $ cli.py --help
   Usage: cli.py [OPTIONS] COMMAND [ARGS]...

   Options:
     --help    Show this message and exit.

   Command:
     foo*
     bar

   $ cli.py
   foo
   $ cli.py foo
   foo
   $ cli.py bar
   bar

"""
import warnings

import click


__all__ = ["DefaultGroup"]
__version__ = "1.2.2"


class DefaultGroup(click.Group):
    """Invokes a subcommand marked with `default=True` if any subcommand not
    chosen.

    :param default_if_no_args: resolves to the default command if no arguments
                               passed.

    """

    def __init__(self, *args, **kwargs):
        # To resolve as the default command.
        if not kwargs.get("ignore_unknown_options", True):
            raise ValueError("Default group accepts unknown options")
        self.ignore_unknown_options = True
        self.default_cmd_name = kwargs.pop("default", None)
        self.default_if_no_args = kwargs.pop("default_if_no_args", False)
        super().__init__(*args, **kwargs)

    def set_default_command(self, command):
        """Sets a command function as the default command."""
        cmd_name = command.name
        self.add_command(command)
        self.default_cmd_name = cmd_name

    def parse_args(self, ctx, args):
        if not args and self.default_if_no_args:
            args.insert(0, self.default_cmd_name)
        return super().parse_args(ctx, args)

    def get_command(self, ctx, cmd_name):
        if cmd_name not in self.commands:
            # No command name matched.
            ctx.arg0 = cmd_name
            cmd_name = self.default_cmd_name
        return super().get_command(ctx, cmd_name)

    def resolve_command(self, ctx, args):
        base = super()
        cmd_name, cmd, args = base.resolve_command(ctx, args)
        if hasattr(ctx, "arg0"):
            args.insert(0, ctx.arg0)
            cmd_name = cmd.name
        return cmd_name, cmd, args

    def format_commands(self, ctx, formatter):
        formatter = DefaultCommandFormatter(self, formatter, mark="*")
        return super().format_commands(ctx, formatter)

    def command(self, *args, **kwargs):
        default = kwargs.pop("default", False)
        decorator = super().command(*args, **kwargs)
        if not default:
            return decorator
        warnings.warn(
            "Use default param of DefaultGroup or " "set_default_command() instead",
            DeprecationWarning,
        )

        def _decorator(f):
            cmd = decorator(f)
            self.set_default_command(cmd)
            return cmd

        return _decorator


class DefaultCommandFormatter:
    """Wraps a formatter to mark a default command."""

    def __init__(self, group, formatter, mark="*"):
        self.group = group
        self.formatter = formatter
        self.mark = mark

    def __getattr__(self, attr):
        return getattr(self.formatter, attr)

    def write_dl(self, rows, *args, **kwargs):
        rows_ = []
        for cmd_name, help in rows:
            if cmd_name == self.group.default_cmd_name:
                rows_.insert(0, (cmd_name + self.mark, help))
            else:
                rows_.append((cmd_name, help))
        return self.formatter.write_dl(rows_, *args, **kwargs)

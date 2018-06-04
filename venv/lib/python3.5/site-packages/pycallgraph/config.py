import argparse
import sys

from .output import outputters
from .globbing_filter import GlobbingFilter


class Config(object):
    '''Handles configuration settings for pycallgraph, tracer, and each output
    module.  It also handles command line arguments.
    '''

    def __init__(self, **kwargs):
        '''
        You can set defaults in the constructor, e.g. Config(verbose=True)
        '''
        self.output = None
        self.verbose = False
        self.debug = False
        self.groups = True
        self.threaded = False
        self.memory = False

        # Filtering
        self.include_stdlib = False
        self.include_pycallgraph = False
        self.max_depth = 99999

        self.trace_filter = GlobbingFilter(
            exclude=['pycallgraph.*'],
            include=['*'],
        )

        self.did_init = True

        # Update the defaults with anything from kwargs
        [setattr(self, k, v) for k, v in kwargs.items()]

        self.create_parser()

    def log_verbose(self, text):
        if self.verbose:
            print(text)

    def log_debug(self, text):
        if self.debug:
            print(text)

    def add_module_arguments(self, usage):
        subparsers = self.parser.add_subparsers(
            help='OUTPUT_TYPE', dest='output')
        parent_parser = self.create_parent_parser()

        for name, cls in list(outputters.items()):
            cls.add_arguments(subparsers, parent_parser, usage)

    def get_output(self):
        if not self.output:
            return
        output = outputters[self.output]()
        output.set_config(self)
        return output

    def parse_args(self, args=None):
        self.parser.parse_args(args, namespace=self)
        self.convert_filter_args()

    def strip_argv(self):
        sys.argv = [self.command] + self.command_args

    def convert_filter_args(self):
        if not self.include:
            self.include = ['*']

        if not self.include_pycallgraph:
            self.exclude.append('pycallgraph.*')

        self.trace_filter = GlobbingFilter(
            include=self.include,
            exclude=self.exclude,
        )

    def create_parser(self):
        '''Used by the pycallgraph command line interface to parse
        arguments.
        '''
        usage = 'pycallgraph [options] OUTPUT_TYPE [output_options] -- ' \
            'SCRIPT.py [ARG ...]'

        self.parser = argparse.ArgumentParser(
            description='Python Call Graph profiles a Python script and '
            'generates a call graph visualization.', usage=usage,
        )

        self.add_ungrouped_arguments()
        self.add_filter_arguments()
        self.add_module_arguments(usage)

    def create_parent_parser(self):
        '''Mixing subparsers with positional arguments can be done with a
        parents option. Found via: http://stackoverflow.com/a/11109863/11125
        '''
        parent_parser = argparse.ArgumentParser(add_help=False)
        parent_parser.add_argument(
            'command', metavar='SCRIPT',
            help='The Python script file to profile',
        )
        parent_parser.add_argument(
            'command_args', metavar='ARG', nargs='*',
            help='Python script arguments.'
        )
        return parent_parser

    def add_ungrouped_arguments(self):
        self.parser.add_argument(
            '-v', '--verbose', action='store_true', default=self.verbose,
            help='Display informative messages while running')

        self.parser.add_argument(
            '-d', '--debug', action='store_true', default=self.debug,
            help='Display debugging messages while running')

        self.parser.add_argument(
            '-t', '--threaded', action='store_true', default=self.threaded,
            help='Process traces asyncronously (Experimental)')

        self.parser.add_argument(
            '-ng', '--no-groups', dest='groups', action='store_false',
            default=self.groups, help='Do not group functions by module')

        self.parser.add_argument(
            '-s', '--stdlib', dest='include_stdlib', action='store_true',
            default=self.include_stdlib,
            help='Include standard library functions in the trace')

        self.parser.add_argument(
            '-m', '--memory', action='store_true', default=self.memory,
            help='(Experimental) Track memory usage')

    def add_filter_arguments(self):
        group = self.parser.add_argument_group('filtering')
        group.add_argument(
            '-i', '--include', default=[], action='append',
            help='Wildcard pattern of modules to include in the output. '
            'You can have multiple include arguments.'
        )

        group.add_argument(
            '-e', '--exclude', default=[], action='append',
            help='Wildcard pattern of modules to exclude in the output. '
            'You can have multiple exclude arguments.'
        )

        group.add_argument(
            '--include-pycallgraph', default=self.include_pycallgraph,
            action='store_true',
            help='Do not automatically filter out pycallgraph',
        )

        group.add_argument(
            '--max-depth', default=self.max_depth, type=int,
            help='Maximum stack depth to trace',
        )

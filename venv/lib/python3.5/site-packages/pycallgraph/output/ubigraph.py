try:
    from xmlrpc.client import Server
except ImportError:
    from xmlrpc.client import Server


# from ..exceptions import PyCallGraphException
from .output import Output


class UbigraphOutput(Output):

    def __init__(self, **kwargs):
        self.fp = None
        self.server_url = 'http://127.0.0.1:20738/RPC2'
        Output.__init__(self, **kwargs)

    def start(self):
        server = Server(self.server_url)
        self.graph = server.ubigraph

        # Create a graph
        for i in range(0, 10):
            self.graph.new_vertex_w_id(i)

        # Make some edges
        for i in range(0, 10):
            self.graph.new_edge(i, (i + 1) % 10)

    def should_update(self):
        return True

    def update(self):
        pass

    @classmethod
    def add_arguments(cls, subparsers, parent_parser, usage):
        defaults = cls()

        subparser = subparsers.add_parser(
            'ubigraph',
            help='Update an Ubigraph visualization in real time',
            parents=[parent_parser], usage=usage,
        )

        subparser.add_argument(
            '-s', '--server-url', type=str, default=defaults.server_url,
            help='The Ubigraph server',
        )

        return subparser

    def done(self):
        pass

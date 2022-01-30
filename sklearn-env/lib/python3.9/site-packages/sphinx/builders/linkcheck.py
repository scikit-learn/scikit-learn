"""
    sphinx.builders.linkcheck
    ~~~~~~~~~~~~~~~~~~~~~~~~~

    The CheckExternalLinksBuilder class.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import json
import re
import socket
import time
import warnings
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from html.parser import HTMLParser
from os import path
from queue import PriorityQueue, Queue
from threading import Thread
from typing import (Any, Dict, Generator, List, NamedTuple, Optional, Pattern, Set, Tuple,
                    Union, cast)
from urllib.parse import unquote, urlparse, urlunparse

from docutils import nodes
from docutils.nodes import Element
from requests import Response
from requests.exceptions import ConnectionError, HTTPError, TooManyRedirects

from sphinx.application import Sphinx
from sphinx.builders.dummy import DummyBuilder
from sphinx.config import Config
from sphinx.deprecation import RemovedInSphinx50Warning
from sphinx.environment import BuildEnvironment
from sphinx.locale import __
from sphinx.transforms.post_transforms import SphinxPostTransform
from sphinx.util import encode_uri, logging, requests
from sphinx.util.console import darkgray, darkgreen, purple, red, turquoise  # type: ignore
from sphinx.util.nodes import get_node_line

logger = logging.getLogger(__name__)

uri_re = re.compile('([a-z]+:)?//')  # matches to foo:// and // (a protocol relative URL)


class Hyperlink(NamedTuple):
    uri: str
    docname: str
    lineno: Optional[int]


class CheckRequest(NamedTuple):
    next_check: float
    hyperlink: Optional[Hyperlink]


class CheckResult(NamedTuple):
    uri: str
    docname: str
    lineno: int
    status: str
    message: str
    code: int


class RateLimit(NamedTuple):
    delay: float
    next_check: float


# Tuple is old styled CheckRequest
CheckRequestType = Union[CheckRequest, Tuple[float, str, str, int]]

DEFAULT_REQUEST_HEADERS = {
    'Accept': 'text/html,application/xhtml+xml;q=0.9,*/*;q=0.8',
}
CHECK_IMMEDIATELY = 0
QUEUE_POLL_SECS = 1
DEFAULT_DELAY = 60.0


def node_line_or_0(node: Element) -> int:
    """
    PriorityQueue items must be comparable. The line number is part of the
    tuple used by the PriorityQueue, keep an homogeneous type for comparison.
    """
    warnings.warn('node_line_or_0() is deprecated.',
                  RemovedInSphinx50Warning, stacklevel=2)
    return get_node_line(node) or 0


class AnchorCheckParser(HTMLParser):
    """Specialized HTML parser that looks for a specific anchor."""

    def __init__(self, search_anchor: str) -> None:
        super().__init__()

        self.search_anchor = search_anchor
        self.found = False

    def handle_starttag(self, tag: Any, attrs: Any) -> None:
        for key, value in attrs:
            if key in ('id', 'name') and value == self.search_anchor:
                self.found = True
                break


def check_anchor(response: requests.requests.Response, anchor: str) -> bool:
    """Reads HTML data from a response object `response` searching for `anchor`.
    Returns True if anchor was found, False otherwise.
    """
    parser = AnchorCheckParser(anchor)
    # Read file in chunks. If we find a matching anchor, we break
    # the loop early in hopes not to have to download the whole thing.
    for chunk in response.iter_content(chunk_size=4096, decode_unicode=True):
        if isinstance(chunk, bytes):    # requests failed to decode
            chunk = chunk.decode()      # manually try to decode it

        parser.feed(chunk)
        if parser.found:
            break
    parser.close()
    return parser.found


class CheckExternalLinksBuilder(DummyBuilder):
    """
    Checks for broken external links.
    """
    name = 'linkcheck'
    epilog = __('Look for any errors in the above output or in '
                '%(outdir)s/output.txt')

    def init(self) -> None:
        self.hyperlinks: Dict[str, Hyperlink] = {}
        self._good: Set[str] = set()
        self._broken: Dict[str, str] = {}
        self._redirected: Dict[str, Tuple[str, int]] = {}
        # set a timeout for non-responding servers
        socket.setdefaulttimeout(5.0)

        # create queues and worker threads
        self._wqueue: PriorityQueue[CheckRequestType] = PriorityQueue()
        self._rqueue: Queue[CheckResult] = Queue()

    @property
    def anchors_ignore(self) -> List[Pattern]:
        warnings.warn(
            "%s.%s is deprecated." % (self.__class__.__name__, "anchors_ignore"),
            RemovedInSphinx50Warning,
            stacklevel=2,
        )
        return [re.compile(x) for x in self.config.linkcheck_anchors_ignore]

    @property
    def auth(self) -> List[Tuple[Pattern, Any]]:
        warnings.warn(
            "%s.%s is deprecated." % (self.__class__.__name__, "auth"),
            RemovedInSphinx50Warning,
            stacklevel=2,
        )
        return [(re.compile(pattern), auth_info) for pattern, auth_info
                in self.config.linkcheck_auth]

    @property
    def to_ignore(self) -> List[Pattern]:
        warnings.warn(
            "%s.%s is deprecated." % (self.__class__.__name__, "to_ignore"),
            RemovedInSphinx50Warning,
            stacklevel=2,
        )
        return [re.compile(x) for x in self.config.linkcheck_ignore]

    @property
    def good(self) -> Set[str]:
        warnings.warn(
            "%s.%s is deprecated." % (self.__class__.__name__, "good"),
            RemovedInSphinx50Warning,
            stacklevel=2,
        )
        return self._good

    @property
    def broken(self) -> Dict[str, str]:
        warnings.warn(
            "%s.%s is deprecated." % (self.__class__.__name__, "broken"),
            RemovedInSphinx50Warning,
            stacklevel=2,
        )
        return self._broken

    @property
    def redirected(self) -> Dict[str, Tuple[str, int]]:
        warnings.warn(
            "%s.%s is deprecated." % (self.__class__.__name__, "redirected"),
            RemovedInSphinx50Warning,
            stacklevel=2,
        )
        return self._redirected

    def check_thread(self) -> None:
        warnings.warn(
            "%s.%s is deprecated." % (self.__class__.__name__, "check_thread"),
            RemovedInSphinx50Warning,
            stacklevel=2,
        )
        # do nothing.

    def limit_rate(self, response: Response) -> Optional[float]:
        warnings.warn(
            "%s.%s is deprecated." % (self.__class__.__name__, "limit_rate"),
            RemovedInSphinx50Warning,
            stacklevel=2,
        )
        worker = HyperlinkAvailabilityCheckWorker(self.env, self.config,
                                                  None, None, {})
        return worker.limit_rate(response)

    def rqueue(self, response: Response) -> Queue:
        warnings.warn(
            "%s.%s is deprecated." % (self.__class__.__name__, "rqueue"),
            RemovedInSphinx50Warning,
            stacklevel=2,
        )
        return self._rqueue

    def workers(self, response: Response) -> List[Thread]:
        warnings.warn(
            "%s.%s is deprecated." % (self.__class__.__name__, "workers"),
            RemovedInSphinx50Warning,
            stacklevel=2,
        )
        return []

    def wqueue(self, response: Response) -> Queue:
        warnings.warn(
            "%s.%s is deprecated." % (self.__class__.__name__, "wqueue"),
            RemovedInSphinx50Warning,
            stacklevel=2,
        )
        return self._wqueue

    def process_result(self, result: CheckResult) -> None:
        filename = self.env.doc2path(result.docname, None)

        linkstat = dict(filename=filename, lineno=result.lineno,
                        status=result.status, code=result.code, uri=result.uri,
                        info=result.message)
        self.write_linkstat(linkstat)

        if result.status == 'unchecked':
            return
        if result.status == 'working' and result.message == 'old':
            return
        if result.lineno:
            logger.info('(%16s: line %4d) ', result.docname, result.lineno, nonl=True)
        if result.status == 'ignored':
            if result.message:
                logger.info(darkgray('-ignored- ') + result.uri + ': ' + result.message)
            else:
                logger.info(darkgray('-ignored- ') + result.uri)
        elif result.status == 'local':
            logger.info(darkgray('-local-   ') + result.uri)
            self.write_entry('local', result.docname, filename, result.lineno, result.uri)
        elif result.status == 'working':
            logger.info(darkgreen('ok        ') + result.uri + result.message)
        elif result.status == 'broken':
            if self.app.quiet or self.app.warningiserror:
                logger.warning(__('broken link: %s (%s)'), result.uri, result.message,
                               location=(result.docname, result.lineno))
            else:
                logger.info(red('broken    ') + result.uri + red(' - ' + result.message))
            self.write_entry('broken', result.docname, filename, result.lineno,
                             result.uri + ': ' + result.message)
        elif result.status == 'redirected':
            try:
                text, color = {
                    301: ('permanently', purple),
                    302: ('with Found', purple),
                    303: ('with See Other', purple),
                    307: ('temporarily', turquoise),
                    308: ('permanently', purple),
                }[result.code]
            except KeyError:
                text, color = ('with unknown code', purple)
            linkstat['text'] = text
            if self.config.linkcheck_allowed_redirects:
                logger.warning('redirect  ' + result.uri + ' - ' + text + ' to ' +
                               result.message, location=(result.docname, result.lineno))
            else:
                logger.info(color('redirect  ') + result.uri +
                            color(' - ' + text + ' to ' + result.message))
            self.write_entry('redirected ' + text, result.docname, filename,
                             result.lineno, result.uri + ' to ' + result.message)
        else:
            raise ValueError("Unknown status %s." % result.status)

    def write_entry(self, what: str, docname: str, filename: str, line: int,
                    uri: str) -> None:
        self.txt_outfile.write("%s:%s: [%s] %s\n" % (filename, line, what, uri))

    def write_linkstat(self, data: dict) -> None:
        self.json_outfile.write(json.dumps(data))
        self.json_outfile.write('\n')

    def finish(self) -> None:
        checker = HyperlinkAvailabilityChecker(self.env, self.config, self)
        logger.info('')

        with open(path.join(self.outdir, 'output.txt'), 'w') as self.txt_outfile,\
             open(path.join(self.outdir, 'output.json'), 'w') as self.json_outfile:
            for result in checker.check(self.hyperlinks):
                self.process_result(result)

        if self._broken:
            self.app.statuscode = 1


class HyperlinkAvailabilityChecker:
    def __init__(self, env: BuildEnvironment, config: Config,
                 builder: CheckExternalLinksBuilder = None) -> None:
        # Warning: builder argument will be removed in the sphinx-5.0.
        # Don't use it from extensions.
        # tag: RemovedInSphinx50Warning
        self.builder = builder
        self.config = config
        self.env = env
        self.rate_limits: Dict[str, RateLimit] = {}
        self.workers: List[Thread] = []

        self.to_ignore = [re.compile(x) for x in self.config.linkcheck_ignore]

        if builder:
            self.rqueue = builder._rqueue
            self.wqueue = builder._wqueue
        else:
            self.rqueue = Queue()
            self.wqueue = PriorityQueue()

    def invoke_threads(self) -> None:
        for _i in range(self.config.linkcheck_workers):
            thread = HyperlinkAvailabilityCheckWorker(self.env, self.config,
                                                      self.rqueue, self.wqueue,
                                                      self.rate_limits, self.builder)
            thread.start()
            self.workers.append(thread)

    def shutdown_threads(self) -> None:
        self.wqueue.join()
        for _worker in self.workers:
            self.wqueue.put(CheckRequest(CHECK_IMMEDIATELY, None), False)

    def check(self, hyperlinks: Dict[str, Hyperlink]) -> Generator[CheckResult, None, None]:
        self.invoke_threads()

        total_links = 0
        for hyperlink in hyperlinks.values():
            if self.is_ignored_uri(hyperlink.uri):
                yield CheckResult(hyperlink.uri, hyperlink.docname, hyperlink.lineno,
                                  'ignored', '', 0)
            else:
                self.wqueue.put(CheckRequest(CHECK_IMMEDIATELY, hyperlink), False)
                total_links += 1

        done = 0
        while done < total_links:
            yield self.rqueue.get()
            done += 1

        self.shutdown_threads()

    def is_ignored_uri(self, uri: str) -> bool:
        return any(pat.match(uri) for pat in self.to_ignore)


class HyperlinkAvailabilityCheckWorker(Thread):
    """A worker class for checking the availability of hyperlinks."""

    def __init__(self, env: BuildEnvironment, config: Config, rqueue: Queue,
                 wqueue: Queue, rate_limits: Dict[str, RateLimit],
                 builder: CheckExternalLinksBuilder = None) -> None:
        # Warning: builder argument will be removed in the sphinx-5.0.
        # Don't use it from extensions.
        # tag: RemovedInSphinx50Warning
        self.config = config
        self.env = env
        self.rate_limits = rate_limits
        self.rqueue = rqueue
        self.wqueue = wqueue

        self.anchors_ignore = [re.compile(x)
                               for x in self.config.linkcheck_anchors_ignore]
        self.documents_exclude = [re.compile(doc)
                                  for doc in self.config.linkcheck_exclude_documents]
        self.auth = [(re.compile(pattern), auth_info) for pattern, auth_info
                     in self.config.linkcheck_auth]

        if builder:
            # if given, fill the result of checks as cache
            self._good = builder._good
            self._broken = builder._broken
            self._redirected = builder._redirected
        else:
            # only for compatibility. Will be removed in Sphinx-5.0
            self._good = set()
            self._broken = {}
            self._redirected = {}

        super().__init__(daemon=True)

    def run(self) -> None:
        kwargs = {}
        if self.config.linkcheck_timeout:
            kwargs['timeout'] = self.config.linkcheck_timeout

        def get_request_headers() -> Dict:
            url = urlparse(uri)
            candidates = ["%s://%s" % (url.scheme, url.netloc),
                          "%s://%s/" % (url.scheme, url.netloc),
                          uri,
                          "*"]

            for u in candidates:
                if u in self.config.linkcheck_request_headers:
                    headers = dict(DEFAULT_REQUEST_HEADERS)
                    headers.update(self.config.linkcheck_request_headers[u])
                    return headers

            return {}

        def check_uri() -> Tuple[str, str, int]:
            # split off anchor
            if '#' in uri:
                req_url, anchor = uri.split('#', 1)
                for rex in self.anchors_ignore:
                    if rex.match(anchor):
                        anchor = None
                        break
            else:
                req_url = uri
                anchor = None

            # handle non-ASCII URIs
            try:
                req_url.encode('ascii')
            except UnicodeError:
                req_url = encode_uri(req_url)

            # Get auth info, if any
            for pattern, auth_info in self.auth:
                if pattern.match(uri):
                    break
            else:
                auth_info = None

            # update request headers for the URL
            kwargs['headers'] = get_request_headers()

            try:
                if anchor and self.config.linkcheck_anchors:
                    # Read the whole document and see if #anchor exists
                    response = requests.get(req_url, stream=True, config=self.config,
                                            auth=auth_info, **kwargs)
                    response.raise_for_status()
                    found = check_anchor(response, unquote(anchor))

                    if not found:
                        raise Exception(__("Anchor '%s' not found") % anchor)
                else:
                    try:
                        # try a HEAD request first, which should be easier on
                        # the server and the network
                        response = requests.head(req_url, allow_redirects=True,
                                                 config=self.config, auth=auth_info,
                                                 **kwargs)
                        response.raise_for_status()
                    # Servers drop the connection on HEAD requests, causing
                    # ConnectionError.
                    except (ConnectionError, HTTPError, TooManyRedirects) as err:
                        if isinstance(err, HTTPError) and err.response.status_code == 429:
                            raise
                        # retry with GET request if that fails, some servers
                        # don't like HEAD requests.
                        response = requests.get(req_url, stream=True,
                                                config=self.config,
                                                auth=auth_info, **kwargs)
                        response.raise_for_status()
            except HTTPError as err:
                if err.response.status_code == 401:
                    # We'll take "Unauthorized" as working.
                    return 'working', ' - unauthorized', 0
                elif err.response.status_code == 429:
                    next_check = self.limit_rate(err.response)
                    if next_check is not None:
                        self.wqueue.put(CheckRequest(next_check, hyperlink), False)
                        return 'rate-limited', '', 0
                    return 'broken', str(err), 0
                elif err.response.status_code == 503:
                    # We'll take "Service Unavailable" as ignored.
                    return 'ignored', str(err), 0
                else:
                    return 'broken', str(err), 0
            except Exception as err:
                return 'broken', str(err), 0
            else:
                netloc = urlparse(req_url).netloc
                try:
                    del self.rate_limits[netloc]
                except KeyError:
                    pass
            if response.url.rstrip('/') == req_url.rstrip('/'):
                return 'working', '', 0
            else:
                new_url = response.url
                if anchor:
                    new_url += '#' + anchor

                if allowed_redirect(req_url, new_url):
                    return 'working', '', 0
                elif response.history:
                    # history contains any redirects, get last
                    code = response.history[-1].status_code
                    return 'redirected', new_url, code
                else:
                    return 'redirected', new_url, 0

        def allowed_redirect(url: str, new_url: str) -> bool:
            for from_url, to_url in self.config.linkcheck_allowed_redirects.items():
                if from_url.match(url) and to_url.match(new_url):
                    return True

            return False

        def check(docname: str) -> Tuple[str, str, int]:
            # check for various conditions without bothering the network

            for doc_matcher in self.documents_exclude:
                if doc_matcher.match(docname):
                    info = (
                        f'{docname} matched {doc_matcher.pattern} from '
                        'linkcheck_exclude_documents'
                    )
                    return 'ignored', info, 0

            if len(uri) == 0 or uri.startswith(('#', 'mailto:', 'tel:')):
                return 'unchecked', '', 0
            elif not uri.startswith(('http:', 'https:')):
                if uri_re.match(uri):
                    # non supported URI schemes (ex. ftp)
                    return 'unchecked', '', 0
                else:
                    srcdir = path.dirname(self.env.doc2path(docname))
                    if path.exists(path.join(srcdir, uri)):
                        return 'working', '', 0
                    else:
                        self._broken[uri] = ''
                        return 'broken', '', 0
            elif uri in self._good:
                return 'working', 'old', 0
            elif uri in self._broken:
                return 'broken', self._broken[uri], 0
            elif uri in self._redirected:
                return 'redirected', self._redirected[uri][0], self._redirected[uri][1]

            # need to actually check the URI
            for _ in range(self.config.linkcheck_retries):
                status, info, code = check_uri()
                if status != "broken":
                    break

            if status == "working":
                self._good.add(uri)
            elif status == "broken":
                self._broken[uri] = info
            elif status == "redirected":
                self._redirected[uri] = (info, code)

            return (status, info, code)

        while True:
            check_request = self.wqueue.get()
            try:
                next_check, hyperlink = check_request
                if hyperlink is None:
                    break

                uri, docname, lineno = hyperlink
            except ValueError:
                # old styled check_request (will be deprecated in Sphinx-5.0)
                next_check, uri, docname, lineno = check_request

            if uri is None:
                break
            netloc = urlparse(uri).netloc
            try:
                # Refresh rate limit.
                # When there are many links in the queue, workers are all stuck waiting
                # for responses, but the builder keeps queuing. Links in the queue may
                # have been queued before rate limits were discovered.
                next_check = self.rate_limits[netloc].next_check
            except KeyError:
                pass
            if next_check > time.time():
                # Sleep before putting message back in the queue to avoid
                # waking up other threads.
                time.sleep(QUEUE_POLL_SECS)
                self.wqueue.put(CheckRequest(next_check, hyperlink), False)
                self.wqueue.task_done()
                continue
            status, info, code = check(docname)
            if status == 'rate-limited':
                logger.info(darkgray('-rate limited-   ') + uri + darkgray(' | sleeping...'))
            else:
                self.rqueue.put(CheckResult(uri, docname, lineno, status, info, code))
            self.wqueue.task_done()

    def limit_rate(self, response: Response) -> Optional[float]:
        next_check = None
        retry_after = response.headers.get("Retry-After")
        if retry_after:
            try:
                # Integer: time to wait before next attempt.
                delay = float(retry_after)
            except ValueError:
                try:
                    # An HTTP-date: time of next attempt.
                    until = parsedate_to_datetime(retry_after)
                except (TypeError, ValueError):
                    # TypeError: Invalid date format.
                    # ValueError: Invalid date, e.g. Oct 52th.
                    pass
                else:
                    next_check = datetime.timestamp(until)
                    delay = (until - datetime.now(timezone.utc)).total_seconds()
            else:
                next_check = time.time() + delay
        netloc = urlparse(response.url).netloc
        if next_check is None:
            max_delay = self.config.linkcheck_rate_limit_timeout
            try:
                rate_limit = self.rate_limits[netloc]
            except KeyError:
                delay = DEFAULT_DELAY
            else:
                last_wait_time = rate_limit.delay
                delay = 2.0 * last_wait_time
                if delay > max_delay and last_wait_time < max_delay:
                    delay = max_delay
            if delay > max_delay:
                return None
            next_check = time.time() + delay
        self.rate_limits[netloc] = RateLimit(delay, next_check)
        return next_check


class HyperlinkCollector(SphinxPostTransform):
    builders = ('linkcheck',)
    default_priority = 800

    def run(self, **kwargs: Any) -> None:
        builder = cast(CheckExternalLinksBuilder, self.app.builder)
        hyperlinks = builder.hyperlinks

        # reference nodes
        for refnode in self.document.findall(nodes.reference):
            if 'refuri' not in refnode:
                continue
            uri = refnode['refuri']
            newuri = self.app.emit_firstresult('linkcheck-process-uri', uri)
            if newuri:
                uri = newuri

            lineno = get_node_line(refnode)
            uri_info = Hyperlink(uri, self.env.docname, lineno)
            if uri not in hyperlinks:
                hyperlinks[uri] = uri_info

        # image nodes
        for imgnode in self.document.findall(nodes.image):
            uri = imgnode['candidates'].get('?')
            if uri and '://' in uri:
                newuri = self.app.emit_firstresult('linkcheck-process-uri', uri)
                if newuri:
                    uri = newuri

                lineno = get_node_line(imgnode)
                uri_info = Hyperlink(uri, self.env.docname, lineno)
                if uri not in hyperlinks:
                    hyperlinks[uri] = uri_info


def rewrite_github_anchor(app: Sphinx, uri: str) -> Optional[str]:
    """Rewrite anchor name of the hyperlink to github.com

    The hyperlink anchors in github.com are dynamically generated.  This rewrites
    them before checking and makes them comparable.
    """
    parsed = urlparse(uri)
    if parsed.hostname == "github.com" and parsed.fragment:
        prefixed = parsed.fragment.startswith('user-content-')
        if not prefixed:
            fragment = f'user-content-{parsed.fragment}'
            return urlunparse(parsed._replace(fragment=fragment))
    return None


def compile_linkcheck_allowed_redirects(app: Sphinx, config: Config) -> None:
    """Compile patterns in linkcheck_allowed_redirects to the regexp objects."""
    for url, pattern in list(app.config.linkcheck_allowed_redirects.items()):
        try:
            app.config.linkcheck_allowed_redirects[re.compile(url)] = re.compile(pattern)
        except re.error as exc:
            logger.warning(__('Failed to compile regex in linkcheck_allowed_redirects: %r %s'),
                           exc.pattern, exc.msg)
        finally:
            # Remove the original regexp-string
            app.config.linkcheck_allowed_redirects.pop(url)


def setup(app: Sphinx) -> Dict[str, Any]:
    app.add_builder(CheckExternalLinksBuilder)
    app.add_post_transform(HyperlinkCollector)

    app.add_config_value('linkcheck_ignore', [], None)
    app.add_config_value('linkcheck_exclude_documents', [], None)
    app.add_config_value('linkcheck_allowed_redirects', {}, None)
    app.add_config_value('linkcheck_auth', [], None)
    app.add_config_value('linkcheck_request_headers', {}, None)
    app.add_config_value('linkcheck_retries', 1, None)
    app.add_config_value('linkcheck_timeout', None, None, [int])
    app.add_config_value('linkcheck_workers', 5, None)
    app.add_config_value('linkcheck_anchors', True, None)
    # Anchors starting with ! are ignored since they are
    # commonly used for dynamic pages
    app.add_config_value('linkcheck_anchors_ignore', ["^!"], None)
    app.add_config_value('linkcheck_rate_limit_timeout', 300.0, None)

    app.add_event('linkcheck-process-uri')

    app.connect('config-inited', compile_linkcheck_allowed_redirects, priority=800)

    # FIXME: Disable URL rewrite handler for github.com temporarily.
    # ref: https://github.com/sphinx-doc/sphinx/issues/9435
    # app.connect('linkcheck-process-uri', rewrite_github_anchor)

    return {
        'version': 'builtin',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }

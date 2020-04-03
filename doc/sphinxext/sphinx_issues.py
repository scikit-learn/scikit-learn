# -*- coding: utf-8 -*-
"""A Sphinx extension for linking to your project's issue tracker.

Copyright 2014 Steven Loria

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""
import re

from docutils import nodes, utils
from sphinx.util.nodes import split_explicit_title

__version__ = "1.2.0"
__author__ = "Steven Loria"
__license__ = "MIT"


def user_role(name, rawtext, text, lineno, inliner, options=None, content=None):
    """Sphinx role for linking to a user profile. Defaults to linking to
    Github profiles, but the profile URIS can be configured via the
    ``issues_user_uri`` config value.
    Examples: ::
        :user:`sloria`
    Anchor text also works: ::
        :user:`Steven Loria <sloria>`
    """
    options = options or {}
    content = content or []
    has_explicit_title, title, target = split_explicit_title(text)

    target = utils.unescape(target).strip()
    title = utils.unescape(title).strip()
    config = inliner.document.settings.env.app.config
    if config.issues_user_uri:
        ref = config.issues_user_uri.format(user=target)
    else:
        ref = "https://github.com/{0}".format(target)
    if has_explicit_title:
        text = title
    else:
        text = "@{0}".format(target)

    link = nodes.reference(text=text, refuri=ref, **options)
    return [link], []


def cve_role(name, rawtext, text, lineno, inliner, options=None, content=None):
    """Sphinx role for linking to a CVE on https://cve.mitre.org.
    Examples: ::
        :cve:`CVE-2018-17175`
    """
    options = options or {}
    content = content or []
    has_explicit_title, title, target = split_explicit_title(text)

    target = utils.unescape(target).strip()
    title = utils.unescape(title).strip()
    ref = "https://cve.mitre.org/cgi-bin/cvename.cgi?name={0}".format(target)
    text = title if has_explicit_title else target
    link = nodes.reference(text=text, refuri=ref, **options)
    return [link], []


class IssueRole(object):

    EXTERNAL_REPO_REGEX = re.compile(r"^(\w+)/(.+)([#@])([\w]+)$")

    def __init__(
        self, uri_config_option, format_kwarg, github_uri_template, format_text=None
    ):
        self.uri_config_option = uri_config_option
        self.format_kwarg = format_kwarg
        self.github_uri_template = github_uri_template
        self.format_text = format_text or self.default_format_text

    @staticmethod
    def default_format_text(issue_no):
        return "#{0}".format(issue_no)

    def make_node(self, name, issue_no, config, options=None):
        name_map = {"pr": "pull", "issue": "issues", "commit": "commit"}
        options = options or {}
        repo_match = self.EXTERNAL_REPO_REGEX.match(issue_no)
        if repo_match:  # External repo
            username, repo, symbol, issue = repo_match.groups()
            if name not in name_map:
                raise ValueError(
                    "External repo linking not supported for :{}:".format(name)
                )
            path = name_map.get(name)
            ref = "https://github.com/{issues_github_path}/{path}/{n}".format(
                issues_github_path="{}/{}".format(username, repo), path=path, n=issue
            )
            formatted_issue = self.format_text(issue).lstrip("#")
            text = "{username}/{repo}{symbol}{formatted_issue}".format(**locals())
            link = nodes.reference(text=text, refuri=ref, **options)
            return link

        if issue_no not in ("-", "0"):
            uri_template = getattr(config, self.uri_config_option, None)
            if uri_template:
                ref = uri_template.format(**{self.format_kwarg: issue_no})
            elif config.issues_github_path:
                ref = self.github_uri_template.format(
                    issues_github_path=config.issues_github_path, n=issue_no
                )
            else:
                raise ValueError(
                    "Neither {} nor issues_github_path "
                    "is set".format(self.uri_config_option)
                )
            issue_text = self.format_text(issue_no)
            link = nodes.reference(text=issue_text, refuri=ref, **options)
        else:
            link = None
        return link

    def __call__(
        self, name, rawtext, text, lineno, inliner, options=None, content=None
    ):
        options = options or {}
        content = content or []
        issue_nos = [each.strip() for each in utils.unescape(text).split(",")]
        config = inliner.document.settings.env.app.config
        ret = []
        for i, issue_no in enumerate(issue_nos):
            node = self.make_node(name, issue_no, config, options=options)
            ret.append(node)
            if i != len(issue_nos) - 1:
                sep = nodes.raw(text=", ", format="html")
                ret.append(sep)
        return ret, []


"""Sphinx role for linking to an issue. Must have
`issues_uri` or `issues_github_path` configured in ``conf.py``.
Examples: ::
    :issue:`123`
    :issue:`42,45`
    :issue:`sloria/konch#123`
"""
issue_role = IssueRole(
    uri_config_option="issues_uri",
    format_kwarg="issue",
    github_uri_template="https://github.com/{issues_github_path}/issues/{n}",
)

"""Sphinx role for linking to a pull request. Must have
`issues_pr_uri` or `issues_github_path` configured in ``conf.py``.
Examples: ::
    :pr:`123`
    :pr:`42,45`
    :pr:`sloria/konch#43`
"""
pr_role = IssueRole(
    uri_config_option="issues_pr_uri",
    format_kwarg="pr",
    github_uri_template="https://github.com/{issues_github_path}/pull/{n}",
)


def format_commit_text(sha):
    return sha[:7]


"""Sphinx role for linking to a commit. Must have
`issues_pr_uri` or `issues_github_path` configured in ``conf.py``.
Examples: ::
    :commit:`123abc456def`
    :commit:`sloria/konch@123abc456def`
"""
commit_role = IssueRole(
    uri_config_option="issues_commit_uri",
    format_kwarg="commit",
    github_uri_template="https://github.com/{issues_github_path}/commit/{n}",
    format_text=format_commit_text,
)


def setup(app):
    # Format template for issues URI
    # e.g. 'https://github.com/sloria/marshmallow/issues/{issue}
    app.add_config_value("issues_uri", default=None, rebuild="html")
    # Format template for PR URI
    # e.g. 'https://github.com/sloria/marshmallow/pull/{issue}
    app.add_config_value("issues_pr_uri", default=None, rebuild="html")
    # Format template for commit URI
    # e.g. 'https://github.com/sloria/marshmallow/commits/{commit}
    app.add_config_value("issues_commit_uri", default=None, rebuild="html")
    # Shortcut for Github, e.g. 'sloria/marshmallow'
    app.add_config_value("issues_github_path", default=None, rebuild="html")
    # Format template for user profile URI
    # e.g. 'https://github.com/{user}'
    app.add_config_value("issues_user_uri", default=None, rebuild="html")
    app.add_role("issue", issue_role)
    app.add_role("pr", pr_role)
    app.add_role("user", user_role)
    app.add_role("commit", commit_role)
    app.add_role("cve", cve_role)
    return {
        "version": __version__,
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }

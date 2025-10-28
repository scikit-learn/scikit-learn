"""This module contains the code for loading intersphinx inventories."""

from __future__ import annotations

import concurrent.futures
import posixpath
import time
from operator import itemgetter
from os import path
from typing import TYPE_CHECKING
from urllib.parse import urlsplit, urlunsplit

from sphinx.builders.html import INVENTORY_FILENAME
from sphinx.errors import ConfigError
from sphinx.ext.intersphinx._shared import LOGGER, InventoryAdapter, _IntersphinxProject
from sphinx.locale import __
from sphinx.util import requests
from sphinx.util.inventory import InventoryFile

if TYPE_CHECKING:
    from pathlib import Path

    from urllib3.response import HTTPResponse

    from sphinx.application import Sphinx
    from sphinx.config import Config
    from sphinx.ext.intersphinx._shared import (
        IntersphinxMapping,
        InventoryCacheEntry,
        InventoryLocation,
        InventoryName,
        InventoryURI,
    )
    from sphinx.util.typing import Inventory, _ReadableStream


def validate_intersphinx_mapping(app: Sphinx, config: Config) -> None:
    """Validate and normalise :confval:`intersphinx_mapping`.

    Ensure that:

    * Keys are non-empty strings.
    * Values are two-element tuples or lists.
    * The first element of each value pair (the target URI)
      is a non-empty string.
    * The second element of each value pair (inventory locations)
      is a tuple of non-empty strings or None.
    """
    # URIs should NOT be duplicated, otherwise different builds may use
    # different project names (and thus, the build are no more reproducible)
    # depending on which one is inserted last in the cache.
    seen: dict[InventoryURI, InventoryName] = {}

    errors = 0
    for name, value in config.intersphinx_mapping.copy().items():
        # ensure that intersphinx projects are always named
        if not isinstance(name, str) or not name:
            errors += 1
            msg = __(
                'Invalid intersphinx project identifier `%r` in intersphinx_mapping. '
                'Project identifiers must be non-empty strings.'
            )
            LOGGER.error(msg, name)
            del config.intersphinx_mapping[name]
            continue

        # ensure values are properly formatted
        if not isinstance(value, (tuple | list)):
            errors += 1
            msg = __(
                'Invalid value `%r` in intersphinx_mapping[%r]. '
                'Expected a two-element tuple or list.'
            )
            LOGGER.error(msg, value, name)
            del config.intersphinx_mapping[name]
            continue
        try:
            uri, inv = value
        except (TypeError, ValueError, Exception):
            errors += 1
            msg = __(
                'Invalid value `%r` in intersphinx_mapping[%r]. '
                'Values must be a (target URI, inventory locations) pair.'
            )
            LOGGER.error(msg, value, name)
            del config.intersphinx_mapping[name]
            continue

        # ensure target URIs are non-empty and unique
        if not uri or not isinstance(uri, str):
            errors += 1
            msg = __(
                'Invalid target URI value `%r` in intersphinx_mapping[%r][0]. '
                'Target URIs must be unique non-empty strings.'
            )
            LOGGER.error(msg, uri, name)
            del config.intersphinx_mapping[name]
            continue
        if uri in seen:
            errors += 1
            msg = __(
                'Invalid target URI value `%r` in intersphinx_mapping[%r][0]. '
                'Target URIs must be unique (other instance in intersphinx_mapping[%r]).'
            )
            LOGGER.error(msg, uri, name, seen[uri])
            del config.intersphinx_mapping[name]
            continue
        seen[uri] = name

        if not isinstance(inv, tuple | list):
            inv = (inv,)

        # ensure inventory locations are None or non-empty
        targets: list[InventoryLocation] = []
        for target in inv:
            if target is None or target and isinstance(target, str):
                targets.append(target)
            else:
                errors += 1
                msg = __(
                    'Invalid inventory location value `%r` in intersphinx_mapping[%r][1]. '
                    'Inventory locations must be non-empty strings or None.'
                )
                LOGGER.error(msg, target, name)
                del config.intersphinx_mapping[name]
                continue

        config.intersphinx_mapping[name] = (name, (uri, tuple(targets)))

    if errors == 1:
        msg = __('Invalid `intersphinx_mapping` configuration (1 error).')
        raise ConfigError(msg)
    if errors > 1:
        msg = __('Invalid `intersphinx_mapping` configuration (%s errors).')
        raise ConfigError(msg % errors)


def load_mappings(app: Sphinx) -> None:
    """Load all intersphinx mappings into the environment.

    The intersphinx mappings are expected to be normalized.
    """
    now = int(time.time())
    inventories = InventoryAdapter(app.builder.env)
    intersphinx_cache: dict[InventoryURI, InventoryCacheEntry] = inventories.cache
    intersphinx_mapping: IntersphinxMapping = app.config.intersphinx_mapping

    projects = []
    for name, (uri, locations) in intersphinx_mapping.values():
        try:
            project = _IntersphinxProject(
                name=name, target_uri=uri, locations=locations
            )
        except ValueError as err:
            msg = __(
                'An invalid intersphinx_mapping entry was added after normalisation.'
            )
            raise ConfigError(msg) from err
        else:
            projects.append(project)

    expected_uris = {project.target_uri for project in projects}
    for uri in frozenset(intersphinx_cache):
        if intersphinx_cache[uri][0] not in intersphinx_mapping:
            # Remove all cached entries that are no longer in `intersphinx_mapping`.
            del intersphinx_cache[uri]
        elif uri not in expected_uris:
            # Remove cached entries with a different target URI
            # than the one in `intersphinx_mapping`.
            # This happens when the URI in `intersphinx_mapping` is changed.
            del intersphinx_cache[uri]

    with concurrent.futures.ThreadPoolExecutor() as pool:
        futures = [
            pool.submit(
                _fetch_inventory_group,
                project=project,
                cache=intersphinx_cache,
                now=now,
                config=app.config,
                srcdir=app.srcdir,
            )
            for project in projects
        ]
        updated = [f.result() for f in concurrent.futures.as_completed(futures)]

    if any(updated):
        # clear the local inventories
        inventories.clear()

        # Duplicate values in different inventories will shadow each
        # other; which one will override which can vary between builds.
        #
        # In an attempt to make this more consistent,
        # we sort the named inventories in the cache
        # by their name and expiry time ``(NAME, EXPIRY)``.
        by_name_and_time = itemgetter(0, 1)  # 0: name, 1: expiry
        cache_values = sorted(intersphinx_cache.values(), key=by_name_and_time)
        for name, _expiry, invdata in cache_values:
            inventories.named_inventory[name] = invdata
            for objtype, objects in invdata.items():
                inventories.main_inventory.setdefault(objtype, {}).update(objects)


def _fetch_inventory_group(
    *,
    project: _IntersphinxProject,
    cache: dict[InventoryURI, InventoryCacheEntry],
    now: int,
    config: Config,
    srcdir: Path,
) -> bool:
    if config.intersphinx_cache_limit >= 0:
        # Positive value: cache is expired if its timestamp is below
        # `now - X days`.
        cache_time = now - config.intersphinx_cache_limit * 86400
    else:
        # Negative value: cache is expired if its timestamp is below
        # zero, which is impossible.
        cache_time = 0

    updated = False
    failures = []

    for location in project.locations:
        # location is either None or a non-empty string
        if location is None:
            inv = posixpath.join(project.target_uri, INVENTORY_FILENAME)
        else:
            inv = location

        # decide whether the inventory must be read: always read local
        # files; remote ones only if the cache time is expired
        if (
            '://' not in inv
            or project.target_uri not in cache
            or cache[project.target_uri][1] < cache_time
        ):
            LOGGER.info(
                __("loading intersphinx inventory '%s' from %s ..."),
                project.name,
                _get_safe_url(inv),
            )

            try:
                invdata = _fetch_inventory(
                    target_uri=project.target_uri,
                    inv_location=inv,
                    config=config,
                    srcdir=srcdir,
                )
            except Exception as err:
                failures.append(err.args)
                continue

            if invdata:
                cache[project.target_uri] = project.name, now, invdata
                updated = True
                break

    if not failures:
        pass
    elif len(failures) < len(project.locations):
        LOGGER.info(
            __(
                'encountered some issues with some of the inventories,'
                ' but they had working alternatives:'
            )
        )
        for fail in failures:
            LOGGER.info(*fail)
    else:
        issues = '\n'.join(f[0] % f[1:] for f in failures)
        LOGGER.warning(
            __('failed to reach any of the inventories ' 'with the following issues:')
            + '\n'
            + issues
        )
    return updated


def fetch_inventory(app: Sphinx, uri: InventoryURI, inv: str) -> Inventory:
    """Fetch, parse and return an intersphinx inventory file."""
    return _fetch_inventory(
        target_uri=uri,
        inv_location=inv,
        config=app.config,
        srcdir=app.srcdir,
    )


def _fetch_inventory(
    *, target_uri: InventoryURI, inv_location: str, config: Config, srcdir: Path
) -> Inventory:
    """Fetch, parse and return an intersphinx inventory file."""
    # both *target_uri* (base URI of the links to generate)
    # and *inv_location* (actual location of the inventory file)
    # can be local or remote URIs
    if '://' in target_uri:
        # case: inv URI points to remote resource; strip any existing auth
        target_uri = _strip_basic_auth(target_uri)
    try:
        if '://' in inv_location:
            f: _ReadableStream[bytes] = _read_from_url(inv_location, config=config)
        else:
            f = open(path.join(srcdir, inv_location), 'rb')  # NoQA: SIM115
    except Exception as err:
        err.args = (
            'intersphinx inventory %r not fetchable due to %s: %s',
            inv_location,
            err.__class__,
            str(err),
        )
        raise
    try:
        if hasattr(f, 'url'):
            new_inv_location = f.url
            if inv_location != new_inv_location:
                msg = __('intersphinx inventory has moved: %s -> %s')
                LOGGER.info(msg, inv_location, new_inv_location)

                if target_uri in {
                    inv_location,
                    path.dirname(inv_location),
                    path.dirname(inv_location) + '/',
                }:
                    target_uri = path.dirname(new_inv_location)
        with f:
            try:
                invdata = InventoryFile.load(f, target_uri, posixpath.join)
            except ValueError as exc:
                msg = f'unknown or unsupported inventory version: {exc!r}'
                raise ValueError(msg) from exc
    except Exception as err:
        err.args = (
            'intersphinx inventory %r not readable due to %s: %s',
            inv_location,
            err.__class__.__name__,
            str(err),
        )
        raise
    else:
        return invdata


def _get_safe_url(url: str) -> str:
    """Gets version of *url* with basic auth passwords obscured. This function
    returns results suitable for printing and logging.

    E.g.: https://user:12345@example.com => https://user@example.com

    :param url: a url
    :type url: ``str``

    :return: *url* with password removed
    :rtype: ``str``
    """
    parts = urlsplit(url)
    if parts.username is None:
        return url
    else:
        frags = list(parts)
        if parts.port:
            frags[1] = f'{parts.username}@{parts.hostname}:{parts.port}'
        else:
            frags[1] = f'{parts.username}@{parts.hostname}'

        return urlunsplit(frags)


def _strip_basic_auth(url: str) -> str:
    """Returns *url* with basic auth credentials removed. Also returns the
    basic auth username and password if they're present in *url*.

    E.g.: https://user:pass@example.com => https://example.com

    *url* need not include basic auth credentials.

    :param url: url which may or may not contain basic auth credentials
    :type url: ``str``

    :return: *url* with any basic auth creds removed
    :rtype: ``str``
    """
    frags = list(urlsplit(url))
    # swap out 'user[:pass]@hostname' for 'hostname'
    if '@' in frags[1]:
        frags[1] = frags[1].split('@')[1]
    return urlunsplit(frags)


def _read_from_url(url: str, *, config: Config) -> HTTPResponse:
    """Reads data from *url* with an HTTP *GET*.

    This function supports fetching from resources which use basic HTTP auth as
    laid out by RFC1738 § 3.1. See § 5 for grammar definitions for URLs.

    .. seealso:

       https://www.ietf.org/rfc/rfc1738.txt

    :param url: URL of an HTTP resource
    :type url: ``str``

    :return: data read from resource described by *url*
    :rtype: ``file``-like object
    """
    r = requests.get(
        url,
        stream=True,
        timeout=config.intersphinx_timeout,
        _user_agent=config.user_agent,
        _tls_info=(config.tls_verify, config.tls_cacerts),
    )
    r.raise_for_status()

    # For inv_location / new_inv_location
    r.raw.url = r.url  # type: ignore[union-attr]

    # Decode content-body based on the header.
    # xref: https://github.com/psf/requests/issues/2155
    r.raw.decode_content = True
    return r.raw

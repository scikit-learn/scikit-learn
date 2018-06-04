"""Dependency Resolution

The dependency resolution in pip is performed as follows:

for top-level requirements:
    a. only one spec allowed per project, regardless of conflicts or not.
       otherwise a "double requirement" exception is raised
    b. they override sub-dependency requirements.
for sub-dependencies
    a. "first found, wins" (where the order is breadth first)
"""

import logging
from collections import defaultdict
from itertools import chain

from pip._internal.exceptions import (
    BestVersionAlreadyInstalled, DistributionNotFound, HashError, HashErrors,
    UnsupportedPythonVersion,
)
from pip._internal.req.req_install import InstallRequirement
from pip._internal.utils.logging import indent_log
from pip._internal.utils.misc import dist_in_usersite, ensure_dir
from pip._internal.utils.packaging import check_dist_requires_python

logger = logging.getLogger(__name__)


class Resolver(object):
    """Resolves which packages need to be installed/uninstalled to perform \
    the requested operation without breaking the requirements of any package.
    """

    _allowed_strategies = {"eager", "only-if-needed", "to-satisfy-only"}

    def __init__(self, preparer, session, finder, wheel_cache, use_user_site,
                 ignore_dependencies, ignore_installed, ignore_requires_python,
                 force_reinstall, isolated, upgrade_strategy):
        super(Resolver, self).__init__()
        assert upgrade_strategy in self._allowed_strategies

        self.preparer = preparer
        self.finder = finder
        self.session = session

        # NOTE: This would eventually be replaced with a cache that can give
        #       information about both sdist and wheels transparently.
        self.wheel_cache = wheel_cache

        self.require_hashes = None  # This is set in resolve

        self.upgrade_strategy = upgrade_strategy
        self.force_reinstall = force_reinstall
        self.isolated = isolated
        self.ignore_dependencies = ignore_dependencies
        self.ignore_installed = ignore_installed
        self.ignore_requires_python = ignore_requires_python
        self.use_user_site = use_user_site

        self._discovered_dependencies = defaultdict(list)

    def resolve(self, requirement_set):
        """Resolve what operations need to be done

        As a side-effect of this method, the packages (and their dependencies)
        are downloaded, unpacked and prepared for installation. This
        preparation is done by ``pip.operations.prepare``.

        Once PyPI has static dependency metadata available, it would be
        possible to move the preparation to become a step separated from
        dependency resolution.
        """
        # make the wheelhouse
        if self.preparer.wheel_download_dir:
            ensure_dir(self.preparer.wheel_download_dir)

        # If any top-level requirement has a hash specified, enter
        # hash-checking mode, which requires hashes from all.
        root_reqs = (
            requirement_set.unnamed_requirements +
            list(requirement_set.requirements.values())
        )
        self.require_hashes = (
            requirement_set.require_hashes or
            any(req.has_hash_options for req in root_reqs)
        )

        # Display where finder is looking for packages
        locations = self.finder.get_formatted_locations()
        if locations:
            logger.info(locations)

        # Actually prepare the files, and collect any exceptions. Most hash
        # exceptions cannot be checked ahead of time, because
        # req.populate_link() needs to be called before we can make decisions
        # based on link type.
        discovered_reqs = []
        hash_errors = HashErrors()
        for req in chain(root_reqs, discovered_reqs):
            try:
                discovered_reqs.extend(
                    self._resolve_one(requirement_set, req)
                )
            except HashError as exc:
                exc.req = req
                hash_errors.append(exc)

        if hash_errors:
            raise hash_errors

    def _is_upgrade_allowed(self, req):
        if self.upgrade_strategy == "to-satisfy-only":
            return False
        elif self.upgrade_strategy == "eager":
            return True
        else:
            assert self.upgrade_strategy == "only-if-needed"
            return req.is_direct

    def _set_req_to_reinstall(self, req):
        """
        Set a requirement to be installed.
        """
        # Don't uninstall the conflict if doing a user install and the
        # conflict is not a user install.
        if not self.use_user_site or dist_in_usersite(req.satisfied_by):
            req.conflicts_with = req.satisfied_by
        req.satisfied_by = None

    # XXX: Stop passing requirement_set for options
    def _check_skip_installed(self, req_to_install):
        """Check if req_to_install should be skipped.

        This will check if the req is installed, and whether we should upgrade
        or reinstall it, taking into account all the relevant user options.

        After calling this req_to_install will only have satisfied_by set to
        None if the req_to_install is to be upgraded/reinstalled etc. Any
        other value will be a dist recording the current thing installed that
        satisfies the requirement.

        Note that for vcs urls and the like we can't assess skipping in this
        routine - we simply identify that we need to pull the thing down,
        then later on it is pulled down and introspected to assess upgrade/
        reinstalls etc.

        :return: A text reason for why it was skipped, or None.
        """
        if self.ignore_installed:
            return None

        req_to_install.check_if_exists(self.use_user_site)
        if not req_to_install.satisfied_by:
            return None

        if self.force_reinstall:
            self._set_req_to_reinstall(req_to_install)
            return None

        if not self._is_upgrade_allowed(req_to_install):
            if self.upgrade_strategy == "only-if-needed":
                return 'already satisfied, skipping upgrade'
            return 'already satisfied'

        # Check for the possibility of an upgrade.  For link-based
        # requirements we have to pull the tree down and inspect to assess
        # the version #, so it's handled way down.
        if not req_to_install.link:
            try:
                self.finder.find_requirement(req_to_install, upgrade=True)
            except BestVersionAlreadyInstalled:
                # Then the best version is installed.
                return 'already up-to-date'
            except DistributionNotFound:
                # No distribution found, so we squash the error.  It will
                # be raised later when we re-try later to do the install.
                # Why don't we just raise here?
                pass

        self._set_req_to_reinstall(req_to_install)
        return None

    def _get_abstract_dist_for(self, req):
        """Takes a InstallRequirement and returns a single AbstractDist \
        representing a prepared variant of the same.
        """
        assert self.require_hashes is not None, (
            "require_hashes should have been set in Resolver.resolve()"
        )

        if req.editable:
            return self.preparer.prepare_editable_requirement(
                req, self.require_hashes, self.use_user_site, self.finder,
            )

        # satisfied_by is only evaluated by calling _check_skip_installed,
        # so it must be None here.
        assert req.satisfied_by is None
        skip_reason = self._check_skip_installed(req)

        if req.satisfied_by:
            return self.preparer.prepare_installed_requirement(
                req, self.require_hashes, skip_reason
            )

        upgrade_allowed = self._is_upgrade_allowed(req)
        abstract_dist = self.preparer.prepare_linked_requirement(
            req, self.session, self.finder, upgrade_allowed,
            self.require_hashes
        )

        # NOTE
        # The following portion is for determining if a certain package is
        # going to be re-installed/upgraded or not and reporting to the user.
        # This should probably get cleaned up in a future refactor.

        # req.req is only avail after unpack for URL
        # pkgs repeat check_if_exists to uninstall-on-upgrade
        # (#14)
        if not self.ignore_installed:
            req.check_if_exists(self.use_user_site)

        if req.satisfied_by:
            should_modify = (
                self.upgrade_strategy != "to-satisfy-only" or
                self.force_reinstall or
                self.ignore_installed or
                req.link.scheme == 'file'
            )
            if should_modify:
                self._set_req_to_reinstall(req)
            else:
                logger.info(
                    'Requirement already satisfied (use --upgrade to upgrade):'
                    ' %s', req,
                )

        return abstract_dist

    def _resolve_one(self, requirement_set, req_to_install):
        """Prepare a single requirements file.

        :return: A list of additional InstallRequirements to also install.
        """
        # Tell user what we are doing for this requirement:
        # obtain (editable), skipping, processing (local url), collecting
        # (remote url or package name)
        if req_to_install.constraint or req_to_install.prepared:
            return []

        req_to_install.prepared = True

        # register tmp src for cleanup in case something goes wrong
        requirement_set.reqs_to_cleanup.append(req_to_install)

        abstract_dist = self._get_abstract_dist_for(req_to_install)

        # Parse and return dependencies
        dist = abstract_dist.dist(self.finder)
        try:
            check_dist_requires_python(dist)
        except UnsupportedPythonVersion as err:
            if self.ignore_requires_python:
                logger.warning(err.args[0])
            else:
                raise

        more_reqs = []

        def add_req(subreq, extras_requested):
            sub_install_req = InstallRequirement.from_req(
                str(subreq),
                req_to_install,
                isolated=self.isolated,
                wheel_cache=self.wheel_cache,
            )
            parent_req_name = req_to_install.name
            to_scan_again, add_to_parent = requirement_set.add_requirement(
                sub_install_req,
                parent_req_name=parent_req_name,
                extras_requested=extras_requested,
            )
            if parent_req_name and add_to_parent:
                self._discovered_dependencies[parent_req_name].append(
                    add_to_parent
                )
            more_reqs.extend(to_scan_again)

        with indent_log():
            # We add req_to_install before its dependencies, so that we
            # can refer to it when adding dependencies.
            if not requirement_set.has_requirement(req_to_install.name):
                # 'unnamed' requirements will get added here
                req_to_install.is_direct = True
                requirement_set.add_requirement(
                    req_to_install, parent_req_name=None,
                )

            if not self.ignore_dependencies:
                if req_to_install.extras:
                    logger.debug(
                        "Installing extra requirements: %r",
                        ','.join(req_to_install.extras),
                    )
                missing_requested = sorted(
                    set(req_to_install.extras) - set(dist.extras)
                )
                for missing in missing_requested:
                    logger.warning(
                        '%s does not provide the extra \'%s\'',
                        dist, missing
                    )

                available_requested = sorted(
                    set(dist.extras) & set(req_to_install.extras)
                )
                for subreq in dist.requires(available_requested):
                    add_req(subreq, extras_requested=available_requested)

            if not req_to_install.editable and not req_to_install.satisfied_by:
                # XXX: --no-install leads this to report 'Successfully
                # downloaded' for only non-editable reqs, even though we took
                # action on them.
                requirement_set.successfully_downloaded.append(req_to_install)

        return more_reqs

    def get_installation_order(self, req_set):
        """Create the installation order.

        The installation order is topological - requirements are installed
        before the requiring thing. We break cycles at an arbitrary point,
        and make no other guarantees.
        """
        # The current implementation, which we may change at any point
        # installs the user specified things in the order given, except when
        # dependencies must come earlier to achieve topological order.
        order = []
        ordered_reqs = set()

        def schedule(req):
            if req.satisfied_by or req in ordered_reqs:
                return
            if req.constraint:
                return
            ordered_reqs.add(req)
            for dep in self._discovered_dependencies[req.name]:
                schedule(dep)
            order.append(req)

        for install_req in req_set.requirements.values():
            schedule(install_req)
        return order

from __future__ import absolute_import

import logging
from collections import OrderedDict

from pip._internal.exceptions import InstallationError
from pip._internal.utils.logging import indent_log
from pip._internal.wheel import Wheel

logger = logging.getLogger(__name__)


class RequirementSet(object):

    def __init__(self, require_hashes=False):
        """Create a RequirementSet.
        """

        self.requirements = OrderedDict()
        self.require_hashes = require_hashes

        # Mapping of alias: real_name
        self.requirement_aliases = {}
        self.unnamed_requirements = []
        self.successfully_downloaded = []
        self.reqs_to_cleanup = []

    def __str__(self):
        reqs = [req for req in self.requirements.values()
                if not req.comes_from]
        reqs.sort(key=lambda req: req.name.lower())
        return ' '.join([str(req.req) for req in reqs])

    def __repr__(self):
        reqs = [req for req in self.requirements.values()]
        reqs.sort(key=lambda req: req.name.lower())
        reqs_str = ', '.join([str(req.req) for req in reqs])
        return ('<%s object; %d requirement(s): %s>'
                % (self.__class__.__name__, len(reqs), reqs_str))

    def add_requirement(self, install_req, parent_req_name=None,
                        extras_requested=None):
        """Add install_req as a requirement to install.

        :param parent_req_name: The name of the requirement that needed this
            added. The name is used because when multiple unnamed requirements
            resolve to the same name, we could otherwise end up with dependency
            links that point outside the Requirements set. parent_req must
            already be added. Note that None implies that this is a user
            supplied requirement, vs an inferred one.
        :param extras_requested: an iterable of extras used to evaluate the
            environment markers.
        :return: Additional requirements to scan. That is either [] if
            the requirement is not applicable, or [install_req] if the
            requirement is applicable and has just been added.
        """
        name = install_req.name
        if not install_req.match_markers(extras_requested):
            logger.info("Ignoring %s: markers '%s' don't match your "
                        "environment", install_req.name,
                        install_req.markers)
            return [], None

        # This check has to come after we filter requirements with the
        # environment markers.
        if install_req.link and install_req.link.is_wheel:
            wheel = Wheel(install_req.link.filename)
            if not wheel.supported():
                raise InstallationError(
                    "%s is not a supported wheel on this platform." %
                    wheel.filename
                )

        # This next bit is really a sanity check.
        assert install_req.is_direct == (parent_req_name is None), (
            "a direct req shouldn't have a parent and also, "
            "a non direct req should have a parent"
        )

        if not name:
            # url or path requirement w/o an egg fragment
            self.unnamed_requirements.append(install_req)
            return [install_req], None
        else:
            try:
                existing_req = self.get_requirement(name)
            except KeyError:
                existing_req = None
            if (parent_req_name is None and existing_req and not
                    existing_req.constraint and
                    existing_req.extras == install_req.extras and not
                    existing_req.req.specifier == install_req.req.specifier):
                raise InstallationError(
                    'Double requirement given: %s (already in %s, name=%r)'
                    % (install_req, existing_req, name))
            if not existing_req:
                # Add requirement
                self.requirements[name] = install_req
                # FIXME: what about other normalizations?  E.g., _ vs. -?
                if name.lower() != name:
                    self.requirement_aliases[name.lower()] = name
                result = [install_req]
            else:
                # Assume there's no need to scan, and that we've already
                # encountered this for scanning.
                result = []
                if not install_req.constraint and existing_req.constraint:
                    if (install_req.link and not (existing_req.link and
                       install_req.link.path == existing_req.link.path)):
                        self.reqs_to_cleanup.append(install_req)
                        raise InstallationError(
                            "Could not satisfy constraints for '%s': "
                            "installation from path or url cannot be "
                            "constrained to a version" % name,
                        )
                    # If we're now installing a constraint, mark the existing
                    # object for real installation.
                    existing_req.constraint = False
                    existing_req.extras = tuple(
                        sorted(set(existing_req.extras).union(
                               set(install_req.extras))))
                    logger.debug("Setting %s extras to: %s",
                                 existing_req, existing_req.extras)
                    # And now we need to scan this.
                    result = [existing_req]
                # Canonicalise to the already-added object for the backref
                # check below.
                install_req = existing_req

            # We return install_req here to allow for the caller to add it to
            # the dependency information for the parent package.
            return result, install_req

    def has_requirement(self, project_name):
        name = project_name.lower()
        if (name in self.requirements and
           not self.requirements[name].constraint or
           name in self.requirement_aliases and
           not self.requirements[self.requirement_aliases[name]].constraint):
            return True
        return False

    @property
    def has_requirements(self):
        return list(req for req in self.requirements.values() if not
                    req.constraint) or self.unnamed_requirements

    def get_requirement(self, project_name):
        for name in project_name, project_name.lower():
            if name in self.requirements:
                return self.requirements[name]
            if name in self.requirement_aliases:
                return self.requirements[self.requirement_aliases[name]]
        raise KeyError("No project with the name %r" % project_name)

    def cleanup_files(self):
        """Clean up files, remove builds."""
        logger.debug('Cleaning up...')
        with indent_log():
            for req in self.reqs_to_cleanup:
                req.remove_temporary_source()

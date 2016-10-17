"""Check whether we or not we should build the documentation

If the last commit message has a "[doc skip]" marker, do not build
the doc. On the contrary if a "[doc build]" marker is found, build the doc
instead of relying on the subsequent rules.

We always build the documentation for jobs that are not related to a specific
PR (e.g. a merge to master or a maintenance branch).

If this is a PR, check that if there are some files in this PR that are under
the "doc/" or "examples/" folders, otherwise skip.

If the introspection of the current commit fails for any reason, the default
behavior is to build the documentation.

"""
import sys
import os
from subprocess import check_output, CalledProcessError


def exit(msg="", skip=False):
    print("%s: %s" % ("SKIP" if skip else "BUILD", msg))
    sys.exit(0)

# Introspect the message for the commit that triggered the build
commit = os.environ.get('CIRCLE_SHA1')
if not commit:
    exit("undefined CIRCLE_SHA1 variable")
try:
    commit_msg = check_output("git log --format=%B -n 1".split() + [commit])
    commit_msg = commit_msg.decode('utf-8')
except CalledProcessError:
    exit("failed to introspect commit message for %s" % commit)

if "[doc skip]" in commit_msg:
    exit("[doc skip] marker found", skip=True)
elif "[doc build]" in commit_msg:
    exit("[doc build] marker found")

# Check whether this commit is part of a pull request or not
pr_url = os.environ.get('CI_PULL_REQUEST')
if not pr_url:
    # The documentation should be always built when executed from one of the
    # main branches
    exit("not a pull request")

# Introspect the list of files changed by all the commits in this PR.
# Hardcode the assumption that this is a PR to origin/master of this repo
# as apparently there is way to reliably get the target of a PR with circle
# ci
git_range = "origin/master...%s" % commit
try:
    check_output("git fetch origin master".split())
    filenames = check_output("git diff --name-only".split() + [git_range])
except CalledProcessError:
    exit("git introspection failed.")
filenames = filenames.decode('utf-8').split()
for filename in filenames:
    if filename.startswith(u'doc/') or filename.startswith(u'examples/'):
        exit("detected doc impacting file modified by PR in range %s: %s"
             % (git_range, filename))

# This PR does not seem to have any documentation related file changed.
msg = "no doc impacting files detected:\n" + u"\n".join(filenames)
exit(msg, skip=True)

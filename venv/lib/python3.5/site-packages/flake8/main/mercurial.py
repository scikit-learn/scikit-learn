"""Module containing the main mecurial hook interface and helpers.

.. autofunction:: hook
.. autofunction:: install

"""
import configparser
import os
import subprocess

from flake8 import exceptions as exc

__all__ = ('hook', 'install')


def hook(ui, repo, **kwargs):
    """Execute Flake8 on the repository provided by Mercurial.

    To understand the parameters read more of the Mercurial documentation
    around Hooks: https://www.mercurial-scm.org/wiki/Hook.

    We avoid using the ``ui`` attribute because it can cause issues with
    the GPL license tha Mercurial is under. We don't import it, but we
    avoid using it all the same.
    """
    from flake8.main import application
    hgrc = find_hgrc(create_if_missing=False)
    if hgrc is None:
        print('Cannot locate your root mercurial repository.')
        raise SystemExit(True)

    hgconfig = configparser_for(hgrc)
    strict = hgconfig.get('flake8', 'strict', fallback=True)

    filenames = list(get_filenames_from(repo, kwargs))

    app = application.Application()
    app.initialize(filenames)
    app.options._running_from_vcs = True
    app.run_checks()
    app.report()

    if strict:
        return app.result_count
    return 0


def install():
    """Ensure that the mercurial hooks are installed.

    This searches for the ``.hg/hgrc`` configuration file and will add commit
    and qrefresh hooks to it, if they do not already exist.

    It will also print a message to stdout about how to configure the hook.

    :returns:
        True if successful, False if the ``.hg/hgrc`` file doesn't exist.
    :rtype:
        bool
    :raises:
        flake8.exceptions.MercurialCommitHookAlreadyExists
    :raises:
        flake8.exceptions.MercurialQRefreshHookAlreadyExists
    """
    hgrc = find_hgrc(create_if_missing=True)
    if hgrc is None:
        return False

    hgconfig = configparser_for(hgrc)

    if not hgconfig.has_section('hooks'):
        hgconfig.add_section('hooks')

    if hgconfig.has_option('hooks', 'commit'):
        raise exc.MercurialCommitHookAlreadyExists(
            path=hgrc,
            value=hgconfig.get('hooks', 'commit'),
        )

    if hgconfig.has_option('hooks', 'qrefresh'):
        raise exc.MercurialQRefreshHookAlreadyExists(
            path=hgrc,
            value=hgconfig.get('hooks', 'qrefresh'),
        )

    hgconfig.set('hooks', 'commit', 'python:flake8.main.mercurial.hook')
    hgconfig.set('hooks', 'qrefresh', 'python:flake8.main.mercurial.hook')

    if not hgconfig.has_section('flake8'):
        hgconfig.add_section('flake8')

    if not hgconfig.has_option('flake8', 'strict'):
        hgconfig.set('flake8', 'strict', False)

    with open(hgrc, 'w') as fd:
        hgconfig.write(fd)

    print('mercurial hooks installed, for configuration options see')
    print('http://flake8.pycqa.org/en/latest/user/using-hooks.html')

    return True


def get_filenames_from(repository, kwargs):
    seen_filenames = set()
    node = kwargs['node']
    for revision in range(repository[node], len(repository)):
        for filename in repository[revision].files():
            full_filename = os.path.join(repository.root, filename)
            have_seen_filename = full_filename in seen_filenames
            filename_does_not_exist = not os.path.exists(full_filename)
            if have_seen_filename or filename_does_not_exist:
                continue

            seen_filenames.add(full_filename)
            if full_filename.endswith('.py'):
                yield full_filename


def find_hgrc(create_if_missing=False):
    root = subprocess.Popen(
        ['hg', 'root'],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    (hg_directory, _) = root.communicate()
    if callable(getattr(hg_directory, 'decode', None)):
        hg_directory = hg_directory.decode('utf-8')

    if not os.path.isdir(hg_directory):
        return None

    hgrc = os.path.abspath(
        os.path.join(hg_directory, '.hg', 'hgrc')
    )
    if not os.path.exists(hgrc):
        if create_if_missing:
            open(hgrc, 'w').close()
        else:
            return None

    return hgrc


def configparser_for(path):
    parser = configparser.ConfigParser(interpolation=None)
    parser.read(path)
    return parser

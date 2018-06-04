"""
Password generation for the IPython notebook.
"""
#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
# Stdlib
import getpass
import hashlib
import random

# Our own
from IPython.core.error import UsageError
from IPython.utils.py3compat import encode

#-----------------------------------------------------------------------------
# Globals
#-----------------------------------------------------------------------------

# Length of the salt in nr of hex chars, which implies salt_len * 4
# bits of randomness.
salt_len = 12

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------

def passwd(passphrase=None, algorithm='sha1'):
    """Generate hashed password and salt for use in notebook configuration.

    In the notebook configuration, set `c.NotebookApp.password` to
    the generated string.

    Parameters
    ----------
    passphrase : str
        Password to hash.  If unspecified, the user is asked to input
        and verify a password.
    algorithm : str
        Hashing algorithm to use (e.g, 'sha1' or any argument supported
        by :func:`hashlib.new`).

    Returns
    -------
    hashed_passphrase : str
        Hashed password, in the format 'hash_algorithm:salt:passphrase_hash'.

    Examples
    --------
    >>> passwd('mypassword')
    'sha1:7cf3:b7d6da294ea9592a9480c8f52e63cd42cfb9dd12'

    """
    if passphrase is None:
        for i in range(3):
            p0 = getpass.getpass('Enter password: ')
            p1 = getpass.getpass('Verify password: ')
            if p0 == p1:
                passphrase = p0
                break
            else:
                print('Passwords do not match.')
        else:
            raise UsageError('No matching passwords found. Giving up.')

    h = hashlib.new(algorithm)
    salt = ('%0' + str(salt_len) + 'x') % random.getrandbits(4 * salt_len)
    h.update(encode(passphrase, 'utf-8') + encode(salt, 'ascii'))

    return ':'.join((algorithm, salt, h.hexdigest()))


def passwd_check(hashed_passphrase, passphrase):
    """Verify that a given passphrase matches its hashed version.

    Parameters
    ----------
    hashed_passphrase : str
        Hashed password, in the format returned by `passwd`.
    passphrase : str
        Passphrase to validate.

    Returns
    -------
    valid : bool
        True if the passphrase matches the hash.

    Examples
    --------
    >>> from IPython.lib.security import passwd_check
    >>> passwd_check('sha1:0e112c3ddfce:a68df677475c2b47b6e86d0467eec97ac5f4b85a',
    ...              'mypassword')
    True

    >>> passwd_check('sha1:0e112c3ddfce:a68df677475c2b47b6e86d0467eec97ac5f4b85a',
    ...              'anotherpassword')
    False
    """
    try:
        algorithm, salt, pw_digest = hashed_passphrase.split(':', 2)
    except (ValueError, TypeError):
        return False

    try:
        h = hashlib.new(algorithm)
    except ValueError:
        return False

    if len(pw_digest) == 0:
        return False

    h.update(encode(passphrase, 'utf-8') + encode(salt, 'ascii'))

    return h.hexdigest() == pw_digest

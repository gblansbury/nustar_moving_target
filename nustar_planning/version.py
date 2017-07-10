# Autogenerated by Astropy-affiliated package nustar_planning's setup.py on 2017-07-10 12:35:55.001534
from __future__ import unicode_literals
import datetime


import locale
import os
import subprocess
import warnings


def _decode_stdio(stream):
    try:
        stdio_encoding = locale.getdefaultlocale()[1] or 'utf-8'
    except ValueError:
        stdio_encoding = 'utf-8'

    try:
        text = stream.decode(stdio_encoding)
    except UnicodeDecodeError:
        # Final fallback
        text = stream.decode('latin1')

    return text


def update_git_devstr(version, path=None):
    """
    Updates the git revision string if and only if the path is being imported
    directly from a git working copy.  This ensures that the revision number in
    the version string is accurate.
    """

    try:
        # Quick way to determine if we're in git or not - returns '' if not
        devstr = get_git_devstr(sha=True, show_warning=False, path=path)
    except OSError:
        return version

    if not devstr:
        # Probably not in git so just pass silently
        return version

    if 'dev' in version:  # update to the current git revision
        version_base = version.split('.dev', 1)[0]
        devstr = get_git_devstr(sha=False, show_warning=False, path=path)

        return version_base + '.dev' + devstr
    else:
        # otherwise it's already the true/release version
        return version


def get_git_devstr(sha=False, show_warning=True, path=None):
    """
    Determines the number of revisions in this repository.

    Parameters
    ----------
    sha : bool
        If True, the full SHA1 hash will be returned. Otherwise, the total
        count of commits in the repository will be used as a "revision
        number".

    show_warning : bool
        If True, issue a warning if git returns an error code, otherwise errors
        pass silently.

    path : str or None
        If a string, specifies the directory to look in to find the git
        repository.  If `None`, the current working directory is used, and must
        be the root of the git repository.
        If given a filename it uses the directory containing that file.

    Returns
    -------
    devversion : str
        Either a string with the revision number (if `sha` is False), the
        SHA1 hash of the current commit (if `sha` is True), or an empty string
        if git version info could not be identified.

    """

    if path is None:
        path = os.getcwd()
        if not _get_repo_path(path, levels=0):
            return ''

    if not os.path.isdir(path):
        path = os.path.abspath(os.path.dirname(path))

    if sha:
        # Faster for getting just the hash of HEAD
        cmd = ['rev-parse', 'HEAD']
    else:
        cmd = ['rev-list', '--count', 'HEAD']

    def run_git(cmd):
        try:
            p = subprocess.Popen(['git'] + cmd, cwd=path,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 stdin=subprocess.PIPE)
            stdout, stderr = p.communicate()
        except OSError as e:
            if show_warning:
                warnings.warn('Error running git: ' + str(e))
            return (None, b'', b'')

        if p.returncode == 128:
            if show_warning:
                warnings.warn('No git repository present at {0!r}! Using '
                              'default dev version.'.format(path))
            return (p.returncode, b'', b'')
        if p.returncode == 129:
            if show_warning:
                warnings.warn('Your git looks old (does it support {0}?); '
                              'consider upgrading to v1.7.2 or '
                              'later.'.format(cmd[0]))
            return (p.returncode, stdout, stderr)
        elif p.returncode != 0:
            if show_warning:
                warnings.warn('Git failed while determining revision '
                              'count: {0}'.format(_decode_stdio(stderr)))
            return (p.returncode, stdout, stderr)

        return p.returncode, stdout, stderr

    returncode, stdout, stderr = run_git(cmd)

    if not sha and returncode == 129:
        # git returns 129 if a command option failed to parse; in
        # particular this could happen in git versions older than 1.7.2
        # where the --count option is not supported
        # Also use --abbrev-commit and --abbrev=0 to display the minimum
        # number of characters needed per-commit (rather than the full hash)
        cmd = ['rev-list', '--abbrev-commit', '--abbrev=0', 'HEAD']
        returncode, stdout, stderr = run_git(cmd)
        # Fall back on the old method of getting all revisions and counting
        # the lines
        if returncode == 0:
            return str(stdout.count(b'\n'))
        else:
            return ''
    elif sha:
        return _decode_stdio(stdout)[:40]
    else:
        return _decode_stdio(stdout).strip()


def _get_repo_path(pathname, levels=None):
    """
    Given a file or directory name, determine the root of the git repository
    this path is under.  If given, this won't look any higher than ``levels``
    (that is, if ``levels=0`` then the given path must be the root of the git
    repository and is returned if so.

    Returns `None` if the given path could not be determined to belong to a git
    repo.
    """

    if os.path.isfile(pathname):
        current_dir = os.path.abspath(os.path.dirname(pathname))
    elif os.path.isdir(pathname):
        current_dir = os.path.abspath(pathname)
    else:
        return None

    current_level = 0

    while levels is None or current_level <= levels:
        if os.path.exists(os.path.join(current_dir, '.git')):
            return current_dir

        current_level += 1
        if current_dir == os.path.dirname(current_dir):
            break

        current_dir = os.path.dirname(current_dir)

    return None

_packagename = "nustar_planning"
_last_generated_version = "0.1.dev38"
_last_githash = "bf3a0ea963fef64f7cb05334aabf9412b1f7e529"

# Determine where the source code for this module
# lives.  If __file__ is not a filesystem path then
# it is assumed not to live in a git repo at all.
if _get_repo_path(__file__, levels=len(_packagename.split('.'))):
    version = update_git_devstr(_last_generated_version, path=__file__)
    githash = get_git_devstr(sha=True, show_warning=False,
                             path=__file__) or _last_githash
else:
    # The file does not appear to live in a git repo so don't bother
    # invoking git
    version = _last_generated_version
    githash = _last_githash


major = 0
minor = 1
bugfix = 0

release = False
timestamp = datetime.datetime(2017, 7, 10, 12, 35, 55, 1534)
debug = False

try:
    from ._compiler import compiler
except ImportError:
    compiler = "unknown"

try:
    from .cython_version import cython_version
except ImportError:
    cython_version = "unknown"

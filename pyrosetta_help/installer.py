__doc__ = """
This module is aimed at helping to install pyrosetta.
If PyRosetta is installed but if it segfaults on import due to wrong GNU lib C (``glic``)
or there's a missing dependency that is a different matter.
"""

import argparse
import importlib
import importlib.util
import os
import re
import requests
import site
import sys
from typing import (Optional, Union)


def check_pyrosetta() -> bool:
    """
    this will return none if it is not installed regardless if it segfaults or there's a missing dependency.
    """
    return importlib.util.find_spec('pyrosetta') is not None


def install_pyrosetta(username: Optional[str] = None,
                      password: Optional[str] = None,
                      path: str = '.',
                      hash_comparison_required: bool = True):
    """
    If there's a folder in ``path`` with PyRosetta release it will install that.

    If username and password arent provided
    the evironmental variables ``PYROSETTA_USERNAME`` or ``PYROSETTA_PASSWORD`` are used.
    It checks to see if hashes of the password is correct or the Rosetta one.
    """
    # ## local version present?
    path = _get_release_path(path)
    if path:
        assert not os.system(f'yes | pip3 install -e {path}/setup/')
        site.main()
        return
    # ## install by download
    username = _parse_environmental(username, 'PYROSETTA_USERNAME')
    password = _parse_environmental(password, 'PYROSETTA_PASSWORD')
    # check if hash is not the vanilla rosetta
    _check_not_rosetta(username, password)
    if hash_comparison_required:
        _check_correct(username, password)
    # get latest release
    latest_url = get_latest_release_url(username, password)
    assert not os.system(f'yes | pip3 install {latest_url}')
    site.main()  # refresh
    return


def parse():
    if check_pyrosetta():
        print('PyRosetta is already installed.')
        SystemExit(0)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-u', '--username',
                        type=str, help='PyRosetta username')
    parser.add_argument('-p', '--password',
                        type=str, help='PyRosetta password')
    parser.add_argument('-f', '--file',
                        type=str,
                        help='Path where a release was release')
    args = parser.parse_args()
    install_pyrosetta(username=args.username,
                      password=args.password,
                      path=args.file,
                      hash_comparison_required=True)


# ===================================================================================================================


def get_os_name() -> str:
    import platform
    if platform.system() == 'Windows':
        raise SystemError('Please Google what one does on Windows or install a real OS')
    elif platform.system() == 'Darwin' and platform.processor() in ('i386', 'x86_64'):
        return 'mac'
    elif platform.system() == 'Darwin' and platform.processor() == 'arm':
        return 'm1'
    elif platform.system() == 'Linux' and platform.processor() in ('i386', 'x86_64'):
        return 'ubuntu'
        # I could read /etc/os-release but I know pyrosetta works on CentOS.
        # Humbug is actually CentOS.
    elif platform.system() == 'Linux' and platform.machine() == 'armv7l':
        raise SystemError('You need to compile PyRosetta yourself for a Raspberry Pi')
    else:
        raise SystemError('Cannot detect OS type')


def get_latest_release_url(username: str, password: str) -> str:
    """
    Doing ``pip install xx:xxx@https://xx:xxx@graylab.jhu.edu/download/PyRos.../latest.html``
    results in authetical loop drama in following the 302.
    """
    # assumes the system is Ubuntu
    py_version = str(sys.version_info.major) + str(sys.version_info.minor)
    machine = get_os_name()
    os_specific = f'PyRosetta4.Release.python{py_version}.{machine}.wheel'
    base_url = f'https://{username}:{password}@graylab.jhu.edu/download/PyRosetta4/archive/release/{os_specific}'
    url_to_latest = f'{base_url}/latest.html'
    latest_response = requests.get(url_to_latest,
                                   # auth=requests.auth.HTTPBasicAuth(username, password)
                                   )
    if latest_response.status_code == 401:
        raise ValueError('Incorrect username or password!')
    elif latest_response.status_code not in (200, 300, 301, 302, 303, 304, 305, 306, 307, 308):
        from IPython.display import display, HTML
        display(HTML(latest_response.text))
        raise ValueError(f'Something is wrong with the url {url_to_latest}')
    return base_url + '/' + re.search(r'[uU][rR][lL]=(.*?)["\s]', latest_response.text).group(1)


# ------------------------------------------------------------------------------------------------------------------

def _parse_environmental(provided, env_key) -> str:
    """
    return provided or env_key from os.environ
    PYROSETTA_USERNAME or PYROSETTA_PASSWORD
    """
    if isinstance(provided, str):
        return provided.strip()
    if provided is not None:
        raise TypeError(f'{env_key} of type {type(provided).__name__}')
    elif provided is None and env_key in os.environ:
        return os.environ[env_key].strip()
    else:
        raise ValueError(f'no {env_key} specified')


def _sha256_hash(value: str) -> str:
    """You cannot unhash a hash. Do not both."""
    import hashlib
    return hashlib.sha256(value.encode()).hexdigest()


def _check_not_rosetta(username: str, password: str) -> None:
    """
    A common issue is users trying to install PyRosetta with Rosetta credentials.
    All operations are hash based. No plain text passwords: this aint talktalk!
    """
    msg = 'Rosetta and PyRosetta use different credentials. ' + \
          'You appear to have tried Rosetta not PyRosetta credentials'
    rosetta_username_sha256 = '26f5f1e313f549389c5159162ea27a50cf519d1701312f3e3429bdb7b6b6bb26'
    rosetta_password_sha256 = '0c49db1c6501f859e11c6546e20b55541245474d88935a9e9ffb5041ef82f5b3'
    assert _sha256_hash(username) != rosetta_username_sha256, msg
    assert _sha256_hash(password) != rosetta_password_sha256, msg


def _check_correct(username: str, password: str) -> None:
    """
    Verify the username and password are correct without actually knowing them.
    I worry that there may be a large number of idiots that try to guess the passwords
    thus getting Colab fail2ban jailed.
    All operations are hash based. No plain text passwords: this aint talktalk!
    """
    hashed_username = _sha256_hash(username)
    hashed_password = _sha256_hash(password)
    expected_hashed_username = 'cf6f296b8145262b22721e52e2edec13ce57af8c6fc990c8ae1a4aa3e50ae40e'
    expected_hashed_password = '45066dd976d8bf0c05dc8dd4d58727945c3437e6eb361ba9870097968db7a0da'
    msg = ('The hash of the {} is not as expected: ' +
           'if your username and password combo are correct and just not the academic ones, ' +
           'set `hash_comparision_required` to False.' +
           'If you dont know the password visit the Rosetta Commons site. ' +
           'The password and username for PyRosetta are DIFFERENT than Rosetta or FoldIt. ' +
           'Please do not google for username:password for this or for anything ' +
           'as stats like number of registered users is mighty important for grants and stuff.')
    assert hashed_username == expected_hashed_username, msg.format('username')
    assert hashed_password == expected_hashed_password, msg.format('password')


def _get_release_path(path) -> Union[None, str]:
    """get the folder with the given folder...
    The idiocy is due to the long name a downloaded release may have
    """
    if path is None:
        return None
    elif not os.path.isdir(path):
        return None
        # the path provided is a pyrosetta release itself
    elif 'PyRosetta4.Release' in path and 'setup' in os.listdir(path):
        return path
    for filename in os.listdir(path):
        subpath = _get_release_path(os.path.join(path, filename))
        if subpath:
            return subpath


if __name__ == '__main__':
    parse()

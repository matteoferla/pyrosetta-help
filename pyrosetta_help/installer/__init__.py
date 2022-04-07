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
from typing import (Optional)
from ._aux import *


def check_pyrosetta() -> bool:
    """
    this will return none if it is not installed regardless if it segfaults or there's a missing dependency.
    """
    return importlib.util.find_spec('pyrosetta') is not None


def download_pyrosetta(username: Optional[str] = None,
                       password: Optional[str] = None,
                       path: str = '.',
                       hash_comparison_required: bool = True,
                       on_preexisting:str='return') -> str:
    """
    ``path`` is wither it is saved.
    ``on_preexisting`` can be 'return', 'raise', 'delete' or 'ignore' (& download anyway).
    """
    existing_path = get_release_path(path)
    if existing_path is None:
        return existing_path
    elif on_preexisting == 'ignore':
        pass
    elif on_preexisting == 'delete':
        os.rmdir(existing_path)
    else: # error or raise
        raise FileExistsError(f'There seems to be a release already downloaded ({existing_path})')
    # get url and download:
    latest_url = get_latest_release_url(username, password,
                                        wheel=True,
                                        hash_comparison_required=hash_comparison_required)
    assert not os.system(f'curl {latest_url} | tar -xj -C {path}')
    pyrosetta_folder = re.sub(r'.tar.\w+$', '', os.path.split(latest_url)[1])
    return os.path.join(path, pyrosetta_folder)


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
    path = get_release_path(path)
    if path:
        assert not os.system(f'yes | pip3 install -e {path}/setup/')
        site.main()  # refresh
        return
    # ## install by download
    # get latest release
    latest_url = get_latest_release_url(username, password,
                                        wheel=True,
                                        hash_comparison_required=hash_comparison_required)
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


def get_latest_release_url(username: str,
                           password: str,
                           wheel: bool = True,
                           hash_comparison_required: bool = True) -> str:
    """
    Doing ``pip install xx:xxx@https://xx:xxx@graylab.jhu.edu/download/PyRos.../latest.html``
    results in authetical loop drama in following the 302.
    """

    username = parse_environmental(username, 'PYROSETTA_USERNAME')
    password = parse_environmental(password, 'PYROSETTA_PASSWORD')
    # check if hash is not the vanilla rosetta
    check_not_rosetta(username, password)
    if hash_comparison_required:
        check_correct(username, password)
    # get specifics:
    py_version = str(sys.version_info.major) + str(sys.version_info.minor)
    machine = get_os_name()
    # assemble
    os_specific = f'PyRosetta4.Release.python{py_version}.{machine}'
    if wheel:
        os_specific += '.wheel'
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


if __name__ == '__main__':
    parse()

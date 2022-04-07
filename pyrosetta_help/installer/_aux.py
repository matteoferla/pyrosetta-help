import os
from typing import (Union)


def parse_environmental(provided, env_key) -> str:
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


def sha256_hash(value: str) -> str:
    """You cannot unhash a hash. Do not both."""
    import hashlib
    return hashlib.sha256(value.encode()).hexdigest()


def check_not_rosetta(username: str, password: str) -> None:
    """
    A common issue is users trying to install PyRosetta with Rosetta credentials.
    All operations are hash based. No plain text passwords: this aint talktalk!
    """
    msg = 'Rosetta and PyRosetta use different credentials. ' + \
          'You appear to have tried Rosetta not PyRosetta credentials'
    rosetta_username_sha256 = '26f5f1e313f549389c5159162ea27a50cf519d1701312f3e3429bdb7b6b6bb26'
    rosetta_password_sha256 = '0c49db1c6501f859e11c6546e20b55541245474d88935a9e9ffb5041ef82f5b3'
    assert sha256_hash(username) != rosetta_username_sha256, msg
    assert sha256_hash(password) != rosetta_password_sha256, msg


def check_correct(username: str, password: str) -> None:
    """
    Verify the username and password are correct without actually knowing them.
    I worry that there may be a large number of idiots that try to guess the passwords
    thus getting Colab fail2ban jailed.
    All operations are hash based. No plain text passwords: this aint talktalk!
    """
    hashed_username = sha256_hash(username)
    hashed_password = sha256_hash(password)
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


def get_release_path(path) -> Union[None, str]:
    """get the folder with the given folder...
    The idiocy is due to the long name a downloaded release may have
    """
    if path is None:
        return None
    elif not os.path.isdir(path):
        return None
        # the path provided is a pyrosetta release itself
    elif 'pyrosetta' in path.lower() and 'setup' in os.listdir(path):
        return path
    for filename in os.listdir(path):
        subpath = get_release_path(os.path.join(path, filename))
        if subpath:
            return subpath

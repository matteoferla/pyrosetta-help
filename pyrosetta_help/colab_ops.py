import os
from typing import Optional
import importlib

__all__ = ['get_shell_mode', 'assert_notebook', 'mount_google_drive', 'install_and_import']

def get_shell_mode() -> str:
    """
    The ``get_ipython().__class__.__name__`` is very nasty, replacing it with
    """
    # get_ipython is a standard libary in ipython, but not script
    from IPython import get_ipython
    shell_name = get_ipython().__class__.__name__
    nicer_names = {'Shell': 'colab',
                   'ZMQInteractiveShell': 'jupyter',
                   'TerminalInteractiveShell': 'ipython',
                   None: 'terminal'}
    if shell_name in nicer_names:
        return nicer_names[shell_name]
    else:
        return shell_name


def assert_notebook():
    shell_name: str = get_shell_mode()
    if shell_name not in ('colab', 'jupyter'):
        raise RuntimeError(f'This is a colabs notebook. Why are you running it in {shell_name}?')


def mount_google_drive(use_drive: bool) -> str:
    """
    Mount drive if needed and return the working directory.
    Change the working path and
    """
    # ## Use drive?
    modality: str = get_shell_mode()
    if modality != 'colab':
        # jupyter --> stay
        return './'
    elif use_drive:
        from google.colab import drive  # noqa the module google cannot be installed by mere mortals AFAIK
        drive.mount('/content/drive')
        path = '/content/drive/MyDrive'
        os.chdir(path)
        return path
    else:  # --> stay
        return '/content'


def install_and_import(package_name: str,
                       pypi_name: Optional[str] = None,
                       alias_name: Optional[str] = None):
    """
    This is an import / installer in case the !pip commands do not work.
    Which they do unless there are issues with the kernel.

    If the module has a different name in pypi (`pypi_name`)
    than its import name (`package_name`), specify it.

         pip install pypi_name
         import package_name as alias_name
    """
    # I will go to hell for this, but shmeh:
    from pip._internal.cli.main import main as pip_main
    if pypi_name is None:
        pypi_name = package_name
    if alias_name is None:
        alias_name = package_name
    try:
        importlib.import_module(package_name)
    except ImportError as error:
        if error.name != package_name:
            # these are not the droids we are looking for
            raise ImportError(f'Import of {package_name} requires module {error.name}...',
                              name=error.name)
        pip_main(['install', pypi_name])
    globals()[alias_name] = importlib.import_module(package_name)

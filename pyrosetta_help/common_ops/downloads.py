__all__ = ['download_map',
           'download_cif',
           'download_pdb',
           'download_opm']

import requests
import shutil
import urllib.request as request
from contextlib import closing

def download_map(code: str):
    """
    This is two functions. EMD database for cryoEM files

    :param code:
    :return:
    """
    if isinstance(code, int) or "EMD-" in code or code.isdigit():
        return _download_mrc_map(str(code))
    else:
        return _download_ccp4_map(code)

def _download_mrc_map(code:str):
    # alt:
    # https://ftp.wwpdb.org/pub/emdb/structures/EMD-20808/map/emd_20808.map.gz
    code = code.replace("EMD-", "")
    ftp_path = f'ftp://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-{code}/map/emd_{code}.map.gz'
    file_path = f'EMD-{code}.map.gz'
    with closing(request.urlopen(ftp_path)) as r, open(file_path, 'wb') as f:
        shutil.copyfileobj(r, f)
    return file_path

def _download_ccp4_map(code:str):
    code = code.replace("pdb:", "")
    assert len(code) == 4, f'PDB code {code} is not 4 alphanumericals long'
    url = f'https://www.ebi.ac.uk/pdbe/coordinates/files/{code.lower()}.ccp4'
    file_path = f'{code}.ccp4'
    reply = requests.get(url, stream=True)
    assert reply.status_code == 200
    with open(file_path, 'wb') as f:
        reply.raw.decode_content = True
        shutil.copyfileobj(reply.raw, f)
    return file_path

def download_cif(code: str):
    """
    Download CIF. Pyrosetta has issues importing Cifs (hetatms).
    This is just a note to self as PyMOL fetch can do it.
    """
    http_path = f'http://www.ebi.ac.uk/pdbe/static/entry/download/{code.lower()}-assembly-1.cif.gz'
    file_path = f'{code}.cif.gz'
    r = requests.get(http_path, stream=True)
    if r.status_code == 200:
        with open(file_path, 'wb') as f:
            r.raw.decode_content = True
            shutil.copyfileobj(r.raw, f)
        return True
    else:
        return False

def download_pdb(pdbcode: str) -> str:
    filename = f'{pdbcode.lower()}.pdb'

    # r = requests.get(f'https://www.ebi.ac.uk/pdbe/coordinates/files/{pdbcode.lower()}.ccp4', stream=True)
    r = requests.get(f'https://www.ebi.ac.uk/pdbe/entry-files/download/pdb{pdbcode.lower()}.ent', stream=True)
    if r.status_code == 200:
        with open(filename, 'wb') as f:
            r.raw.decode_content = True
            shutil.copyfileobj(r.raw, f)
    return filename

def download_opm(code: str):
    """
    Download OPM PDB
    """
    import shutil, requests
    http_path = f'https://opm-assets.storage.googleapis.com/pdb/{code.lower()}.pdb'
    file_path = f'{code}_OMP.pdb'
    r = requests.get(http_path, stream=True)
    if r.status_code == 200:
        with open(file_path, 'wb') as f:
            r.raw.decode_content = True
            shutil.copyfileobj(r.raw, f)
        return True
    else:
        return False
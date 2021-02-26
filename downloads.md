## Download off PDBe

    import requests, shutil
    from functools import partial
    from typing import *
    
    def download_from_pdbe(code: str,
                           base_url: str,
                           suffix: str,
                           save_suffix: Optional[str]=None):
        http_path = f'{base_url}/{code.lower()}{suffix}'
        if save_suffix:
            file_path = f'{code}{save_suffix}'
        else:
            file_path = f'{code}{suffix}'
        r = requests.get(http_path, stream=True)
        print(http_path)
        if r.status_code == 200:
            with open(file_path, 'wb') as f:
                r.raw.decode_content = True
                shutil.copyfileobj(r.raw, f)
            return True
        else:
            return False
            
    download_ccp4 = partial(download_from_pdbe, 
                            base_url='https://www.ebi.ac.uk/pdbe/coordinates/files', 
                            suffix='.ccp4')
    # download_ent = partial(download_from_pdbe,
    #                        base_url='https://www.ebi.ac.uk/pdbe/entry-files/download', 
    #                        suffix='v6.ent', <---- Wrong for many
    #                       save_suffix='.pdb')
    download_cif = partial(download_from_pdbe,
                           base_url='https://www.ebi.ac.uk/pdbe/entry-files/download', 
                           suffix='.cif')
                           
## Download off PDB
                           
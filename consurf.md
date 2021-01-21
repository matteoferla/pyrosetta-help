## Consurf

Migrate [Consurf](https://consurf.tau.ac.il/) scores from the file `consurf.grades` to a Pyrosetta pose.

    from typing import *
    
    def add_bfactor_from_consurf(pose: pyrosetta.Pose, grades_filename:str, strict:bool=True):
        """
        Adds the bfactors from a consurf run to a pose based on the PDB residues.
        ``replace_res_remap_bfactors`` or ``set_bfactors_from_atom_id_map``
        were not used but may have been a cleaner strategy. This was quicker to write.
    
        The Consurf download folder, something like ``Consurf_Outputs_123456789``,
        contains a file ``consurf.grades``: this path is the ``grades_filename`` argument.
        
        when ``strict`` is true an error is raised in the residues do not match.
        """
        conscores = read_conscores(grades_filename)
        add_bfactors_from_conscores(pose, conscores, strict)
        return conscores
    
    def read_conscores(grades_filename: str) -> dict:
        conscores = {} # key is MET1:A
        with open(grades_filename) as r:
            for row in r:
                row = row.strip()
                if not len(row) or not row[0].isdigit():
                    continue
                parts = row.split()
                conscores[parts[2]] = float(parts[3])
        return conscores
        
    def add_bfactors_from_conscores(pose: pyrosetta.Pose, conscores:dict, strict:bool=True):
        """
        Adds the conscores dictionary to bfactor in place.
        """
        pdb_info = pose.pdb_info()
        pdb2pose = pdb_info.pdb2pose
        for con_res in conscores:
            pose_res = pdb2pose(res=int(con_res[3:-2]), chain=con_res[-1])
            if pose_res == 0:
                continue
            residue = pose.residue(pose_res)
            if strict:
                assert residue.name3() == con_res[:3], f'{residue.name3()} ≠ {con_res[:3]}' 
            for i in range(residue.natoms()):
                pdb_info.bfactor(pose_res,i, conscores[con_res])
            
    def remap_conscores_chains(conscores:dict, chain_map: Dict[str, str]) -> dict:
        """
        Say conscores is chain A based, but the model is not.
        """
        return {res[:-2]+':'+chain_map[res[-1]]: score for res, score in conscores.items()}
        
Basic example:

    add_bfactor_from_consurf(med27, 'Consurf_Outputs_1606141517/consurf.grades')
    
Chain is A in consurf but V in pose example:

    conscores = read_conscores('Consurf_Outputs_1606141517/consurf.grades')
    remap_conscores_chains(conscores, {'A': 'V'})
    add_bfactors_from_conscores(pose, conscores)
        
## PyMOL

b-factor putty preset is in rainbow. This is a better option:

    cartoon putty
    spectrum b, blue_white_red
    color atomic, not element C
    
Alternatively, if not for publication green to red may be better: `green_white_red`.

To disable putty:

    cartoon automatic
    color atomic
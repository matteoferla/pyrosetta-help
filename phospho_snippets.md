## PTMs

> This has not been converted into a class, it's basically a boilerplate snippet.

[Phosphosite plus](https://www.phosphosite.org/) contains an automatically generated dataset built upon curated data listing
detected post-translational modifications (PTMs).

I could, and have for [VENUS](https://venus.sgc.ox.ac.uk), use the downloaded dataset, but inspecting the values by human never hurt
so the starting point is a copy-paste, which is something like:

    '''Y222-p	
    kGYNENVyTEDGkLD 
    0	1	
    K227-ub	
    NVyTEDGkLDIWsks 
    0	3	'''


This can be converted using `ptms = convert_psp(raw)` from below:

        
    from typing import List
    import re
    
    def convert_psp(copypasted: str):
        return split_ptms(psp2list(copypasted))
    
    def psp2list(copypasted: str):
        return [part for part in copypasted.split() if '-' in part]

    def split_ptm(ptm: str):
        return {'resn': ptm[0],
               'resi': int(re.search('(\d+)', ptm).group(1)),
               'modification': ptm.split('-')[1]}
    
    def split_ptms(ptms: List[str]):
        return [split_ptm(ptm)  for ptm in ptms]
        
Apply the PTMs to a pose:

    # mod from michelanglo_protein/analyse/pyrosetta_modifier.py
    modification_name = {'p': 'phosphorylated',
                         'ac': 'acetylated',
                         'm1': 'monomethylated',
                         'm2': 'dimethylated',
                         'm3': 'trimethylated'
                        }
    # no ubiquitination or gal
    
    def apply_ptm(pose: pyrosetta.Pose, ptm:dict, chain:str='A') -> None:
        MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
        r = pose.pdb_info().pdb2pose(res=ptm['resi'], chain=chain)
        if r == 0:
            return
        residue = pose.residue(r)
        assert residue.name1() == ptm['resn'], f"Unexpected {ptm['resn']} for {residue.name1()}"
        #new_res = pyrosetta.rosetta.core.conformation.get_residue_from_name1(resn).name3()
        if ptm['modification'] not in modification_name:
            print(f'Not doing {ptm["modification"]}')
            return
        mod_name = modification_name[ptm['modification']]
        new_res = f'{residue.name3()}:{mod_name}'
        MutateResidue(target=r, new_res=new_res).apply(pose)

So using the above and minimising (do note that `MutateResidue` can repack if requested to)

    for ptm in ptms: # not too logical applying them all, but shhh.
        apply_ptm(pose, ptm, 'V')
    scorefxn = pyrosetta.get_fa_scorefxn()
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 5).apply(pose)
    pose.dump_scored_pdb('👾👾👾.r.phospho.pdb', scorefxn)
    
To inspect (if using a Jupyter notebook) via `nglview`
and a comically long list of residue names from [blog post about PTMs](http://blog.matteoferla.com/2019/01/phosphorylated-pdb-files.html)

    import nglview
    view = nglview.show_rosetta(pose)
    # this overkill list
    view.add_representation('hyperball', 'SEP or TPO or PTR or MLZ or MLY or M3L or ORN or SEC or MSE or HCS or SLZ or NLE or ALO or ALN or 2MR or CIR or ALY')
    view
    
Now, a nice thing to know is how close is a given residue to a PTM.

    import itertools, functools
    
    # v = pyrosetta.rosetta.utility.vector1_unsigned_long()
    # v.extend(resis)
    # resi_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(v)
    
    def get_distance(pose: pyrosetta.Pose, ptms: dict, chain:str='A') -> float:
        # ptm
        # PEP says its bad practise to assign an anonymous fx to variable. Schmeh.
        ptm2pose = lambda ptm: pose.pdb_info().pdb2pose(res=ptm['resi'], chain=chain)
        ptm_resis = list(filter(bool, map(ptm2pose, ptms)))
        # common xyz
        get_xyzs = lambda residue: [residue.xyz(atom_i) for atom_i in range(1, 1+residue.natoms())]
        # mut xyzs
        mut_res = pose.pdb_info().pdb2pose(res=int(row.mutation[1:-1]), chain='V')
        if mut_res == 0:
            return float('nan')
        mut_xyzs = get_xyzs(pose.residue(mut_res))
        # ptm xyzs
        ptm_xyzs = []
        for ptm_res in ptm_resis:
            ptm_xyzs.extend(get_xyzs(pose.residue(ptm_res)))
        # distances
        distance = [(mut_xyz - ptm_xyz).norm() for ptm_xyz in ptm_xyzs for mut_xyz in mut_xyzs]
        return min(distance)

Using a pandas dataframe `scores`:

    def get_pose(pdb_filename: str, params_filenames: Optional[List[str]]=None) -> pyrosetta.Pose:
            pose = pyrosetta.Pose()
            if params_filenames and isinstance(params_filenames, pyrosetta.rosetta.utility.vector1_string):
                pyrosetta.generate_nonstandard_residue_set(pose, params_filenames)
            if params_filenames and isinstance(params_filenames, list):
                params_filenames2 = pyrosetta.rosetta.utility.vector1_string()
                params_filenames2.extend(params_filenames)
                pyrosetta.generate_nonstandard_residue_set(pose, params_filenames2)
            else:
                pass
            pyrosetta.rosetta.core.import_pose.pose_from_file(pose, pdb_filename)
            return pose
             
    def filename2distance(filename:str) -> float:
        return get_distance(get_pose(filename))
        
    scores['distance_to_PTM'] = scores.filename.apply(filename2distance)
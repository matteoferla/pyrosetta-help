## Distances

Rosetta is rather fast at calculating distances between atoms —faster than PyMOL by a lot.
Whereas RDKit quickly gives a numpy matrix of distances, but loading PDBs in RDKit is very slow.

I use a list of xyzVectors of dim 3, which is a lot easier to deal with rather than a flat xyzVector of dim 3&times;_n_.
The other way would require ResidueSelector -> AtomIDs and using ``conformation.batch_get_xyz``.

    from typing import *
        
    def get_mindist(a_xyzs: List[pyrosetta.rosetta.numeric.xyzVector_double_t],
                    b_xyzs: List[pyrosetta.rosetta.numeric.xyzVector_double_t]) -> float:
        distance = [(a_xyz - b_xyz).norm() for a_xyz in a_xyzs for b_xyz in b_xyzs]
        min_d = min(distance)
        return min_d
    
    def residuevector2xyzs(pose: pyrosetta.Pose, 
                residues: pyrosetta.rosetta.core.select.residue_selector.ResidueVector) \
                -> List[pyrosetta.rosetta.numeric.xyzVector_double_t]:
        """
        conformation.batch_get_xyz returns a vector. I want a matrix.
        And it requires AtomID. Too much hassle.
        """
        return [xyz for res in residues for xyz in residue2xyzs(pose.residue(res))]
        
    def residue2xyzs(residue: pyrosetta.rosetta.core.conformation.Residue) -> List[pyrosetta.rosetta.numeric.xyzVector_double_t]:
        return [residue.xyz(atom_i) for atom_i in range(1, 1+residue.natoms())]
    

            def chain_id2xyzs(pose: pyrosetta.Pose, chain_id: int):
        rv = chain_id2resivector(pose, chain_id)
        return residuevector2xyzs(pose, rv)  
  
So `get_mindist` fed the output of `residuevector2xyzs` for each of the two two `ResidueVector`s to be compared
returns the closest distance between them.
(Remember that a `ResidueVector` is a vector of residue indices, while `apply(pose)` of a `ResidueSelector` is a bool vector: see `chain_id2resivector` below)
Most of the time it's closest distance between chains I want. So here is that code:
  
    def mindist_bewteen_chain_ids(pose, chain_id_a: int, chain_id_b: int) -> float:
        """
        Given a pose and two chain ids, get the closest distance between them.
        Calls ``chain_id2xyzs``, which in turn calls ``chain_id2resivector`` and ``residuevector2xyzs``.
        """
        return get_mindist(chain_id2xyzs(pose, chain_id_a),
                           chain_id2xyzs(pose, chain_id_b)
                          )
    
    def mindist_bewteen_chains(pose, chain_a: dict, chain_b: dict) -> float:
        """
        Accepts the chain dict from chain_ops and then calls ``mindist_bewteen_chain_ids``.
        """
        return mindist_bewteen_chain_ids(pose, chain_a['number'], chain_b['number'])
   
    def chain_id2resivector(pose: pyrosetta.Pose, chain_id: int) -> pyrosetta.rosetta.core.select.residue_selector.ResidueVector:
        sele = pyrosetta.rosetta.core.select.residue_selector.ChainSelector(chain_id)
        bv = sele.apply(pose)
        return pyrosetta.rosetta.core.select.residue_selector.ResidueVector(bv)

To get a long-form pandas DataFrame using a ChainOps object:

    import itertools
    import pandas as pd

    distances = [(first['gene_name'], second['gene_name'], mindist_bewteen_chains(pose, first, second)) for first, second in itertools.combinations(chain_ops.chains, 2)]
    distatable = pd.DataFrame(distances)
    
NB. remember humans like wide-form DataFrames!
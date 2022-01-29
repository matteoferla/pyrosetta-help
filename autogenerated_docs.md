# pyrosetta_help package

## Subpackages
* pyrosetta_help.alphafold package
 * Submodules
 * pyrosetta_help.alphafold.constraints module
 * pyrosetta_help.alphafold.multimodel module
 * pyrosetta_help.alphafold.plot module
 * pyrosetta_help.alphafold.retrieval module
 * pyrosetta_help.alphafold.superimpose module
* pyrosetta_help.blueprint_maker package
* pyrosetta_help.chain_ops package
 * Submodules
 * pyrosetta_help.chain_ops.chain_ops module
 * pyrosetta_help.chain_ops.transmogrifier module
* pyrosetta_help.common_ops package
 * Submodules
 * pyrosetta_help.common_ops.constraints module
 * pyrosetta_help.common_ops.downloads module
 * pyrosetta_help.common_ops.faux_selectors module
 * pyrosetta_help.common_ops.minimise module
 * pyrosetta_help.common_ops.nglview module
 * pyrosetta_help.common_ops.ss_changes module
 * pyrosetta_help.common_ops.utils module
* pyrosetta_help.init_ops package
 * Submodules
 * pyrosetta_help.init_ops.log module
 * pyrosetta_help.init_ops.make_options module
* pyrosetta_help.ligands package
 * Submodules
 * pyrosetta_help.ligands.hunter module
 * pyrosetta_help.ligands.load module
 * pyrosetta_help.ligands.nick module
* pyrosetta_help.per_atom package
* pyrosetta_help.residue_decription package
 * Submodules
 * pyrosetta_help.residue_decription.betaturn module
 * pyrosetta_help.residue_decription.cis module
 * pyrosetta_help.residue_decription.hbonds module
 * pyrosetta_help.residue_decription.ss module
* pyrosetta_help.score_mutants package
 * Submodules
 * pyrosetta_help.score_mutants.mutation module
 * pyrosetta_help.score_mutants.scores module
 * pyrosetta_help.score_mutants.variant module
* pyrosetta_help.threading package
* pyrosetta_help.weights package
 * Submodules
 * pyrosetta_help.weights.terms module
# pyrosetta_help.init_ops package

## Submodules

## pyrosetta_help.init_ops.log module


### pyrosetta_help.init_ops.log.configure_logger()
The function get_logger, simply adds a stringIO handler to the log and captures the log,
thus making it easier to use.
The function get_log_entries, spits out entries of a given level.
* **Returns**

    logger



### pyrosetta_help.init_ops.log.get_all_log_entries()

### pyrosetta_help.init_ops.log.get_log_entries(levelname: Union[str, int] = 20, query: Optional[str] = None)
Get a list of all entries in log at a given level.
levelname can be either an int (`logging.INFO` etc. are numbers multiples of 10 in increasing severity)
or a string of the level.
Note that it is very crude: if INFO is requested, ERROR is not shown!
* **Parameters**

**levelname** – int for the level number or str of the name
* **Returns**

    List of str


## pyrosetta_help.init_ops.make_options module


### pyrosetta_help.init_ops.make_options.make_option_string(\*\*options)
This just converts the key:value pairs to a command line string for the pyrosetta init.

>> make_option_string(no_optH=False, ex1=None, remodel=dict(blueprint=’mod.blue’))

Bools are converted,
None results in a value argument,
Tuples are converted to xx:xx type arguments
Dictionaries are converted to xx:xx type arguments (multiple, if multiple keys in the nested dictionary)

Also… Full option list: [https://www.rosettacommons.org/docs/latest/full-options-list](https://www.rosettacommons.org/docs/latest/full-options-list)
# pyrosetta_help.common_ops package

## Submodules

## pyrosetta_help.common_ops.constraints module


### pyrosetta_help.common_ops.constraints.get_NGL_selection_from_AtomID(pose: pyrosetta.rosetta.core.pose.Pose, atom_id: pyrosetta.rosetta.core.id.AtomID, named: bool = False)
Given a pyrosetta AtomID give an NGL selection.
NB `named` gives the residue name (`'[SER]3:A.CA'`) but is not a valid selection.
* **Parameters**

    
 * **pose** – 
 * **atom_id** – 
 * **named** – 
* **Returns**

    


### pyrosetta_help.common_ops.constraints.print_constraint_score(pose: pyrosetta.rosetta.core.pose.Pose, con)
Print the constraint details. atoms and score
* **Parameters**

    
 * **pose** – 
 * **con** – 
* **Returns**

    


### pyrosetta_help.common_ops.constraints.print_constraint_scores(pose: pyrosetta.rosetta.core.pose.Pose)
Prints the scores for each constraint in the pose
* **Parameters**

**pose** – 
* **Returns**

    


### pyrosetta_help.common_ops.constraints.get_AtomID(pose: pyrosetta.rosetta.core.pose.Pose, chain: str, resi: int, atomname: str)

### pyrosetta_help.common_ops.constraints.get_AtomID_by_NGL_sele(pose: pyrosetta.rosetta.core.pose.Pose, selection: str)
23:A.CA


### pyrosetta_help.common_ops.constraints.get_AtomID_from_pymol_line(pose: pyrosetta.rosetta.core.pose.Pose, line: Optional[str] = None)
Given a copypaste from the console in pymol following an atom selection in edit mode:
(`You clicked /1amq/B/A/PMP\`413/N1 -> (pk2)`) returns that atom in PyRosetta.

If the line argument is blank the clipboard is read.
* **Parameters**

    
 * **pose** – 
 * **line** – 
* **Returns**

    


### pyrosetta_help.common_ops.constraints.make_constraint_from_pymol_line(pose: pyrosetta.rosetta.core.pose.Pose, lines: str)
two atoms clicked…
You clicked /1amq/A/A/ASP\`222/OD2 -> (pk1)
You clicked /1amq/B/A/PMP\`413/N1 -> (pk2)
distance measured in PyRosetta hence the pose.
* **Parameters**

    
 * **lines** – 
 * **pose** – 
* **Returns**

    


### pyrosetta_help.common_ops.constraints.print_bad_constraint_scores(pose, cutoff=0.5)

### pyrosetta_help.common_ops.constraints.constraints2pandas(pose)
## pyrosetta_help.common_ops.downloads module


### pyrosetta_help.common_ops.downloads.download_map(code: str)
This is two functions. EMD database for cryoEM files
* **Parameters**

**code** – 
* **Returns**

    


### pyrosetta_help.common_ops.downloads.download_cif(code: str)
Download CIF. Pyrosetta has issues importing Cifs (hetatms).
This is just a note to self as PyMOL fetch can do it.


### pyrosetta_help.common_ops.downloads.download_pdb(pdbcode: str)

### pyrosetta_help.common_ops.downloads.download_opm(code: str)
Download OPM PDB

## pyrosetta_help.common_ops.faux_selectors module


### class pyrosetta_help.common_ops.faux_selectors.RingSelector(radius=12)
Bases: `object`

Select all residues in the “ring”
based upon 12A from origin.
There is probably a saner way.

NB> This is not actually a residue selector. The logical selectors will not accept it.


#### \__init__(radius=12)
Initialize self.  See help(type(self)) for accurate signature.


#### apply(pose: pyrosetta.rosetta.core.pose.Pose)

### class pyrosetta_help.common_ops.faux_selectors.AlteredSelector(threader: pyrosetta.rosetta.protocols.comparative_modeling.ThreadingMover)
Bases: `object`

Select residues that were altered in the threading.
NB> This is not actually a residue selector. The logical selectors will not accept it.


#### \__init__(threader: pyrosetta.rosetta.protocols.comparative_modeling.ThreadingMover)
Initialize self.  See help(type(self)) for accurate signature.


#### apply(pose: pyrosetta.rosetta.core.pose.Pose)

### class pyrosetta_help.common_ops.faux_selectors.UnalteredSelector(threader: pyrosetta.rosetta.protocols.comparative_modeling.ThreadingMover)
Bases: `object`

Select residues that were unaltered in the threading.
NB> This is not actually a residue selector. The logical selectors will not accept it.


#### \__init__(threader: pyrosetta.rosetta.protocols.comparative_modeling.ThreadingMover)
Initialize self.  See help(type(self)) for accurate signature.


#### apply(pose: pyrosetta.rosetta.core.pose.Pose)

### pyrosetta_help.common_ops.faux_selectors.OrListSelector(\*selectors)
OrResidueSelector but 2+
(not a class, but returns a Or
:param selectors:
:return:


### pyrosetta_help.common_ops.faux_selectors.get_bfactor_vector(pose: pyrosetta.rosetta.core.pose.Pose, cutoff: float, above=True)
Return a selection vector based on b-factors.
above = get all above. So to select bad b-factors above is `True`,
but to select AF2 bad ones. above is `False`

## pyrosetta_help.common_ops.minimise module


### pyrosetta_help.common_ops.minimise.get_local_scorefxn()

### pyrosetta_help.common_ops.minimise.prep_ED(pose: pyrosetta.rosetta.core.pose.Pose, map_filename: str)

### pyrosetta_help.common_ops.minimise.get_local_relax(scorefxn: Optional[pyrosetta.rosetta.core.scoring.ScoreFunction] = None, ncyc: int = 3, nexp: int = 2)

### pyrosetta_help.common_ops.minimise.do_local_relax(pose: pyrosetta.rosetta.core.pose.Pose, scorefxn: Optional[pyrosetta.rosetta.core.scoring.ScoreFunction] = None)

### pyrosetta_help.common_ops.minimise.do_chainwise_relax(pose: pyrosetta.rosetta.core.pose.Pose, scorefxn: Optional[pyrosetta.rosetta.core.scoring.ScoreFunction] = None, cycles: int = 5)
## pyrosetta_help.common_ops.nglview module

Adds to the NGLView Widget the method `.add_selector`

```python
import nglview as nv
view = nv.view_rosetta(pose)
view.add_selector(pose, selector)
view
```


### pyrosetta_help.common_ops.nglview.add_constraints(self: nglview.widget.NGLWidget, pose: pyrosetta.rosetta.core.pose.Pose, component=0, color='skyblue')

### pyrosetta_help.common_ops.nglview.add_rosetta(self: nglview.widget.NGLWidget, pose: pyrosetta.rosetta.core.pose.Pose, color: Optional[str] = None)
The module method `show_rosetta` creates an NGLWidget
This is a monkeypatched bound method.
* **Parameters**

    
 * **self** – 
 * **pose** – 
* **Returns**

    component



### pyrosetta_help.common_ops.nglview.add_selector(self: nglview.widget.NGLWidget, pose: pyrosetta.rosetta.core.pose.Pose, selector: pyrosetta.rosetta.core.select.residue_selector.ResidueSelector, representation_name: str = 'hyperball', color: str = 'grey', \*\*other)
Add a representation of type `representation_name` (def. ‘hyperball’) and center
based upon the Pyrosetta residue selector
* **Parameters**

    
 * **self** – 
 * **pose** – 
 * **selector** – 
 * **representation_name** – 
 * **color** – 
 * **other** – 
* **Returns**

    


### pyrosetta_help.common_ops.nglview.make_pose_comparison(self: nglview.widget.NGLWidget, first_pose: pyrosetta.rosetta.core.pose.Pose, second_pose: pyrosetta.rosetta.core.pose.Pose, first_color: str = '#00B4C4', second_color: str = '#F8766D')
Adds two objects, the first colored by default in #00B4C4, which is turquoise,

    while the second #F8766D, which is salmon.
    The poses are assumed aligned.
* **Parameters**

    
 * **self** – 
 * **first_pose** – 
 * **second_pose** – 
 * **first_color** – 
 * **second_color** – 
* **Returns**

    


### pyrosetta_help.common_ops.nglview.selector_to_ngl(self: nglview.widget.NGLWidget, pose: pyrosetta.rosetta.core.pose.Pose, selector: pyrosetta.rosetta.core.select.residue_selector.ResidueSelector)
Given a pose and a selector return the selection string for NGL.
* **Parameters**

    
 * **self** – 
 * **pose** – 
 * **selector** – 
* **Returns**

    

## pyrosetta_help.common_ops.ss_changes module

There is likely a mover that does this already and better…


### pyrosetta_help.common_ops.ss_changes.make_310_helical(pose, begin: int = 1, end: int = - 1)
> Given a pose, likely one made via `pyrosetta.pose_from_sequence`,
> make the fortran-indexed range [begin, end] (inclusive) into that SS.

In this case 3.10-helix with phi=-74.0, psi=-4.0


### pyrosetta_help.common_ops.ss_changes.make_alpha_helical(pose, begin: int = 1, end: int = - 1)
> Given a pose, likely one made via `pyrosetta.pose_from_sequence`,
> make the fortran-indexed range [begin, end] (inclusive) into that SS.

In this case alpha-helix with phi=-57.8, psi=-47.0


### pyrosetta_help.common_ops.ss_changes.make_pi_helical(pose, begin: int = 1, end: int = - 1)
> Given a pose, likely one made via `pyrosetta.pose_from_sequence`,
> make the fortran-indexed range [begin, end] (inclusive) into that SS.

In this case pi-helix with phi=-57.1, psi=-69.7


### pyrosetta_help.common_ops.ss_changes.make_sheet(pose, begin: int = 1, end: int = - 1)
> Given a pose, likely one made via `pyrosetta.pose_from_sequence`,
> make the fortran-indexed range [begin, end] (inclusive) into that SS.

In this case sheet with phi=-139, psi=+135


### pyrosetta_help.common_ops.ss_changes.make_ss(pose, begin: int = 1, end: int = - 1, phi: float = 180, psi: float = 180)
Given a pose, likely one made via `pyrosetta.pose_from_sequence`,
make the fortran-indexed range [begin, end] (inclusive) into that SS.

## pyrosetta_help.common_ops.utils module


### pyrosetta_help.common_ops.utils.pose_from_file(pdb_filename: str, params_filenames: Union[pyrosetta.rosetta.utility.vector1_std_string, List[str], None] = None)
Return a pose like pose_from_file but with params.
* **Parameters**

    
 * **pdb_filename** – 
 * **params_filenames** – 
* **Returns**

    


### pyrosetta_help.common_ops.utils.pose2pandas(pose: pyrosetta.rosetta.core.pose.Pose, scorefxn: pyrosetta.rosetta.core.scoring.ScoreFunction)
Return a pandas dataframe from the scores of the pose
* **Parameters**

**pose** – 
* **Returns**

    


### pyrosetta_help.common_ops.utils.make_blank_pose(params_filenames: Union[pyrosetta.rosetta.utility.vector1_std_string, List[str], None] = None)
Returns an empty pose, but with params from files.
* **Parameters**

**params_filenames** – 
* **Returns**

    


### pyrosetta_help.common_ops.utils.add_bfactor_from_score(pose: pyrosetta.rosetta.core.pose.Pose)
Adds the bfactors from total_score.
Snippet for testing in Jupyter

```python
import nglview as nv
view = nv.show_rosetta(pose)
# view = nv.show_file('test.cif')
view.clear_representations()
view.add_tube(radiusType="bfactor", color="bfactor", radiusScale=0.10, colorScale="RdYlBu")
view
```

`replace_res_remap_bfactors` may have been a cleaner strategy. This was quicker to write.

If this fails, it may be because the pose was not scored first.


### pyrosetta_help.common_ops.utils.get_last_res_in_chain(pose, chain='A')
Returns last residue index in a chain. THere is probably a mover that does this.
* **Parameters**

    
 * **pose** – 
 * **chain** – letter or number
* **Returns**

    


### pyrosetta_help.common_ops.utils.clarify_selector(selector: pyrosetta.rosetta.core.select.residue_selector.ResidueSelector, pose: pyrosetta.rosetta.core.pose.Pose)
Given a selector and pose return a list of residues in NGL selection format
Example, [CMP]787:H
* **Parameters**

    
 * **selector** – 
 * **pose** – 
* **Returns**

    list of residues in NGL selection format



### pyrosetta_help.common_ops.utils.count_ligands(pose: pyrosetta.rosetta.core.pose.Pose)

### pyrosetta_help.common_ops.utils.correct_numbering(pose)
A fresh PDBInfo has the PDB residue number the same as the pose one
as opposed to restarting per chain.
* **Parameters**

**pose** – 
* **Returns**

    


### pyrosetta_help.common_ops.utils.get_pdbstr(pose)

### pyrosetta_help.common_ops.utils.pose_range(pose: pyrosetta.rosetta.core.pose.Pose, protein_only=True)
range but for the pose…
For now naive of chains, non-amino acid residues etc.
# pyrosetta_help.chain_ops package

## Submodules

## pyrosetta_help.chain_ops.chain_ops module


### class pyrosetta_help.chain_ops.chain_ops.ChainOps(chains: List[dict])
Bases: `object`

Works on the list of dict (metadata.json) with among others keys
* number (pose number)
* chain (chain letter)
* gene_name (gene name)


#### \__init__(chains: List[dict])
Initialize self.  See help(type(self)) for accurate signature.


#### dump(json_filename: str = 'metadata.json')

#### get_entry(value)
guess key…
* **Parameters**

**value** – 
* **Returns**

    


#### get_entry_of_key(value, key)

#### get_pose_of_chain(pose, value)

#### load(json_filename: str = 'metadata.json')
## pyrosetta_help.chain_ops.transmogrifier module


### class pyrosetta_help.chain_ops.transmogrifier.Murinizer(chains: List[dict])
Bases: `pyrosetta_help.chain_ops.transmogrifier.Transmogrifier`

Human –> Mouse


#### \__init__(chains: List[dict])
Initialize self.  See help(type(self)) for accurate signature.


### class pyrosetta_help.chain_ops.transmogrifier.Transmogrifier(chains: List[dict], wanted_label: str, owned_label: str)
Bases: `pyrosetta_help.chain_ops.chain_ops.ChainOps`

This is to convert a mutation from one species (`owned_label`) to another (`wanted_label`)
based on provided sequences in the chains list of dict


#### \__init__(chains: List[dict], wanted_label: str, owned_label: str)
Initialize self.  See help(type(self)) for accurate signature.


#### align_seqs(chain_selection)
Align the human seq to the mouse and store it the chain dict


#### covert_A2B(seqA: str, seqB: str, resiA: int)
Given seqA and seqB as two gap aligned sequences,
and an off-by-one residue number of seqA without counting gaps,
return an off-by-one residue number of seqB without counting gaps.


#### classmethod from_chain_ops(chain_ops: pyrosetta_help.chain_ops.chain_ops.ChainOps, wanted_label: str, owned_label: str)

#### transmogrify(mutation: str, chain_selection: Union[str, int, dict])
# pyrosetta_help.alphafold package

## Submodules

## pyrosetta_help.alphafold.constraints module


### pyrosetta_help.alphafold.constraints.add_pae_constraints(pose: pyrosetta.rosetta.core.pose.Pose, errors: numpy.ndarray, cutoff: float = 12, tolerance: Optional[float] = None, weight: float = 1, adjecency_threshold=5)
Add constrains to the pose based on the errors matrix.
NB. this matrix is a reshaped version of what AF2 returns.

A harmonic function is added to CA atoms that are in residues with the error under a specified cutoff.
The mu is the current distance and the standard deviation of the harmonic is the error times `weight`.

To find out how many were added:

```python
len(pose.constraint_set().get_all_constraints())
```
* **Parameters**

    
 * **pose** – 
 * **errors** – 
 * **cutoff** – 
 * **tolerance** – if None Harmonic, if value, tollerance of FlatHarmonic
 * **weight** – this is added to the SD part so squared inverse.
 * **adjecency_threshold** – min residue separation of sequence neighbours
* **Returns**

    


### pyrosetta_help.alphafold.constraints.add_interchain_pae_constraints(pose, errors, cutoff=15)

### pyrosetta_help.alphafold.constraints.add_stretch_constraint(pose: pyrosetta.rosetta.core.pose.Pose, weight: float = 5, slope_in: float = - 0.05, residue_index_A: int = 1, residue_index_B: int = - 1, distance: Optional[float] = None, sigmoid: bool = True)
Add a constraint to “stretch out” the model, because `slope_in` is negative.
* **Parameters**

    
 * **pose** – Pose to add constraint to
 * **weight** – how strength of constraint (max of 0.5 for `SigmoidFunc`)
 * **slope_in** – negative number to stretch
 * **residue_index_A** – first residue?
 * **residue_index_B** – last residue is “-1”
 * **distance** – if omitted, the midpoint of Sigmoid will be the current distance
 * **sigmoid** – use sigmoid or identity/linear (bad idea)
* **Returns**

    


### pyrosetta_help.alphafold.constraints.get_distance_matrix(pose)
Note the distance matrix is zero indexed as it would be confusing using numpy with one indexed data.
* **Parameters**

**pose** – 
* **Returns**

    


### pyrosetta_help.alphafold.constraints.make_pae_constraint(pose, residue1_pose_idx: int, residue2_pose_idx: int, error: float, tolerance: Optional[float] = None, weight: float = 1)
## pyrosetta_help.alphafold.multimodel module


### class pyrosetta_help.alphafold.multimodel.AF2NotebookAnalyser(folder: str, load_poses: bool = True)
Bases: `object`


#### \__init__(folder: str, load_poses: bool = True)
Initialize self.  See help(type(self)) for accurate signature.


#### calculate_interface(interface='A_B')

#### constrain(groupname: str = 'relaxed', \*\*add_pae_constraints_arguments)

#### constrain_and_relax(cycles: int = 3)

#### dump(folder: Optional[str] = None)

#### dump_pdbs(groupname: str = 'relaxed', folder: Optional[str] = None, prefix: Optional[str] = '')

#### find_interface_residues()

#### get_errors()

#### get_interactions(pose: pyrosetta.rosetta.core.pose.Pose, chain_id: int, threshold: float = 3.0)

#### get_median_interface_bfactors()

#### get_poses()
This used to be a `pyrosetta.rosetta.utility.vector1_core_pose_Pose`
but it turns out that the getitem of `vector1_core_pose_Pose`
returns a clone of the pose, not the pose itself.
* **Returns**

    


#### classmethod load(folder: str, load_poses=False, params=)

#### make_AF2_dataframe()
Given a folder form ColabsFold return a dictionary
with key rank index
and value a dictionary of details

This is convoluted, but it may have been altered by a human.


#### make_phosphorylated(pdb_ptms: dict, chain: str = 'A', cycles: int = 3)

#### classmethod parse_phosphosite(raw: str, minimum: int = 1, maximum: int = - 1)
A rubbish method to convert copy-pasted Phosphosite web table.
* **Parameters**

    
 * **raw** – 
 * **minimum** – 
 * **maximum** – 
* **Returns**

    


#### property poses()

#### relax(cycles: int = 3)

#### sidechain_relax(cycles: int = 5)
## pyrosetta_help.alphafold.plot module


### pyrosetta_help.alphafold.plot.make_pae_plot(errors: numpy.ndarray)
Make AlphaFold2-EBI like PAE plot

## pyrosetta_help.alphafold.retrieval module


### pyrosetta_help.alphafold.retrieval.pose_from_alphafold2(uniprot: str)

### pyrosetta_help.alphafold.retrieval.get_alphafold2_error(uniprot: str, reshaped=True)
Returns the distances errors either as numpy matrix (`reshaped=True`)
or as the weird format from AF2-EBI —see `help(pyrosetta_help.alphafold.retrieval.reshape_errors)` for more.

Remember that the matrix is zero indexed and that these values are in Ångström
and are not pLDDT, which are stored as b-factors.


### pyrosetta_help.alphafold.retrieval.reshape_errors(errors: List[Dict[str, list]])
The JSON from AF2 has a single element list.
the sole element is a dictionary with keys ‘residue1’, ‘residue2’ and ‘distance’.
This method returns a matrix of distances reshaped based on the stated residue indices.
This is rather unlikely to differ from a regular reshape…
but idiotically I am not taking changes assuming it is always sorted.

## pyrosetta_help.alphafold.superimpose module


### pyrosetta_help.alphafold.superimpose.superimpose_by_pLDDT(pose: pyrosetta.rosetta.core.pose.Pose, original: pyrosetta.rosetta.core.pose.Pose, cutoff=70, pose_range=None)
Superimpose two poses, based on residues with pLDDT above a given threshold.
* **Parameters**

    
 * **pose** – 
 * **original** – 
 * **cutoff** – %
 * **pose_range** – optional argument to subset (start:int, end:int)
* **Returns**
# pyrosetta_help.score_mutants package

## Submodules

## pyrosetta_help.score_mutants.mutation module


### class pyrosetta_help.score_mutants.mutation.Mutation(mutation_name: str, chain: str, pose: pyrosetta.rosetta.core.pose.Pose)
Bases: `object`

A mutation is an object that has all the details of the mutation.
A variant, as interpreted in `Model.score_mutations` is a pose with a mutation.


#### \__init__(mutation_name: str, chain: str, pose: pyrosetta.rosetta.core.pose.Pose)
Initialize self.  See help(type(self)) for accurate signature.


#### assert_valid()

#### is_valid()

#### parse_mutation(mutation: str)
## pyrosetta_help.score_mutants.scores module


### pyrosetta_help.score_mutants.scores.extend_scores(scores: pandas.core.frame.DataFrame)
Adds the following fields:
* highest/lowest_contributor
* highest/lowest_contributor_value
* highest/lowest_contributor_wordy
* **Parameters**

**scores** – pd.DataFrame(output_of_variants)
* **Returns**

    


### pyrosetta_help.score_mutants.scores.get_highest_contributor(row: pandas.core.series.Series)
Not largest in abs contribution, but highest number (i.e. most positive)


### pyrosetta_help.score_mutants.scores.get_largest_contributor(row: pandas.core.series.Series)
Largest in abs amount as per the confusing fact that a very low negative number is large.


### pyrosetta_help.score_mutants.scores.get_lowest_contributor(row: pandas.core.series.Series)
Not smallest/infinitesimal in abs contribution, but lowest number (i.e. most negative)

## pyrosetta_help.score_mutants.variant module


### class pyrosetta_help.score_mutants.variant.MutantScorer(pose: pyrosetta.rosetta.core.pose.Pose, modelname: str, scorefxn: Optional[pyrosetta.rosetta.core.scoring.ScoreFunction] = None, strict_about_starting_residue: bool = True, verbose: bool = False)
Bases: `object`

Copy pasted from PI4KA <- GNB2 <- SnoopCatcher


#### CA_RMSD(poseA: pyrosetta.rosetta.core.pose.Pose, poseB: pyrosetta.rosetta.core.pose.Pose, resi: int, chain: Optional[str], distance: int)

#### FA_RMSD(poseA: pyrosetta.rosetta.core.pose.Pose, poseB: pyrosetta.rosetta.core.pose.Pose, resi: int, chain: Optional[str], distance: int)

#### \__init__(pose: pyrosetta.rosetta.core.pose.Pose, modelname: str, scorefxn: Optional[pyrosetta.rosetta.core.scoring.ScoreFunction] = None, strict_about_starting_residue: bool = True, verbose: bool = False)
Initialize self.  See help(type(self)) for accurate signature.


#### delta_scoredict(minuend: Dict[str, float], subtrahend: Dict[str, float])
minuend - subtrahend = difference
given two dict return the difference -without using pandas.


#### does_contain(mutation: Union[pyrosetta_help.score_mutants.mutation.Mutation, str], chain: Optional[str] = None)

#### classmethod from_file(filename: str, params_filenames: Optional[List[str]] = None, \*\*kwargs)

#### get_neighbour_vector(pose: pyrosetta.rosetta.core.pose.Pose, resi: int, chain: str, distance: int, include_focus_in_subset: bool = True, own_chain_only: bool = False)

#### get_present_chains(pose: Optional[pyrosetta.rosetta.core.pose.Pose] = None)

#### get_scoredict(pose: pyrosetta.rosetta.core.pose.Pose)
Given a pose get the global scores.


#### get_unweighted_scorefxn()

#### get_wscoredict(pose: pyrosetta.rosetta.core.pose.Pose)

#### has_interface(pose: pyrosetta.rosetta.core.pose.Pose, interface: str)

#### has_residue(pose: pyrosetta.rosetta.core.pose.Pose, resi: int, chain: str)

#### make_mutant(pose: pyrosetta.rosetta.core.pose.Pose, mutation: Union[str, pyrosetta_help.score_mutants.mutation.Mutation], chain='A', distance: int = 10, cycles: int = 5, inplace: bool = False)
Make a point mutant (`A23D`).
* **Parameters**

    
 * **pose** – pose
 * **mutation** – 
 * **chain** – 
 * **inplace** – 
* **Returns**

    a copy if inplace is false



#### make_output_folder()

#### movement(original: pyrosetta.rosetta.core.pose.Pose, resi: int, chain: str, distance: int, trials: int = 50, temperature: int = 1.0, replicate_number: int = 10)
This method adapted from a notebook of mine, but not from an official source, is not well written.
It should be a filter and score combo.

It returns the largest bb_rmsd of the pdb residue resi following backrub.


#### parse_mutation(mutation: Union[str, pyrosetta_help.score_mutants.mutation.Mutation], chain, pose: Optional[pyrosetta.rosetta.core.pose.Pose] = None)

#### prefix_dict(dex: Dict[str, Any], prefix: str)

#### relax_around_mover(pose: pyrosetta.rosetta.core.pose.Pose, mutation: Optional[pyrosetta_help.score_mutants.mutation.Mutation] = None, resi: Optional[int] = None, chain: Optional[str] = None, cycles=5, distance=5, own_chain_only=False)
Relaxes pose `distance` around resi:chain or mutation
* **Parameters**

    
 * **resi** – PDB residue number.
 * **chain** – 
 * **pose** – 
 * **cycles** – of relax (3 quick, 15 thorough)
 * **distance** – 
* **Returns**

    


#### score_interface(pose: pyrosetta.rosetta.core.pose.Pose, interface: str)

#### score_mutation(mutation_name: str, chains: str, distance: int, cycles: int, interfaces, ref_interface_dG: Dict[str, float], final_func: Optional[Callable] = None, preminimise: bool = False, movement: bool = False)
Scores the mutation `mutation_name` (str or Mutation instance)
returning three objects: a dict of scores, the wt (may differ from pose if preminimise=True) and mutant pose
* **Parameters**

    
 * **mutation_name** – 
 * **chains** – 
 * **distance** – 
 * **cycles** – 
 * **interfaces** – 
 * **ref_interface_dG** – premade if no preminimise.
 * **final_func** – 
 * **preminimise** – 
* **Returns**

    


#### score_mutations(mutations, chains='A', interfaces=, preminimise=False, distance=10, cycles=5, final_func: Optional[Callable] = None)

#### score_only(variant: pyrosetta.rosetta.core.pose.Pose, reference: pyrosetta.rosetta.core.pose.Pose, mutation: pyrosetta_help.score_mutants.mutation.Mutation, chains: str, distance: int, interfaces: List[Tuple[str, str]], ref_interface_dG: Dict, final_func: Optional[Callable] = None, movement: bool = False)

#### vector2list(vector: pyrosetta.rosetta.utility.vector1_bool)
# pyrosetta_help.residue_decription package

## Submodules

## pyrosetta_help.residue_decription.betaturn module


### pyrosetta_help.residue_decription.betaturn.get_betaturns(pose: pyrosetta.rosetta.core.pose.Pose)
There is a wee problem in that alpha-helices get classified as beta-turns
And I am not sure all are beta-turns…

## pyrosetta_help.residue_decription.cis module


### pyrosetta_help.residue_decription.cis.get_cis_residues(pose: pyrosetta.rosetta.core.pose.Pose)
Returns the pose indices of residues in cis (omega of zero)

## pyrosetta_help.residue_decription.hbonds module


### class pyrosetta_help.residue_decription.hbonds.BondDataType()
Bases: `dict`


#### acc_atm_name( = None)

#### acc_resi( = None)

#### acc_resn( = None)

#### distance( = None)

#### don_atm_name( = None)

#### don_resi( = None)

#### don_resn( = None)

#### energy( = None)

### pyrosetta_help.residue_decription.hbonds.get_hbond_dicts(pose: pyrosetta.rosetta.core.pose.Pose)

### pyrosetta_help.residue_decription.hbonds.hbond2dict(pose: pyrosetta.rosetta.core.pose.Pose, bond: pyrosetta.rosetta.core.scoring.hbonds.HBond)
## pyrosetta_help.residue_decription.ss module


### pyrosetta_help.residue_decription.ss.get_ss(pose: pyrosetta.rosetta.core.pose.Pose)
returns a string of the SS of types H, S, L
This is not the ALBEGO classifier
# pyrosetta_help.weights package


### class pyrosetta_help.weights.WeightWatcher()
Bases: `object`

A class to easily look at the weights for the terms of different scorefunctions.
It has various functionalities.

List available scorefunctions:

Find a scorefunction that mentions a word:

Get the comment block of a scorefunction:

Get a scorefunction by name by calling the appropriate options first.

Get weights (including ref) for a scorefunction or a name of one

Get a pandas table of the weights


#### compare(scorefxn_names: List[str], different_only: bool = False)

#### find_metion(word: str)
Find all the scorefxn names whose files contain the string word (case insensitive).

```python
find_metion('spades')
```

returns `['hydrate_score12']`


#### folder( = '/Users/matteo/miniconda3/lib/python3.7/site-packages/pyrosetta/database/scoring/weights')

#### get_ref_values_badly(scorefxn: pyrosetta.rosetta.core.scoring.ScoreFunction, prefix: Optional[bool] = True)
I could not find a direct way to get_method_weights.
Here this horrid way.


#### get_scorefxn(scorefxn_name: str)
Gets the scorefxn with appropriate corrections.


#### get_scorefxn_block(scorefxn_name: str)
Read the file given a scorefxn name


#### get_scorefxn_comments(scorefxn_name: str)
Get only the comments


#### get_weight(scorefxn: pyrosetta.rosetta.core.scoring.ScoreFunction, score_type: pyrosetta.rosetta.core.scoring.ScoreType)

#### get_weights(scorefxn: Union[str, pyrosetta.rosetta.core.scoring.ScoreFunction])
Gets weights for scorefxn


#### property possible_scorefxn_names()
Returns the scorefxn names


#### term_meanings( = {'METHOD_WEIGHTS': 'Not an energy term itself, but the parameters for each amino acid used by the ref energy term.', 'dslf_fa13': 'Disulfide geometry potential.', 'fa_atr': 'Lennard-Jones attractive between atoms in different residues (r^6 term, London dispersion forces).', 'fa_dun': "Internal energy of sidechain rotamers as derived from Dunbrack's statistics (2010 Rotamer Library used in Talaris2013).", 'fa_dun_semi': "Internal energy of sidechain semi-rotamers as derived from Dunbrack's statistics (2010 Rotamer Library used in Talaris2013).", 'fa_elec': 'Coulombic electrostatic potential with a distance-dependent dielectric.', 'fa_intra_atr_xover4': 'Intra-residue LJ attraction, counted for the atom-pairs beyond torsion-relationship.', 'fa_intra_elec': 'Intra-residue Coulombic interaction, counted for the atom-pairs beyond torsion-relationship.', 'fa_intra_rep': 'Lennard-Jones repulsive between atoms in the same residue.', 'fa_intra_rep_xover4': 'Intra-residue LJ repulsion, counted for the atom-pairs beyond torsion-relationship.', 'fa_intra_sol_xover4': 'Intra-residue LK solvation, counted for the atom-pairs beyond torsion-relationship.', 'fa_rep': 'Lennard-Jones repulsive between atoms in different residues (r^12 term, Pauli repulsion forces).', 'fa_sol': 'Lazaridis-Karplus solvation energy.', 'hbond_bb_sc': 'Sidechain-backbone hydrogen bond energy.', 'hbond_lr_bb': 'Backbone-backbone hbonds distant in primary sequence.', 'hbond_sc': 'Sidechain-sidechain hydrogen bond energy.', 'hbond_sr_bb': 'Backbone-backbone hbonds close in primary sequence.', 'hxl_tors': 'Sidechain hydroxyl group torsion preference for Ser/Thr/Tyr, supersedes yhh_planarity (that covers L- and D-Tyr only).', 'lk_ball': 'Anisotropic contribution to the solvation.', 'lk_ball_bridge': "Bonus to solvation coming from bridging waters, measured by overlap of the 'balls' from two interacting polar atoms.", 'lk_ball_bridge_uncpl': 'Same as lk_ball_bridge, but the value is uncoupled with dGfree (i.e. constant bonus, whereas lk_ball_bridge is proportional to dGfree values).', 'lk_ball_iso': 'Same as fa_sol; see below.', 'lk_ball_wtd': 'weighted sum of lk_ball & lk_ball_iso (w1\*lk_ball + w2\*lk_ball_iso); w2 is negative so that anisotropic contribution(lk_ball) replaces some portion of isotropic contribution (fa_sol=lk_ball_iso).', 'omega': 'Omega dihedral in the backbone. A Harmonic constraint on planarity with standard deviation of ~6 deg.', 'p_aa_pp': 'Probability of amino acid at Φ/Ψ.', 'pro_close': 'Proline ring closure energy and energy of psi angle of preceding residue.', 'rama': 'Ramachandran preferences.', 'rama_prepro': 'Backbone torsion preference term that takes into account of whether preceding amono acid is Proline or not.', 'ref': 'Reference energy for each amino acid. Balances internal energy of amino acid terms.  Plays role in design.', 'yhh_planarity': 'Sidechain hydroxyl group torsion preference for Tyr, superseded by hxl_tors'})
## Submodules

## pyrosetta_help.weights.terms module
# pyrosetta_help.blueprint_maker package


### class pyrosetta_help.blueprint_maker.Blueprinter(seq: str, ss: str = None)
Bases: `pyrosetta_help.blueprint_maker._init.BlueprinterInit`, `pyrosetta_help.blueprint_maker._subscripted.BlueprinterSubscripted`, `pyrosetta_help.blueprint_maker._common.BlueprinterCommon`, `pyrosetta_help.blueprint_maker._expected.BlueprinterExpected`, `pyrosetta_help.blueprint_maker._remodel.Remodel`, `pyrosetta_help.blueprint_maker._pdb_info.BlueCopier`

Make a blueprint file for Rosetta Remodel and load it into the options (`.set(fn)`).
The rows property is a list of lists, this is what the operations manipulate.
Where each row is something like [20, ‘G’, ‘.’]
subscript assignment allows the 4th field to be altered. The SS is based on the pose.
* `.from_pose(pose)` initialises it from pose, while the reg. init is from sequence and ss.

Whereas
* `.del_span(10, 29)` deletes a span
* `.wobble_span(100, 105)` wobbles a span (NATAA)
* `.mutate(5, 'W')` mutates resi i to aa
* `.insert(3, 'PIKAA W')` or `.insert(3, ['PIKAA W', 'PIKAA A'])` inserts.

a more flexible approach is using subscripts and slices:

```python
blue = Blueprinter.from_pose(pose)
blue[10:14] = 'NATAA' # preceding loop
del blue[15:20]
blue[20:25] = 'NATAA' # following loop
blue[22] = 'PIKAA W'
```

blue.set(bluprint_filename)

remember that if remodelmover is already initialised to call `rm.register_options()`

About the subscript setter, there are three (weird) things to note:
* A star will be convered into the native amino acid. PIKAA F\*CK on a glutamine, will result in PIKAA FECK
* The ranges are human style: deletion of resi 10-15 means that 6 residues are missing, including 15. Where normally in python 10:15 is 5, 15 is excluded.
* The value is a string but there are no safegards or checks, even if there are a limited amount of possibilities and it would be really nice seeing if the non-canonical for EMPTY NC XAA


#### class ResInfo(number: int, chain: str, icode: str = ' ', segmentID: str = '    ')
Bases: `object`


#### \__init__(number: int, chain: str, icode: str = ' ', segmentID: str = '    ')
Initialize self.  See help(type(self)) for accurate signature.


#### classmethod get(index: int)

#### pdb_info( = None)

#### set(index: int, pdb_info: Optional[pyrosetta.rosetta.core.pose.PDBInfo] = None)
# pyrosetta_help.ligands package

## Submodules

## pyrosetta_help.ligands.hunter module


### class pyrosetta_help.ligands.hunter.LigandHunter(sequence: str)
Bases: `object`

Given a sequence find homologues (Blast) and find if they have ligands (PDBe API).
The definition of what ligand is a cofactor comes from PDBe and does not count ions or
triphospho-nucleotides as cofactors, but does count NADH adducts.
* `.data` is a list of dictionaries.
* `.to_dataframe()` converts into a pandas dataframe for further analysis.
* `.candidate_ligands` list of ligand residue 3-letter codes


#### \__init__(sequence: str)
* **Parameters**

**sequence** – is assumed clean protein sequence



#### property candidate_ligands()

#### cofactor_codes( = ['ASC', 'F43', 'M43', 'MDO', 'PNS', '0WD', '1DG', '3AA', '3CD', '6V0', '8ID', 'A3D', 'AP0', 'CND', 'DG1', 'DN4', 'EAD', 'ENA', 'LNC', 'N01', 'NA0', 'NAD', 'NAE', 'NAI', 'NAJ', 'NAP', 'NAQ', 'NAX', 'NBD', 'NBP', 'NDC', 'NDE', 'NDO', 'NDP', 'NHD', 'NPW', 'ODP', 'P1H', 'PAD', 'SAD', 'SAE', 'SND', 'TAD', 'TAP', 'TDT', 'TXD', 'TXE', 'TXP', 'ZID', '18W', '29P', 'DPM', '2MD', 'MCN', 'MGD', 'MSS', 'MTE', 'MTQ', 'MTV', 'PCD', 'XAX', 'B12', 'CNC', 'COB', 'COY', '6FA', 'FA8', 'FAA', 'FAB', 'FAD', 'FAE', 'FAO', 'FAS', 'FCG', 'FDA', 'FED', 'FSH', 'P5F', 'RFL', 'SFD', '1YJ', 'C2F', 'FFO', 'FON', 'FOZ', 'THF', 'THG', 'THH', '01A', '01K', '0ET', '1C4', '1CV', '1CZ', '1HA', '1VU', '1XE', '2CP', '2NE', '3CP', '3H9', '3HC', '4CA', '4CO', '8JD', '8Z2', 'ACO', 'AMX', 'BCA', 'BCO', 'BSJ', 'BYC', 'CA3', 'CA5', 'CA6', 'CA8', 'CAA', 'CAJ', 'CAO', 'CIC', 'CMC', 'CMX', 'CO6', 'CO8', 'COA', 'COD', 'COF', 'COO', 'COT', 'COW', 'COZ', 'DCA', 'DCC', 'FAM', 'FCX', 'FRE', 'FYN', 'GRA', 'HAX', 'HMG', 'HSC', 'HXC', 'MCA', 'MCD', 'MDE', 'MLC', 'MYA', 'NHM', 'NHQ', 'NHW', 'NMX', 'OXK', 'S0N', 'SCA', 'SCD', 'SCO', 'SDX', 'SOP', 'T1G', 'TC6', 'WCA', 'YNC', 'ZOZ', 'SHT', 'TP7', 'TPZ', 'TXZ', 'XP8', 'XP9', '4LS', '4LU', 'FMN', 'FNR', 'FNS', 'IRF', 'RBF', 'MQ7', None, 'COM', '6HE', '7HE', 'CCH', 'COH', 'DDH', 'DHE', 'FDE', 'FMI', 'HAS', 'HDD', 'HDE', 'HEA', 'HEB', 'HEC', 'HEM', 'HIF', 'ISW', 'MH0', 'MNH', 'MNR', 'PP9', 'SH0', 'SRM', 'ZEM', 'ZNH', '4AB', '7AP', 'BHS', 'BIO', 'H2B', 'H4B', 'HBI', 'WSD', 'PQQ', 'BTI', 'BTN', 'BYT', 'DTB', 'Y7Y', 'LPA', 'LPB', '4YP', 'AT5', 'UQ1', 'UQ2', 'UQ5', 'UQ6', '0HG', '0HH', '1JO', '1JP', '1R4', '3GC', '48T', '5AU', 'ABY', 'AHE', 'ATA', 'BOB', 'BYG', 'EPY', 'ESG', 'GBI', 'GBP', 'GBX', 'GDN', 'GDS', 'GF5', 'GGC', 'GIP', 'GNB', 'GPR', 'GPS', 'GS8', 'GSB', 'GSF', 'GSH', 'GSM', 'GSN', 'GSO', 'GTB', 'GTD', 'GTS', 'GTX', 'GTY', 'GVX', 'HAG', 'IBG', 'ICY', 'JM2', 'JM5', 'JM7', 'L9X', 'LEE', 'LZ6', 'RGE', 'TGG', 'TS5', 'VWW', 'ZBF', '0AF', 'TOQ', 'TQQ', 'TRQ', '0UM', '0XU', '0Y0', '0Y1', '0Y2', '36A', '37H', '4IK', '62X', '6NR', '76H', '76J', '76K', '76L', '76M', 'EEM', 'K15', 'SA8', 'SAH', 'SAM', 'SFG', 'SMM', 'SX0', 'TT8', '1TP', '1U0', '2TP', '5GY', '8EF', '8EL', '8EO', '8FL', '8PA', 'D7K', 'EN0', 'HTL', 'M6T', 'N1T', 'N3T', 'R1T', 'S1T', 'T5X', 'T6F', 'TD6', 'TD7', 'TD8', 'TD9', 'TDK', 'TDL', 'TDM', 'TDP', 'TDW', 'THD', 'THV', 'THW', 'THY', 'TP8', 'TPP', 'TPU', 'TPW', 'TZD', 'WWF', 'MPL', 'NOP', 'NPL', 'PDP', 'PLP', 'PLR', 'PMP', 'PXP', 'PZP', 'UAH', '1TY', '2TY', 'AGQ', 'G27', 'HCC', 'P2Q', 'P3Q', 'TPQ', 'TYQ', 'TYY', 'ATP', 'GTP', 'CA', 'MG', 'W'])

#### cofactor_reference( = {'Adenosylcobalamin': [{'EC': ['1.1.1.259', '1.13.12.7', '1.14.19.30', '1.17.4.1', '1.17.4.2', '1.2.1.10', '1.2.1.17', '1.2.1.25', '1.2.1.27', '1.2.1.44', '1.2.1.52', '1.2.4.4', '1.2.7.1', '1.2.7.7', '1.21.1.2', '1.21.98.3', '1.21.99.5', '1.3.1.25', '1.3.7.1', '1.8.1.18', '2.1.1.13', '2.1.1.258', '2.1.1.308', '2.1.1.326', '2.1.3.1', '2.3.1.107', '2.3.1.148', '2.3.1.180', '2.3.1.26', '2.3.3.10', '2.3.3.5', '3.4.11.18', '4.1.1.1', '4.2.1.28', '4.2.1.30', '4.3.1.2', '4.3.1.3', '4.3.1.7', '5.4.3.3', '5.4.3.4', '5.4.3.5', '5.4.3.7', '5.4.99.1', '5.4.99.13', '5.4.99.2', '5.4.99.4', '5.4.99.63', '5.4.99.64'], 'cofactors': ['B12', 'CNC', 'COB', 'COY']}], 'Ascorbic acid': [{'EC': ['1.10.99.3', '1.11.1.11', '1.13.11.5', '1.14.11.1', '1.14.11.10', '1.14.11.11', '1.14.11.13', '1.14.11.14', '1.14.11.15', '1.14.11.17', '1.14.11.18', '1.14.11.19', '1.14.11.2', '1.14.11.20', '1.14.11.22', '1.14.11.23', '1.14.11.3', '1.14.11.4', '1.14.11.6', '1.14.11.7', '1.14.11.8', '1.14.11.9', '1.14.13.38', '1.14.14.18', '1.14.15.4', '1.14.17.1', '1.14.17.3', '1.14.17.4', '1.14.18.3', '1.14.19.1', '1.14.99.21', '1.16.5.1', '1.2.2.4', '2.4.1.202', '4.1.2.32', '5.4.99.3'], 'cofactors': ['ASC']}], 'Biopterin': [{'EC': ['1.14.13.12', '1.14.13.165', '1.14.14.47', '1.14.16.1', '1.14.16.2', '1.14.16.3', '1.14.16.4', '1.14.16.5', '1.14.16.6', '1.3.1.6'], 'cofactors': ['4AB', '7AP', 'BHS', 'BIO', 'H2B', 'H4B', 'HBI', 'WSD']}], 'Biotin': [{'EC': ['1.14.16.1', '2.1.3.1', '4.1.1.3', '4.1.1.41', '4.1.1.70', '4.1.1.74', '4.1.1.89', '6.3.4.15', '6.3.4.6', '6.4.1.1', '6.4.1.2', '6.4.1.3', '6.4.1.4', '6.4.1.5', '6.4.1.7'], 'cofactors': ['BTI', 'BTN', 'BYT', 'DTB', 'Y7Y']}], 'Coenzyme A': [{'EC': ['1.1.1.34', '1.1.1.88', '1.14.19.3', '1.17.99.3', '1.2.1.10', '1.2.1.17', '1.2.1.18', '1.2.1.2', '1.2.1.25', '1.2.1.27', '1.2.1.42', '1.2.1.44', '1.2.1.50', '1.2.1.51', '1.2.1.52', '1.2.1.57', '1.2.1.58', '1.2.1.75', '1.2.3.6', '1.2.7.1', '1.2.7.2', '1.2.7.3', '1.2.7.7', '1.2.7.8', '1.8.1.10', '1.8.1.14', '1.8.4.3', '2.3.1.1', '2.3.1.10', '2.3.1.100', '2.3.1.102', '2.3.1.104', '2.3.1.105', '2.3.1.106', '2.3.1.107', '2.3.1.108', '2.3.1.109', '2.3.1.11', '2.3.1.110', '2.3.1.111', '2.3.1.112', '2.3.1.113', '2.3.1.114', '2.3.1.115', '2.3.1.116', '2.3.1.117', '2.3.1.118', '2.3.1.119', '2.3.1.12', '2.3.1.121', '2.3.1.123', '2.3.1.125', '2.3.1.126', '2.3.1.127', '2.3.1.128', '2.3.1.13', '2.3.1.130', '2.3.1.131', '2.3.1.132', '2.3.1.133', '2.3.1.136', '2.3.1.137', '2.3.1.138', '2.3.1.139', '2.3.1.14', '2.3.1.140', '2.3.1.142', '2.3.1.144', '2.3.1.145', '2.3.1.146', '2.3.1.148', '2.3.1.15', '2.3.1.150', '2.3.1.151', '2.3.1.153', '2.3.1.154', '2.3.1.155', '2.3.1.156', '2.3.1.157', '2.3.1.159', '2.3.1.16', '2.3.1.160', '2.3.1.161', '2.3.1.162', '2.3.1.163', '2.3.1.164', '2.3.1.165', '2.3.1.166', '2.3.1.167', '2.3.1.168', '2.3.1.169', '2.3.1.17', '2.3.1.170', '2.3.1.171', '2.3.1.172', '2.3.1.173', '2.3.1.174', '2.3.1.175', '2.3.1.176', '2.3.1.177', '2.3.1.178', '2.3.1.18', '2.3.1.180', '2.3.1.182', '2.3.1.183', '2.3.1.185', '2.3.1.186', '2.3.1.188', '2.3.1.19', '2.3.1.2', '2.3.1.20', '2.3.1.205', '2.3.1.209', '2.3.1.21', '2.3.1.210', '2.3.1.22', '2.3.1.23', '2.3.1.24', '2.3.1.25', '2.3.1.255', '2.3.1.257', '2.3.1.258', '2.3.1.26', '2.3.1.27', '2.3.1.28', '2.3.1.29', '2.3.1.3', '2.3.1.30', '2.3.1.31', '2.3.1.33', '2.3.1.34', '2.3.1.36', '2.3.1.37', '2.3.1.38', '2.3.1.39', '2.3.1.4', '2.3.1.42', '2.3.1.44', '2.3.1.45', '2.3.1.46', '2.3.1.47', '2.3.1.48', '2.3.1.5', '2.3.1.50', '2.3.1.51', '2.3.1.52', '2.3.1.53', '2.3.1.54', '2.3.1.57', '2.3.1.58', '2.3.1.59', '2.3.1.6', '2.3.1.60', '2.3.1.61', '2.3.1.62', '2.3.1.63', '2.3.1.64', '2.3.1.65', '2.3.1.66', '2.3.1.67', '2.3.1.68', '2.3.1.69', '2.3.1.7', '2.3.1.70', '2.3.1.71', '2.3.1.74', '2.3.1.75', '2.3.1.76', '2.3.1.78', '2.3.1.79', '2.3.1.8', '2.3.1.80', '2.3.1.81', '2.3.1.82', '2.3.1.84', '2.3.1.85', '2.3.1.86', '2.3.1.87', '2.3.1.88', '2.3.1.89', '2.3.1.9', '2.3.1.93', '2.3.1.94', '2.3.1.95', '2.3.1.96', '2.3.1.97', '2.3.1.99', '2.3.1.B6', '2.3.3.1', '2.3.3.10', '2.3.3.11', '2.3.3.12', '2.3.3.13', '2.3.3.14', '2.3.3.2', '2.3.3.3', '2.3.3.4', '2.3.3.5', '2.3.3.6', '2.3.3.7', '2.3.3.8', '2.3.3.9', '2.7.7.23', '2.8.3.1', '2.8.3.5', '3.1.2.1', '3.1.2.10', '3.1.2.11', '3.1.2.17', '3.1.2.18', '3.1.2.19', '3.1.2.2', '3.1.2.20', '3.1.2.23', '3.1.2.25', '3.1.2.26', '3.1.2.27', '3.1.2.3', '3.1.2.4', '3.1.2.5', '3.8.1.7', '4.1.1.31', '4.1.3.36', '4.1.3.6', '5.4.99.13', '6.2.1.1', '6.2.1.10', '6.2.1.11', '6.2.1.12', '6.2.1.13', '6.2.1.14', '6.2.1.15', '6.2.1.16', '6.2.1.17', '6.2.1.18', '6.2.1.2', '6.2.1.23', '6.2.1.24', '6.2.1.25', '6.2.1.26', '6.2.1.27', '6.2.1.28', '6.2.1.3', '6.2.1.30', '6.2.1.31', '6.2.1.32', '6.2.1.33', '6.2.1.34', '6.2.1.4', '6.2.1.5', '6.2.1.6', '6.2.1.7', '6.2.1.8', '6.2.1.9'], 'cofactors': ['01A', '01K', '0ET', '1C4', '1CV', '1CZ', '1HA', '1VU', '1XE', '2CP', '2NE', '3CP', '3H9', '3HC', '4CA', '4CO', '8JD', '8Z2', 'ACO', 'AMX', 'BCA', 'BCO', 'BSJ', 'BYC', 'CA3', 'CA5', 'CA6', 'CA8', 'CAA', 'CAJ', 'CAO', 'CIC', 'CMC', 'CMX', 'CO6', 'CO8', 'COA', 'COD', 'COF', 'COO', 'COT', 'COW', 'COZ', 'DCA', 'DCC', 'FAM', 'FCX', 'FRE', 'FYN', 'GRA', 'HAX', 'HMG', 'HSC', 'HXC', 'MCA', 'MCD', 'MDE', 'MLC', 'MYA', 'NHM', 'NHQ', 'NHW', 'NMX', 'OXK', 'S0N', 'SCA', 'SCD', 'SCO', 'SDX', 'SOP', 'T1G', 'TC6', 'WCA', 'YNC', 'ZOZ']}], 'Coenzyme B': [{'EC': ['2.8.4.1'], 'cofactors': ['SHT', 'TP7', 'TPZ', 'TXZ', 'XP8', 'XP9']}], 'Coenzyme M': [{'EC': ['1.8.1.B5', '2.8.4.1'], 'cofactors': ['COM']}], 'Dipyrromethane': [{'EC': ['2.5.1.61'], 'cofactors': ['18W', '29P', 'DPM']}], 'Factor F430': [{'EC': ['2.8.4.1'], 'cofactors': ['F43', 'M43']}], 'Flavin Mononucleotide': [{'EC': ['1.1.1.27', '1.1.1.284', '1.1.1.328', '1.1.1.404', '1.1.2.3', '1.1.2.4', '1.1.3.15', '1.1.3.46', '1.1.3.6', '1.1.5.3', '1.1.99.18', '1.1.99.27', '1.1.99.31', '1.1.99.B9', '1.10.3.16', '1.10.5.1', '1.11.1.23', '1.12.1.2', '1.12.1.4', '1.12.99.6', '1.13.11.32', '1.13.11.52', '1.13.12.16', '1.13.12.4', '1.14.12.12', '1.14.12.13', '1.14.12.7', '1.14.12.8', '1.14.12.9', '1.14.13.104', '1.14.13.108', '1.14.13.109', '1.14.13.131', '1.14.13.156', '1.14.13.162', '1.14.13.165', '1.14.13.187', '1.14.13.215', '1.14.13.22', '1.14.13.28', '1.14.13.39', '1.14.13.58', '1.14.13.76', '1.14.13.B1', '1.14.14.1', '1.14.14.10', '1.14.14.11', '1.14.14.12', '1.14.14.13', '1.14.14.20', '1.14.14.21', '1.14.14.22', '1.14.14.28', '1.14.14.3', '1.14.14.30', '1.14.14.33', '1.14.14.34', '1.14.14.35', '1.14.14.47', '1.14.14.5', '1.14.14.6', '1.14.14.9', '1.14.15.1', '1.14.15.3', '1.14.19.1', '1.14.19.7', '1.14.21.1', '1.14.21.2', '1.14.99.15', '1.14.99.21', '1.14.99.24', '1.14.99.46', '1.16.1.10', '1.16.1.3', '1.16.1.4', '1.16.1.5', '1.16.1.6', '1.16.1.8', '1.16.8.1', '1.17.1.4', '1.17.1.8', '1.17.4.1', '1.17.99.3', '1.18.1.1', '1.18.1.2', '1.18.1.4', '1.18.1.8', '1.2.1.2', '1.2.1.93', '1.2.7.1', '1.21.1.1', '1.22.1.1', '1.3.1.10', '1.3.1.100', '1.3.1.108', '1.3.1.14', '1.3.1.15', '1.3.1.2', '1.3.1.31', '1.3.1.34', '1.3.1.42', '1.3.1.44', '1.3.1.6', '1.3.1.8', '1.3.1.9', '1.3.1.95', '1.3.3.1', '1.3.3.12', '1.3.3.7', '1.3.5.2', '1.3.5.3', '1.3.7.15', '1.3.8.1', '1.3.98.1', '1.3.99.24', '1.3.99.33', '1.4.1.13', '1.4.1.14', '1.4.3.1', '1.4.3.2', '1.4.3.5', '1.4.7.1', '1.5.1.20', '1.5.1.29', '1.5.1.30', '1.5.1.34', '1.5.1.36', '1.5.1.38', '1.5.1.41', '1.5.1.46', '1.5.1.B1', '1.5.3.1', '1.5.3.22', '1.5.5.2', '1.5.8.1', '1.5.8.2', '1.5.8.3', '1.5.99.1', '1.5.99.15', '1.5.99.4', '1.5.99.B2', '1.6.2.4', '1.6.3.1', '1.6.5.10', '1.6.5.11', '1.6.5.2', '1.6.5.3', '1.6.5.4', '1.6.5.5', '1.6.5.6', '1.6.5.7', '1.6.5.8', '1.6.6.9', '1.6.99.1', '1.6.99.3', '1.7.1.0', '1.7.1.16', '1.7.1.2', '1.7.1.3', '1.7.1.4', '1.7.1.6', '1.7.1.B1', '1.7.2.3', '1.7.3.1', '1.7.3.5', '1.8.1.2', '1.8.1.4', '1.8.2.3', '1.8.7.1', '1.8.99.2', '2.1.1.148', '2.1.1.21', '2.2.1.6', '2.3.1.86', '2.5.1.129', '2.5.1.15', '2.7.1.161', '2.7.1.26', '2.7.1.42', '2.7.11.1', '2.7.13.3', '2.7.7.2', '3.13.1.3', '3.6.1.18', '4.1.1.36', '4.1.1.74', '4.1.1.98', '4.1.2.32', '4.1.99.3', '4.2.3.5', '4.3.3.7', '4.99.1.7', '5.2.1.13', '5.3.3.2', '5.4.99.9', '5.5.1.18', '6.3.2.5'], 'cofactors': ['4LS', '4LU', 'FMN', 'FNR', 'FNS', 'IRF', 'RBF']}], 'Flavin adenine dinucleotide': [{'EC': ['1.1.1.158', '1.1.1.215', '1.1.1.28', '1.1.1.289', '1.1.1.47', '1.1.1.94', '1.1.1.B58', '1.1.2.3', '1.1.2.4', '1.1.3.10', '1.1.3.12', '1.1.3.13', '1.1.3.15', '1.1.3.17', '1.1.3.19', '1.1.3.20', '1.1.3.21', '1.1.3.23', '1.1.3.28', '1.1.3.3', '1.1.3.37', '1.1.3.38', '1.1.3.39', '1.1.3.4', '1.1.3.41', '1.1.3.42', '1.1.3.43', '1.1.3.44', '1.1.3.45', '1.1.3.47', '1.1.3.49', '1.1.3.5', '1.1.3.6', '1.1.3.7', '1.1.3.8', '1.1.3.9', '1.1.5.10', '1.1.5.12', '1.1.5.3', '1.1.5.4', '1.1.5.9', '1.1.98.3', '1.1.98.4', '1.1.99.1', '1.1.99.10', '1.1.99.11', '1.1.99.13', '1.1.99.16', '1.1.99.18', '1.1.99.2', '1.1.99.20', '1.1.99.21', '1.1.99.27', '1.1.99.29', '1.1.99.3', '1.1.99.37', '1.1.99.39', '1.1.99.4', '1.1.99.40', '1.1.99.6', '1.1.99.9', '1.1.99.B1', '1.1.99.B3', '1.10.3.4', '1.10.5.1', '1.10.99.2', '1.11.1.1', '1.11.1.17', '1.11.1.23', '1.12.1.2', '1.12.1.3', '1.12.1.4', '1.12.1.5', '1.12.7.2', '1.12.98.1', '1.12.98.4', '1.12.99.6', '1.13.11.17', '1.13.11.20', '1.13.11.32', '1.13.12.1', '1.13.12.16', '1.13.12.2', '1.13.12.3', '1.13.12.4', '1.13.12.9', '1.13.99.1', '1.14.11.B1', '1.14.11.B2', '1.14.12.10', '1.14.12.11', '1.14.12.12', '1.14.12.13', '1.14.12.14', '1.14.12.15', '1.14.12.17', '1.14.12.18', '1.14.12.21', '1.14.12.22', '1.14.12.24', '1.14.12.3', '1.14.12.4', '1.14.12.5', '1.14.12.7', '1.14.12.8', '1.14.13.1', '1.14.13.10', '1.14.13.101', '1.14.13.104', '1.14.13.105', '1.14.13.107', '1.14.13.108', '1.14.13.109', '1.14.13.111', '1.14.13.113', '1.14.13.114', '1.14.13.12', '1.14.13.123', '1.14.13.127', '1.14.13.129', '1.14.13.130', '1.14.13.135', '1.14.13.148', '1.14.13.149', '1.14.13.152', '1.14.13.16', '1.14.13.160', '1.14.13.162', '1.14.13.163', '1.14.13.165', '1.14.13.166', '1.14.13.167', '1.14.13.168', '1.14.13.170', '1.14.13.171', '1.14.13.18', '1.14.13.180', '1.14.13.187', '1.14.13.19', '1.14.13.195', '1.14.13.196', '1.14.13.2', '1.14.13.20', '1.14.13.200', '1.14.13.208', '1.14.13.210', '1.14.13.211', '1.14.13.212', '1.14.13.215', '1.14.13.216', '1.14.13.217', '1.14.13.218', '1.14.13.22', '1.14.13.222', '1.14.13.223', '1.14.13.224', '1.14.13.225', '1.14.13.23', '1.14.13.231', '1.14.13.232', '1.14.13.233', '1.14.13.236', '1.14.13.24', '1.14.13.25', '1.14.13.27', '1.14.13.28', '1.14.13.29', '1.14.13.3', '1.14.13.31', '1.14.13.32', '1.14.13.33', '1.14.13.35', '1.14.13.4', '1.14.13.40', '1.14.13.44', '1.14.13.5', '1.14.13.50', '1.14.13.54', '1.14.13.58', '1.14.13.59', '1.14.13.6', '1.14.13.61', '1.14.13.63', '1.14.13.64', '1.14.13.69', '1.14.13.7', '1.14.13.76', '1.14.13.78', '1.14.13.8', '1.14.13.83', '1.14.13.84', '1.14.13.9', '1.14.13.90', '1.14.13.92', '1.14.13.B30', '1.14.14.1', '1.14.14.11', '1.14.14.12', '1.14.14.15', '1.14.14.17', '1.14.14.20', '1.14.14.21', '1.14.14.27', '1.14.14.3', '1.14.14.30', '1.14.14.47', '1.14.14.5', '1.14.14.6', '1.14.14.7', '1.14.14.8', '1.14.14.9', '1.14.15.1', '1.14.15.21', '1.14.15.3', '1.14.17.3', '1.14.19.1', '1.14.19.49', '1.14.19.7', '1.14.19.9', '1.14.19.B7', '1.14.21.1', '1.14.21.2', '1.14.21.6', '1.14.99.21', '1.14.99.24', '1.14.99.34', '1.14.99.46', '1.14.99.7', '1.16.1.1', '1.16.1.3', '1.16.1.4', '1.16.1.5', '1.16.1.6', '1.16.1.7', '1.16.1.8', '1.16.1.9', '1.16.8.1', '1.17.1.4', '1.17.1.5', '1.17.3.2', '1.17.3.3', '1.17.7.2', '1.17.7.4', '1.17.8.1', '1.17.99.1', '1.17.99.3', '1.17.99.5', '1.18.1.1', '1.18.1.2', '1.18.1.3', '1.18.1.4', '1.18.1.5', '1.18.1.6', '1.18.1.7', '1.18.1.8', '1.19.1.1', '1.2.1.36', '1.2.1.51', '1.2.1.58', '1.2.1.88', '1.2.2.2', '1.2.2.4', '1.2.3.1', '1.2.3.13', '1.2.3.3', '1.2.3.4', '1.2.3.6', '1.2.3.7', '1.2.4.2', '1.2.4.4', '1.2.5.1', '1.2.5.3', '1.2.7.1', '1.2.7.4', '1.2.99.2', '1.2.99.6', '1.2.99.8', '1.2.99.9', '1.21.1.1', '1.21.3.3', '1.21.3.7', '1.21.3.8', '1.21.98.2', '1.21.99.1', '1.22.1.1', '1.3.1.10', '1.3.1.101', '1.3.1.108', '1.3.1.109', '1.3.1.110', '1.3.1.14', '1.3.1.15', '1.3.1.2', '1.3.1.27', '1.3.1.31', '1.3.1.34', '1.3.1.42', '1.3.1.44', '1.3.1.52', '1.3.1.6', '1.3.1.72', '1.3.1.86', '1.3.1.88', '1.3.1.91', '1.3.1.95', '1.3.1.98', '1.3.2.3', '1.3.3.1', '1.3.3.12', '1.3.3.15', '1.3.3.4', '1.3.3.6', '1.3.3.7', '1.3.3.8', '1.3.5.1', '1.3.5.2', '1.3.5.4', '1.3.5.5', '1.3.5.6', '1.3.7.11', '1.3.7.13', '1.3.7.8', '1.3.7.9', '1.3.8.1', '1.3.8.10', '1.3.8.11', '1.3.8.12', '1.3.8.13', '1.3.8.2', '1.3.8.3', '1.3.8.4', '1.3.8.5', '1.3.8.6', '1.3.8.7', '1.3.8.8', '1.3.8.9', '1.3.98.1', '1.3.99.1', '1.3.99.10', '1.3.99.12', '1.3.99.13', '1.3.99.17', '1.3.99.18', '1.3.99.19', '1.3.99.2', '1.3.99.21', '1.3.99.23', '1.3.99.24', '1.3.99.28', '1.3.99.3', '1.3.99.31', '1.3.99.32', '1.3.99.33', '1.3.99.37', '1.3.99.38', '1.3.99.4', '1.3.99.5', '1.3.99.7', '1.3.99.8', '1.3.99.B12', '1.3.99.B2', '1.4.1.13', '1.4.1.14', '1.4.1.5', '1.4.3.1', '1.4.3.10', '1.4.3.11', '1.4.3.12', '1.4.3.14', '1.4.3.15', '1.4.3.16', '1.4.3.19', '1.4.3.2', '1.4.3.21', '1.4.3.22', '1.4.3.23', '1.4.3.24', '1.4.3.25', '1.4.3.3', '1.4.3.4', '1.4.3.5', '1.4.3.7', '1.4.5.1', '1.4.7.1', '1.4.99.1', '1.4.99.5', '1.4.99.6', '1.4.99.B3', '1.5.1.15', '1.5.1.20', '1.5.1.29', '1.5.1.30', '1.5.1.34', '1.5.1.B1', '1.5.3.1', '1.5.3.10', '1.5.3.11', '1.5.3.12', '1.5.3.13', '1.5.3.14', '1.5.3.15', '1.5.3.16', '1.5.3.17', '1.5.3.18', '1.5.3.19', '1.5.3.2', '1.5.3.21', '1.5.3.4', '1.5.3.5', '1.5.3.6', '1.5.3.7', '1.5.5.1', '1.5.5.2', '1.5.7.1', '1.5.7.2', '1.5.8.1', '1.5.8.2', '1.5.8.3', '1.5.8.4', '1.5.99.1', '1.5.99.12', '1.5.99.13', '1.5.99.14', '1.5.99.15', '1.5.99.2', '1.5.99.3', '1.5.99.4', '1.5.99.5', '1.5.99.6', '1.5.99.8', '1.5.99.B2', '1.6.1.1', '1.6.1.2', '1.6.1.3', '1.6.1.4', '1.6.2.2', '1.6.2.4', '1.6.2.5', '1.6.2.6', '1.6.3.1', '1.6.3.2', '1.6.3.3', '1.6.3.4', '1.6.3.5', '1.6.5.10', '1.6.5.11', '1.6.5.12', '1.6.5.2', '1.6.5.3', '1.6.5.4', '1.6.5.6', '1.6.5.7', '1.6.5.8', '1.6.5.9', '1.6.6.9', '1.6.99.1', '1.6.99.3', '1.6.99.6', '1.7.1.1', '1.7.1.10', '1.7.1.15', '1.7.1.2', '1.7.1.3', '1.7.1.4', '1.7.1.5', '1.7.1.6', '1.7.2.1', '1.7.2.3', '1.7.3.1', '1.7.3.5', '1.7.99.1', '1.8.1.10', '1.8.1.11', '1.8.1.12', '1.8.1.13', '1.8.1.14', '1.8.1.15', '1.8.1.16', '1.8.1.18', '1.8.1.19', '1.8.1.2', '1.8.1.4', '1.8.1.5', '1.8.1.7', '1.8.1.8', '1.8.1.9', '1.8.1.B1', '1.8.1.B4', '1.8.2.1', '1.8.2.3', '1.8.3.2', '1.8.3.3', '1.8.3.5', '1.8.3.6', '1.8.4.10', '1.8.5.1', '1.8.5.4', '1.8.98.1', '1.8.99.2', '2.1.1.13', '2.1.1.148', '2.1.1.190', '2.1.1.229', '2.1.1.61', '2.1.1.74', '2.1.1.90', '2.2.1.6', '2.3.1.168', '2.5.1.26', '2.5.1.94', '2.7.1.161', '2.7.1.42', '2.7.11.1', '2.7.13.3', '2.7.7.2', '2.8.1.5', '2.8.1.6', '3.1.3.43', '3.13.1.3', '3.13.1.4', '3.6.1.18', '3.7.1.11', '4.1.1.36', '4.1.1.47', '4.1.1.70', '4.1.2.10', '4.1.2.47', '4.1.99.11', '4.1.99.13', '4.1.99.3', '4.2.1.120', '4.2.1.53', '4.2.1.54', '4.2.1.B25', '4.2.3.5', '4.3.3.7', '4.6.1.15', '4.99.1.2', '5.1.1.1', '5.2.1.13', '5.3.1.12', '5.3.3.13', '5.3.3.2', '5.3.3.3', '5.3.3.B2', '5.3.99.8', '5.4.3.3', '5.4.99.9', '5.5.1.18', '5.5.1.19', '5.5.1.20'], 'cofactors': ['6FA', 'FA8', 'FAA', 'FAB', 'FAD', 'FAE', 'FAO', 'FAS', 'FCG', 'FDA', 'FED', 'FSH', 'P5F', 'RFL', 'SFD']}], 'Glutathione': [{'EC': ['1.11.1.12', '1.11.1.17', '1.11.1.9', '1.13.11.18', '1.14.16.5', '1.2.1.46', '1.20.4.1', '1.20.4.2', '1.21.99.4', '1.5.4.1', '1.8.1.10', '1.8.1.16', '1.8.1.7', '1.8.1.B5', '1.8.3.3', '1.8.4.1', '1.8.4.2', '1.8.4.3', '1.8.4.4', '1.8.4.7', '1.8.4.9', '1.8.5.1', '1.8.98.2', '2.1.1.41', '2.1.1.45', '2.5.1.18', '2.8.1.3', '2.8.2.33', '3.1.2.12', '3.1.2.13', '3.1.2.6', '3.1.2.7', '3.5.1.78', '3.5.3.8', '4.4.1.20', '4.4.1.22', '4.4.1.5', '4.5.1.1', '4.5.1.3', '5.2.1.4', '5.3.4.1', '5.3.99.2', '5.3.99.3', '5.99.1.4', '6.3.2.3'], 'cofactors': ['0HG', '0HH', '1JO', '1JP', '1R4', '3GC', '48T', '5AU', 'ABY', 'AHE', 'ATA', 'BOB', 'BYG', 'EPY', 'ESG', 'GBI', 'GBP', 'GBX', 'GDN', 'GDS', 'GF5', 'GGC', 'GIP', 'GNB', 'GPR', 'GPS', 'GS8', 'GSB', 'GSF', 'GSH', 'GSM', 'GSN', 'GSO', 'GTB', 'GTD', 'GTS', 'GTX', 'GTY', 'GVX', 'HAG', 'IBG', 'ICY', 'JM2', 'JM5', 'JM7', 'L9X', 'LEE', 'LZ6', 'RGE', 'TGG', 'TS5', 'VWW', 'ZBF']}], 'Heme': [{'EC': ['1.1.1.215', '1.1.1.28', '1.1.2.2', '1.1.2.3', '1.1.2.4', '1.1.2.6', '1.1.2.7', '1.1.2.8', '1.1.2.9', '1.1.2.B5', '1.1.3.39', '1.1.3.B4', '1.1.5.11', '1.1.5.5', '1.1.5.6', '1.1.9.1', '1.1.98.1', '1.1.99.11', '1.1.99.12', '1.1.99.14', '1.1.99.18', '1.1.99.20', '1.1.99.21', '1.1.99.3', '1.1.99.4', '1.1.99.6', '1.1.99.B7', '1.10.2.2', '1.10.3.1', '1.10.3.10', '1.10.3.11', '1.10.3.12', '1.10.3.13', '1.10.3.14', '1.10.3.3', '1.10.3.9', '1.10.9.1', '1.11.1.10', '1.11.1.11', '1.11.1.13', '1.11.1.14', '1.11.1.16', '1.11.1.19', '1.11.1.21', '1.11.1.5', '1.11.1.6', '1.11.1.7', '1.11.1.8', '1.11.1.B2', '1.11.1.B7', '1.11.2.1', '1.11.2.2', '1.11.2.3', '1.11.2.4', '1.11.2.5', '1.12.2.1', '1.12.5.1', '1.12.99.6', '1.13.11.11', '1.13.11.26', '1.13.11.44', '1.13.11.49', '1.13.11.50', '1.13.11.52', '1.13.11.60', '1.13.11.9', '1.13.12.16', '1.13.12.17', '1.13.99.3', '1.14.11.9', '1.14.12.17', '1.14.13.1', '1.14.13.100', '1.14.13.103', '1.14.13.104', '1.14.13.106', '1.14.13.108', '1.14.13.109', '1.14.13.11', '1.14.13.110', '1.14.13.112', '1.14.13.115', '1.14.13.119', '1.14.13.120', '1.14.13.129', '1.14.13.13', '1.14.13.137', '1.14.13.138', '1.14.13.139', '1.14.13.140', '1.14.13.141', '1.14.13.142', '1.14.13.146', '1.14.13.147', '1.14.13.15', '1.14.13.150', '1.14.13.151', '1.14.13.152', '1.14.13.154', '1.14.13.156', '1.14.13.161', '1.14.13.17', '1.14.13.177', '1.14.13.181', '1.14.13.184', '1.14.13.188', '1.14.13.197', '1.14.13.198', '1.14.13.199', '1.14.13.2', '1.14.13.202', '1.14.13.203', '1.14.13.21', '1.14.13.221', '1.14.13.25', '1.14.13.28', '1.14.13.30', '1.14.13.37', '1.14.13.39', '1.14.13.4', '1.14.13.41', '1.14.13.42', '1.14.13.47', '1.14.13.48', '1.14.13.49', '1.14.13.5', '1.14.13.52', '1.14.13.53', '1.14.13.55', '1.14.13.56', '1.14.13.57', '1.14.13.6', '1.14.13.60', '1.14.13.67', '1.14.13.68', '1.14.13.7', '1.14.13.70', '1.14.13.71', '1.14.13.72', '1.14.13.73', '1.14.13.74', '1.14.13.75', '1.14.13.76', '1.14.13.77', '1.14.13.79', '1.14.13.80', '1.14.13.81', '1.14.13.85', '1.14.13.86', '1.14.13.87', '1.14.13.88', '1.14.13.89', '1.14.13.91', '1.14.13.93', '1.14.13.94', '1.14.13.95', '1.14.13.96', '1.14.13.97', '1.14.13.98', '1.14.13.99', '1.14.13.B14', '1.14.13.B21', '1.14.13.B3', '1.14.14.1', '1.14.14.14', '1.14.14.16', '1.14.14.17', '1.14.14.18', '1.14.14.19', '1.14.14.23', '1.14.14.24', '1.14.14.25', '1.14.14.29', '1.14.14.32', '1.14.14.36', '1.14.14.40', '1.14.14.41', '1.14.14.42', '1.14.14.44', '1.14.14.45', '1.14.14.47', '1.14.15.1', '1.14.15.10', '1.14.15.13', '1.14.15.15', '1.14.15.16', '1.14.15.18', '1.14.15.19', '1.14.15.20', '1.14.15.3', '1.14.15.4', '1.14.15.5', '1.14.15.6', '1.14.15.8', '1.14.15.9', '1.14.18.1', '1.14.18.6', '1.14.18.7', '1.14.18.8', '1.14.19.1', '1.14.19.12', '1.14.19.16', '1.14.19.17', '1.14.19.18', '1.14.19.19', '1.14.19.20', '1.14.19.22', '1.14.19.25', '1.14.19.26', '1.14.19.29', '1.14.19.3', '1.14.19.30', '1.14.19.31', '1.14.19.33', '1.14.19.34', '1.14.19.38', '1.14.19.4', '1.14.19.44', '1.14.19.47', '1.14.19.51', '1.14.19.6', '1.14.21.1', '1.14.21.10', '1.14.21.11', '1.14.21.2', '1.14.21.3', '1.14.21.4', '1.14.21.5', '1.14.21.7', '1.14.21.9', '1.14.99.1', '1.14.99.10', '1.14.99.14', '1.14.99.15', '1.14.99.19', '1.14.99.22', '1.14.99.24', '1.14.99.28', '1.14.99.37', '1.14.99.39', '1.14.99.4', '1.14.99.45', '1.14.99.49', '1.14.99.57', '1.14.99.9', '1.14.99.B1', '1.15.1.2', '1.16.5.1', '1.16.9.1', '1.17.1.4', '1.17.2.1', '1.17.2.2', '1.17.99.1', '1.17.99.2', '1.18.1.8', '1.2.1.70', '1.2.2.1', '1.2.3.7', '1.2.5.2', '1.20.2.1', '1.20.99.1', '1.21.98.2', '1.3.1.6', '1.3.2.3', '1.3.3.10', '1.3.3.9', '1.3.5.1', '1.3.5.4', '1.3.98.3', '1.3.99.1', '1.4.9.1', '1.4.9.2', '1.5.1.20', '1.5.3.1', '1.5.99.3', '1.5.99.4', '1.5.99.5', '1.5.99.6', '1.6.3.1', '1.6.5.10', '1.7.1.1', '1.7.1.14', '1.7.1.15', '1.7.1.2', '1.7.1.3', '1.7.1.4', '1.7.2.1', '1.7.2.2', '1.7.2.3', '1.7.2.5', '1.7.2.6', '1.7.2.7', '1.7.2.8', '1.7.3.4', '1.7.3.6', '1.7.5.1', '1.7.5.2', '1.7.6.1', '1.7.7.1', '1.7.7.2', '1.7.7.4', '1.7.99.7', '1.8.1.2', '1.8.1.3', '1.8.2.1', '1.8.2.2', '1.8.2.3', '1.8.2.4', '1.8.2.5', '1.8.3.1', '1.8.7.1', '1.8.98.1', '1.8.98.3', '1.8.99.2', '1.8.99.3', '1.8.99.5', '1.8.99.B2', '1.9.3.1', '1.9.6.1', '1.9.98.1', '1.97.1.1', '1.97.1.9', '2.1.1.304', '2.7.13.3', '2.7.7.65', '3.1.4.52', '4.1.1.B11', '4.2.1.121', '4.2.1.22', '4.2.1.92', '4.2.1.B9', '4.3.1.26', '4.6.1.2', '4.99.1.1', '4.99.1.3', '4.99.1.5', '4.99.1.7', '5.3.99.3', '5.3.99.4', '5.3.99.5', '5.3.99.6', '5.4.4.5', '5.4.4.6', '6.6.1.2'], 'cofactors': ['6HE', '7HE', 'CCH', 'COH', 'DDH', 'DHE', 'FDE', 'FMI', 'HAS', 'HDD', 'HDE', 'HEA', 'HEB', 'HEC', 'HEM', 'HIF', 'ISW', 'MH0', 'MNH', 'MNR', 'PP9', 'SH0', 'SRM', 'ZEM', 'ZNH']}], 'Lipoic acid': [{'EC': ['1.1.1.5', '1.17.4.2', '1.2.4.1 ', '1.2.4.2', '1.2.4.4', '1.4.4.2', '1.8.1.4', '2.3.1.61'], 'cofactors': ['LPA', 'LPB']}], 'MIO': [{'EC': ['4.3.1.23', '4.3.1.24', '4.3.1.25', '5.4.3.6'], 'cofactors': ['MDO']}], 'Menaquinone': [{'EC': ['1.1.5.4', '1.1.5.6', '1.17.99.1', '1.21.99.5', '1.3.5.4', '1.7.5.1', '1.8.5.5', '3.4.21.21', '3.4.21.6', '3.4.21.69', '4.1.1.90'], 'cofactors': ['MQ7']}], 'Molybdopterin': [{'EC': ['1.1.1.288', '1.1.5.6', '1.1.99.33', '1.12.98.4', '1.17.1.4', '1.17.1.5', '1.17.2.1', '1.17.3.2', '1.17.3.3', '1.17.5.1', '1.17.99.2', '1.2.1.2', '1.2.2.4', '1.2.3.1', '1.2.3.8', '1.2.5.3', '1.2.7.4', '1.2.7.5', '1.2.7.6', '1.2.7.7', '1.2.99.6', '1.2.99.8', '1.2.99.9', '1.20.9.1', '1.20.98.1', '1.3.7.9', '1.3.99.16', '1.3.99.17', '1.5.99.14', '1.5.99.4', '1.7.2.3', '1.7.5.1', '1.7.7.2', '1.8.2.1', '1.8.2.4', '1.8.3.1', '1.8.5.3', '1.9.6.1', '1.97.1.1', '1.97.1.9', '4.2.1.112'], 'cofactors': ['2MD', 'MCN', 'MGD', 'MSS', 'MTE', 'MTQ', 'MTV', 'PCD', 'XAX']}], 'NA': [{'EC': [None], 'cofactors': [None]}], 'Nicotinamide-adenine dinucleotide': [{'EC': ['1.1.1.1', '1.1.1.10', '1.1.1.100', '1.1.1.101', '1.1.1.102', '1.1.1.103', '1.1.1.104', '1.1.1.105', '1.1.1.106', '1.1.1.107', '1.1.1.108', '1.1.1.11', '1.1.1.110', '1.1.1.111', '1.1.1.112', '1.1.1.113', '1.1.1.114', '1.1.1.115', '1.1.1.116', '1.1.1.117', '1.1.1.118', '1.1.1.119', '1.1.1.12', '1.1.1.120', '1.1.1.121', '1.1.1.122', '1.1.1.123', '1.1.1.124', '1.1.1.125', '1.1.1.126', '1.1.1.127', '1.1.1.128', '1.1.1.129', '1.1.1.13', '1.1.1.130', '1.1.1.131', '1.1.1.132', '1.1.1.133', '1.1.1.134', '1.1.1.135', '1.1.1.136', '1.1.1.137', '1.1.1.138', '1.1.1.14', '1.1.1.140', '1.1.1.141', '1.1.1.142', '1.1.1.143', '1.1.1.144', '1.1.1.145', '1.1.1.146', '1.1.1.147', '1.1.1.148', '1.1.1.149', '1.1.1.15', '1.1.1.150', '1.1.1.151', '1.1.1.152', '1.1.1.153', '1.1.1.154', '1.1.1.156', '1.1.1.157', '1.1.1.158', '1.1.1.159', '1.1.1.16', '1.1.1.160', '1.1.1.161', '1.1.1.162', '1.1.1.163', '1.1.1.164', '1.1.1.165', '1.1.1.166', '1.1.1.167', '1.1.1.168', '1.1.1.169', '1.1.1.17', '1.1.1.170', '1.1.1.172', '1.1.1.173', '1.1.1.174', '1.1.1.175', '1.1.1.176', '1.1.1.177', '1.1.1.178', '1.1.1.179', '1.1.1.18', '1.1.1.181', '1.1.1.183', '1.1.1.184', '1.1.1.185', '1.1.1.186', '1.1.1.187', '1.1.1.188', '1.1.1.189', '1.1.1.19', '1.1.1.190', '1.1.1.191', '1.1.1.192', '1.1.1.193', '1.1.1.194', '1.1.1.195', '1.1.1.196', '1.1.1.197', '1.1.1.198', '1.1.1.199', '1.1.1.2', '1.1.1.20', '1.1.1.200', '1.1.1.201', '1.1.1.202', '1.1.1.203', '1.1.1.205', '1.1.1.206', '1.1.1.207', '1.1.1.208', '1.1.1.209', '1.1.1.21', '1.1.1.210', '1.1.1.211', '1.1.1.212', '1.1.1.213', '1.1.1.214', '1.1.1.215', '1.1.1.216', '1.1.1.217', '1.1.1.218', '1.1.1.219', '1.1.1.22', '1.1.1.220', '1.1.1.221', '1.1.1.222', '1.1.1.223', '1.1.1.224', '1.1.1.225', '1.1.1.226', '1.1.1.227', '1.1.1.228', '1.1.1.229', '1.1.1.23', '1.1.1.230', '1.1.1.231', '1.1.1.232', '1.1.1.233', '1.1.1.234', '1.1.1.235', '1.1.1.236', '1.1.1.237', '1.1.1.238', '1.1.1.239', '1.1.1.24', '1.1.1.240', '1.1.1.241', '1.1.1.243', '1.1.1.244', '1.1.1.245', '1.1.1.246', '1.1.1.247', '1.1.1.248', '1.1.1.25', '1.1.1.250', '1.1.1.251', '1.1.1.252', '1.1.1.254', '1.1.1.255', '1.1.1.256', '1.1.1.257', '1.1.1.258', '1.1.1.259', '1.1.1.26', '1.1.1.260', '1.1.1.261', '1.1.1.262', '1.1.1.263', '1.1.1.264', '1.1.1.265', '1.1.1.266', '1.1.1.267', '1.1.1.268', '1.1.1.269', '1.1.1.27', '1.1.1.270', '1.1.1.271', '1.1.1.272', '1.1.1.273', '1.1.1.274', '1.1.1.275', '1.1.1.276', '1.1.1.277', '1.1.1.278', '1.1.1.279', '1.1.1.28', '1.1.1.280', '1.1.1.281', '1.1.1.282', '1.1.1.283', '1.1.1.284', '1.1.1.285', '1.1.1.286', '1.1.1.287', '1.1.1.288', '1.1.1.289', '1.1.1.29', '1.1.1.290', '1.1.1.291', '1.1.1.292', '1.1.1.294', '1.1.1.295', '1.1.1.296', '1.1.1.297', '1.1.1.298', '1.1.1.299', '1.1.1.3', '1.1.1.30', '1.1.1.300', '1.1.1.301', '1.1.1.302', '1.1.1.303', '1.1.1.304', '1.1.1.305', '1.1.1.306', '1.1.1.307', '1.1.1.308', '1.1.1.309', '1.1.1.31', '1.1.1.311', '1.1.1.312', '1.1.1.314', '1.1.1.315', '1.1.1.316', '1.1.1.32', '1.1.1.320', '1.1.1.321', '1.1.1.322', '1.1.1.323', '1.1.1.324', '1.1.1.326', '1.1.1.327', '1.1.1.328', '1.1.1.329', '1.1.1.33', '1.1.1.331', '1.1.1.332', '1.1.1.333', '1.1.1.334', '1.1.1.335', '1.1.1.336', '1.1.1.337', '1.1.1.34', '1.1.1.340', '1.1.1.342', '1.1.1.343', '1.1.1.344', '1.1.1.345', '1.1.1.346', '1.1.1.347', '1.1.1.348', '1.1.1.35', '1.1.1.350', '1.1.1.351', '1.1.1.352', '1.1.1.354', '1.1.1.356', '1.1.1.357', '1.1.1.359', '1.1.1.36', '1.1.1.360', '1.1.1.361', '1.1.1.363', '1.1.1.364', '1.1.1.365', '1.1.1.366', '1.1.1.367', '1.1.1.368', '1.1.1.369', '1.1.1.37', '1.1.1.370', '1.1.1.371', '1.1.1.373', '1.1.1.375', '1.1.1.376', '1.1.1.377', '1.1.1.378', '1.1.1.379', '1.1.1.38', '1.1.1.380', '1.1.1.381', '1.1.1.382', '1.1.1.383', '1.1.1.384', '1.1.1.386', '1.1.1.387', '1.1.1.388', '1.1.1.389', '1.1.1.39', '1.1.1.390', '1.1.1.391', '1.1.1.394', '1.1.1.395', '1.1.1.396', '1.1.1.398', '1.1.1.399', '1.1.1.4', '1.1.1.40', '1.1.1.400', '1.1.1.401', '1.1.1.402', '1.1.1.403', '1.1.1.404', '1.1.1.406', '1.1.1.41', '1.1.1.410', '1.1.1.411', '1.1.1.42', '1.1.1.43', '1.1.1.44', '1.1.1.45', '1.1.1.46', '1.1.1.47', '1.1.1.48', '1.1.1.49', '1.1.1.5', '1.1.1.50', '1.1.1.51', '1.1.1.52', '1.1.1.53', '1.1.1.54', '1.1.1.55', '1.1.1.56', '1.1.1.57', '1.1.1.58', '1.1.1.59', '1.1.1.6', '1.1.1.60', '1.1.1.61', '1.1.1.62', '1.1.1.63', '1.1.1.64', '1.1.1.65', '1.1.1.66', '1.1.1.67', '1.1.1.69', '1.1.1.7', '1.1.1.71', '1.1.1.72', '1.1.1.73', '1.1.1.75', '1.1.1.76', '1.1.1.77', '1.1.1.78', '1.1.1.79', '1.1.1.8', '1.1.1.80', '1.1.1.81', '1.1.1.82', '1.1.1.83', '1.1.1.84', '1.1.1.85', '1.1.1.86', '1.1.1.87', '1.1.1.88', '1.1.1.9', '1.1.1.90', '1.1.1.91', '1.1.1.92', '1.1.1.93', '1.1.1.94', '1.1.1.95', '1.1.1.96', '1.1.1.97', '1.1.1.98', '1.1.1.99', '1.1.1.B18', '1.1.1.B19', '1.1.1.B20', '1.1.1.B25', '1.1.1.B28', '1.1.1.B3', '1.1.1.B35', '1.1.1.B38', '1.1.1.B4', '1.1.1.B40', '1.1.1.B47', '1.1.1.B52', '1.1.1.B58', '1.1.1.B60', '1.1.1.B61', '1.1.3.16', '1.1.99.14', '1.1.99.2', '1.1.99.21', '1.1.99.24', '1.1.99.28', '1.1.99.36', '1.1.99.37', '1.1.99.7', '1.10.1.1', '1.10.5.1', '1.11.1.1', '1.11.1.2', '1.11.1.23', '1.12.1.2', '1.12.1.3', '1.12.1.4', '1.12.2.1', '1.12.7.2', '1.12.98.4', '1.12.99.6', '1.13.11.20', '1.13.11.29', '1.13.11.30', '1.13.11.79', '1.13.12.14', '1.13.12.17', '1.14.12.1', '1.14.12.10', '1.14.12.11', '1.14.12.12', '1.14.12.13', '1.14.12.14', '1.14.12.15', '1.14.12.16', '1.14.12.17', '1.14.12.18', '1.14.12.19', '1.14.12.20', '1.14.12.21', '1.14.12.22', '1.14.12.23', '1.14.12.24', '1.14.12.3', '1.14.12.4', '1.14.12.5', '1.14.12.7', '1.14.12.8', '1.14.12.9', '1.14.13.1', '1.14.13.10', '1.14.13.100', '1.14.13.101', '1.14.13.102', '1.14.13.103', '1.14.13.104', '1.14.13.105', '1.14.13.106', '1.14.13.107', '1.14.13.108', '1.14.13.109', '1.14.13.11', '1.14.13.110', '1.14.13.111', '1.14.13.112', '1.14.13.113', '1.14.13.114', '1.14.13.115', '1.14.13.12', '1.14.13.123', '1.14.13.127', '1.14.13.128', '1.14.13.129', '1.14.13.13', '1.14.13.130', '1.14.13.131', '1.14.13.135', '1.14.13.14', '1.14.13.141', '1.14.13.142', '1.14.13.15', '1.14.13.151', '1.14.13.158', '1.14.13.16', '1.14.13.161', '1.14.13.162', '1.14.13.163', '1.14.13.165', '1.14.13.166', '1.14.13.167', '1.14.13.17', '1.14.13.172', '1.14.13.178', '1.14.13.179', '1.14.13.18', '1.14.13.180', '1.14.13.182', '1.14.13.19', '1.14.13.196', '1.14.13.2', '1.14.13.20', '1.14.13.202', '1.14.13.205', '1.14.13.209', '1.14.13.21', '1.14.13.210', '1.14.13.211', '1.14.13.212', '1.14.13.215', '1.14.13.216', '1.14.13.217', '1.14.13.218', '1.14.13.22', '1.14.13.220', '1.14.13.222', '1.14.13.224', '1.14.13.225', '1.14.13.23', '1.14.13.230', '1.14.13.24', '1.14.13.25', '1.14.13.26', '1.14.13.27', '1.14.13.28', '1.14.13.29', '1.14.13.3', '1.14.13.30', '1.14.13.31', '1.14.13.32', '1.14.13.33', '1.14.13.34', '1.14.13.35', '1.14.13.36', '1.14.13.37', '1.14.13.38', '1.14.13.39', '1.14.13.4', '1.14.13.40', '1.14.13.41', '1.14.13.42', '1.14.13.43', '1.14.13.44', '1.14.13.46', '1.14.13.47', '1.14.13.48', '1.14.13.49', '1.14.13.5', '1.14.13.50', '1.14.13.51', '1.14.13.52', '1.14.13.53', '1.14.13.54', '1.14.13.55', '1.14.13.56', '1.14.13.57', '1.14.13.58', '1.14.13.59', '1.14.13.6', '1.14.13.60', '1.14.13.61', '1.14.13.62', '1.14.13.63', '1.14.13.64', '1.14.13.66', '1.14.13.67', '1.14.13.68', '1.14.13.69', '1.14.13.7', '1.14.13.70', '1.14.13.71', '1.14.13.72', '1.14.13.73', '1.14.13.74', '1.14.13.75', '1.14.13.76', '1.14.13.77', '1.14.13.78', '1.14.13.79', '1.14.13.8', '1.14.13.80', '1.14.13.81', '1.14.13.82', '1.14.13.83', '1.14.13.84', '1.14.13.85', '1.14.13.86', '1.14.13.87', '1.14.13.88', '1.14.13.89', '1.14.13.9', '1.14.13.90', '1.14.13.91', '1.14.13.92', '1.14.13.93', '1.14.13.94', '1.14.13.95', '1.14.13.96', '1.14.13.97', '1.14.13.98', '1.14.13.99', '1.14.13.B1', '1.14.13.B28', '1.14.14.1', '1.14.14.10', '1.14.14.16', '1.14.14.18', '1.14.14.20', '1.14.14.21', '1.14.14.22', '1.14.14.37', '1.14.14.47', '1.14.14.9', '1.14.15.1', '1.14.15.15', '1.14.15.21', '1.14.15.3', '1.14.15.4', '1.14.16.5', '1.14.18.2', '1.14.18.3', '1.14.18.4', '1.14.18.8', '1.14.19.1', '1.14.19.17', '1.14.19.20', '1.14.19.21', '1.14.19.22', '1.14.19.3', '1.14.19.30', '1.14.19.5', '1.14.19.7', '1.14.21.1', '1.14.21.2', '1.14.21.3', '1.14.21.4', '1.14.21.5', '1.14.21.6', '1.14.21.7', '1.14.99.14', '1.14.99.15', '1.14.99.19', '1.14.99.2', '1.14.99.20', '1.14.99.21', '1.14.99.22', '1.14.99.24', '1.14.99.26', '1.14.99.3', '1.14.99.31', '1.14.99.32', '1.14.99.34', '1.14.99.40', '1.14.99.9', '1.14.99.B1', '1.16.1.1', '1.16.1.10', '1.16.1.2', '1.16.1.3', '1.16.1.4', '1.16.1.5', '1.16.1.6', '1.16.1.7', '1.16.1.8', '1.16.1.9', '1.17.1.1', '1.17.1.2', '1.17.1.3', '1.17.1.4', '1.17.1.5', '1.17.1.8', '1.17.3.2', '1.17.7.4', '1.17.99.1', '1.18.1.1', '1.18.1.2', '1.18.1.3', '1.18.1.4', '1.18.1.5', '1.18.1.6', '1.18.1.7', '1.18.1.8', '1.19.1.1', '1.2.1.10', '1.2.1.11', '1.2.1.12', '1.2.1.13', '1.2.1.15', '1.2.1.16', '1.2.1.17', '1.2.1.18', '1.2.1.19', '1.2.1.2', '1.2.1.20', '1.2.1.21', '1.2.1.22', '1.2.1.23', '1.2.1.24', '1.2.1.25', '1.2.1.26', '1.2.1.27', '1.2.1.28', '1.2.1.29', '1.2.1.3', '1.2.1.30', '1.2.1.31', '1.2.1.32', '1.2.1.33', '1.2.1.36', '1.2.1.38', '1.2.1.39', '1.2.1.4', '1.2.1.40', '1.2.1.41', '1.2.1.42', '1.2.1.43', '1.2.1.44', '1.2.1.45', '1.2.1.46', '1.2.1.47', '1.2.1.48', '1.2.1.49', '1.2.1.5', '1.2.1.50', '1.2.1.51', '1.2.1.52', '1.2.1.53', '1.2.1.54', '1.2.1.57', '1.2.1.58', '1.2.1.59', '1.2.1.60', '1.2.1.61', '1.2.1.62', '1.2.1.63', '1.2.1.64', '1.2.1.65', '1.2.1.66', '1.2.1.67', '1.2.1.68', '1.2.1.69', '1.2.1.7', '1.2.1.70', '1.2.1.71', '1.2.1.72', '1.2.1.73', '1.2.1.74', '1.2.1.75', '1.2.1.76', '1.2.1.77', '1.2.1.78', '1.2.1.79', '1.2.1.8', '1.2.1.80', '1.2.1.81', '1.2.1.83', '1.2.1.85', '1.2.1.87', '1.2.1.88', '1.2.1.89', '1.2.1.9', '1.2.1.90', '1.2.1.91', '1.2.1.92', '1.2.1.93', '1.2.1.94', '1.2.1.96', '1.2.1.97', '1.2.1.98', '1.2.1.99', '1.2.1.B29', '1.2.1.B8', '1.2.3.1', '1.2.4.1', '1.2.4.2', '1.2.4.4', '1.2.5.3', '1.2.98.1', '1.20.1.1', '1.21.1.1', '1.21.1.2', '1.21.4.1', '1.21.4.3', '1.21.4.4', '1.22.1.1', '1.3.1.1', '1.3.1.10', '1.3.1.100', '1.3.1.101', '1.3.1.102', '1.3.1.105', '1.3.1.106', '1.3.1.107', '1.3.1.109', '1.3.1.11', '1.3.1.111', '1.3.1.112', '1.3.1.12', '1.3.1.13', '1.3.1.14', '1.3.1.15', '1.3.1.16', '1.3.1.17', '1.3.1.18', '1.3.1.19', '1.3.1.2', '1.3.1.20', '1.3.1.21', '1.3.1.22', '1.3.1.24', '1.3.1.25', '1.3.1.26', '1.3.1.27', '1.3.1.28', '1.3.1.29', '1.3.1.3', '1.3.1.30', '1.3.1.31', '1.3.1.32', '1.3.1.33', '1.3.1.34', '1.3.1.35', '1.3.1.36', '1.3.1.37', '1.3.1.38', '1.3.1.39', '1.3.1.4', '1.3.1.40', '1.3.1.41', '1.3.1.42', '1.3.1.43', '1.3.1.44', '1.3.1.45', '1.3.1.46', '1.3.1.47', '1.3.1.48', '1.3.1.49', '1.3.1.5', '1.3.1.51', '1.3.1.52', '1.3.1.53', '1.3.1.54', '1.3.1.56', '1.3.1.57', '1.3.1.58', '1.3.1.6', '1.3.1.60', '1.3.1.62', '1.3.1.63', '1.3.1.64', '1.3.1.65', '1.3.1.66', '1.3.1.67', '1.3.1.68', '1.3.1.69', '1.3.1.7', '1.3.1.70', '1.3.1.71', '1.3.1.72', '1.3.1.73', '1.3.1.74', '1.3.1.75', '1.3.1.76', '1.3.1.77', '1.3.1.78', '1.3.1.79', '1.3.1.8', '1.3.1.80', '1.3.1.81', '1.3.1.82', '1.3.1.83', '1.3.1.84', '1.3.1.85', '1.3.1.86', '1.3.1.9', '1.3.1.92', '1.3.1.95', '1.3.1.98', '1.3.1.99', '1.3.3.3', '1.3.3.9', '1.3.5.2', '1.3.7.12', '1.3.7.2', '1.3.7.3', '1.3.8.1', '1.3.98.3', '1.3.99.23', '1.3.99.5', '1.3.99.6', '1.4.1.1', '1.4.1.10', '1.4.1.11', '1.4.1.12', '1.4.1.13', '1.4.1.14', '1.4.1.15', '1.4.1.16', '1.4.1.17', '1.4.1.18', '1.4.1.19', '1.4.1.2', '1.4.1.20', '1.4.1.21', '1.4.1.23', '1.4.1.3', '1.4.1.4', '1.4.1.5', '1.4.1.7', '1.4.1.8', '1.4.1.9', '1.4.1.B2', '1.4.1.B4', '1.4.1.B5', '1.4.99.6', '1.5.1.1', '1.5.1.10', '1.5.1.11', '1.5.1.12', '1.5.1.15', '1.5.1.16', '1.5.1.17', '1.5.1.18', '1.5.1.19', '1.5.1.2', '1.5.1.20', '1.5.1.21', '1.5.1.22', '1.5.1.23', '1.5.1.24', '1.5.1.25', '1.5.1.26', '1.5.1.27', '1.5.1.28', '1.5.1.29', '1.5.1.3', '1.5.1.30', '1.5.1.31', '1.5.1.32', '1.5.1.33', '1.5.1.34', '1.5.1.36', '1.5.1.37', '1.5.1.38', '1.5.1.40', '1.5.1.42', '1.5.1.43', '1.5.1.46', '1.5.1.47', '1.5.1.49', '1.5.1.5', '1.5.1.51', '1.5.1.6', '1.5.1.7', '1.5.1.8', '1.5.1.9', '1.5.1.B1', '1.5.1.B5', '1.5.3.1', '1.5.99.5', '1.6.1.1', '1.6.1.2', '1.6.2.2', '1.6.2.4', '1.6.2.5', '1.6.2.6', '1.6.3.1', '1.6.3.2', '1.6.3.3', '1.6.3.4', '1.6.5.10', '1.6.5.11', '1.6.5.2', '1.6.5.3', '1.6.5.4', '1.6.5.5', '1.6.5.6', '1.6.5.7', '1.6.5.8', '1.6.5.9', '1.6.6.9', '1.6.99.1', '1.6.99.3', '1.6.99.5', '1.6.99.6', '1.7.1.1', '1.7.1.10', '1.7.1.11', '1.7.1.12', '1.7.1.13', '1.7.1.14', '1.7.1.15', '1.7.1.16', '1.7.1.2', '1.7.1.3', '1.7.1.4', '1.7.1.5', '1.7.1.6', '1.7.1.7', '1.7.1.9', '1.7.1.B1', '1.7.1.B4', '1.7.2.3', '1.7.2.5', '1.8.1.10', '1.8.1.11', '1.8.1.12', '1.8.1.13', '1.8.1.14', '1.8.1.15', '1.8.1.16', '1.8.1.17', '1.8.1.18', '1.8.1.2', '1.8.1.3', '1.8.1.4', '1.8.1.5', '1.8.1.6', '1.8.1.7', '1.8.1.8', '1.8.1.9', '1.8.1.B1', '1.8.4.13', '1.8.4.14', '2.1.1.147', '2.1.1.148', '2.1.1.201', '2.1.1.341', '2.1.1.74', '2.1.2.13', '2.3.1.111', '2.3.1.119', '2.3.1.161', '2.3.1.165', '2.3.1.170', '2.3.1.190', '2.3.1.201', '2.3.1.203', '2.3.1.85', '2.3.1.86', '2.3.1.94', '2.4.2.19', '2.4.2.30', '2.4.2.31', '2.4.2.36', '2.4.2.B12', '2.4.2.B14', '2.4.2.B15', '2.4.2.B16', '2.5.1.21', '2.5.1.44', '2.5.1.45', '2.5.1.46', '2.7.1.160', '2.7.1.19', '2.7.1.23', '2.7.1.86', '2.7.3.2', '2.7.8.11', '2.8.1.6', '3.1.2.1', '3.13.1.1', '3.2.1.122', '3.2.1.20', '3.2.1.22', '3.2.1.49', '3.2.1.86', '3.3.1.1', '3.5.1.10', '3.5.1.98', '3.5.4.9', '3.7.1.11', '4.1.1.101', '4.1.1.3', '4.1.1.35', '4.1.1.71', '4.1.1.73', '4.1.1.75', '4.1.1.B4', '4.1.2.32', '4.1.3.30', '4.1.3.39', '4.1.99.5', '4.2.1.10', '4.2.1.115', '4.2.1.135', '4.2.1.17', '4.2.1.45', '4.2.1.46', '4.2.1.47', '4.2.1.49', '4.2.1.76', '4.2.1.93', '4.2.3.124', '4.2.3.152', '4.2.3.154', '4.2.3.155', '4.2.3.4', '4.3.1.12', '4.3.1.22', '4.3.1.28', '4.4.1.17', '5.1.2.3', '5.1.3.10', '5.1.3.12', '5.1.3.13', '5.1.3.16', '5.1.3.18', '5.1.3.2', '5.1.3.20', '5.1.3.5', '5.1.3.6', '5.1.3.7', '5.2.1.6', '5.3.1.12', '5.3.3.1', '5.3.3.2', '5.4.99.5', '5.4.99.9', '5.5.1.3', '5.5.1.4', '6.2.1.13', '6.3.2.11', '6.5.1.1', '6.5.1.2'], 'cofactors': ['0WD', '1DG', '3AA', '3CD', '6V0', '8ID', 'A3D', 'AP0', 'CND', 'DG1', 'DN4', 'EAD', 'ENA', 'LNC', 'N01', 'NA0', 'NAD', 'NAE', 'NAI', 'NAJ', 'NAP', 'NAQ', 'NAX', 'NBD', 'NBP', 'NDC', 'NDE', 'NDO', 'NDP', 'NHD', 'NPW', 'ODP', 'P1H', 'PAD', 'SAD', 'SAE', 'SND', 'TAD', 'TAP', 'TDT', 'TXD', 'TXE', 'TXP', 'ZID']}], 'Orthoquinone residues (LTQ, TTQ, CTQ)': [{'EC': ['1.4.99.3'], 'cofactors': ['0AF', 'TOQ', 'TQQ', 'TRQ']}], 'Phosphopantetheine': [{'EC': ['1.1.1.100', '1.3.1.10', '2.3.1.244', '2.3.1.38', '2.3.1.39', '2.3.1.41', '2.3.1.85', '2.3.1.86', '2.7.8.7', '3.1.2.14', '4.2.1.61', '5.1.1.11', '6.2.1.10'], 'cofactors': ['PNS']}], "Pyridoxal 5'-phosphate": [{'EC': ['1.17.1.1', '1.4.1.12', '1.4.3.21', '1.4.3.22', '1.4.4.2', '2.1.2.1', '2.1.2.10', '2.1.2.5', '2.1.2.7', '2.2.1.8', '2.3.1.201', '2.3.1.263', '2.3.1.29', '2.3.1.30', '2.3.1.37', '2.3.1.47', '2.3.1.50', '2.4.1.1', '2.4.1.8', '2.5.1.113', '2.5.1.140', '2.5.1.47', '2.5.1.48', '2.5.1.49', '2.5.1.51', '2.5.1.52', '2.5.1.53', '2.5.1.65', '2.5.1.73', '2.5.1.76', '2.6.1.1', '2.6.1.100', '2.6.1.101', '2.6.1.102', '2.6.1.103', '2.6.1.104', '2.6.1.105', '2.6.1.106', '2.6.1.107', '2.6.1.108', '2.6.1.109', '2.6.1.11', '2.6.1.110', '2.6.1.111', '2.6.1.112', '2.6.1.113', '2.6.1.12', '2.6.1.13', '2.6.1.14', '2.6.1.15', '2.6.1.17', '2.6.1.18', '2.6.1.19', '2.6.1.2', '2.6.1.21', '2.6.1.22', '2.6.1.24', '2.6.1.26', '2.6.1.27', '2.6.1.28', '2.6.1.29', '2.6.1.3', '2.6.1.30', '2.6.1.33', '2.6.1.34', '2.6.1.35', '2.6.1.36', '2.6.1.37', '2.6.1.38', '2.6.1.39', '2.6.1.4', '2.6.1.40', '2.6.1.41', '2.6.1.42', '2.6.1.43', '2.6.1.44', '2.6.1.45', '2.6.1.46', '2.6.1.47', '2.6.1.48', '2.6.1.49', '2.6.1.5', '2.6.1.50', '2.6.1.51', '2.6.1.52', '2.6.1.55', '2.6.1.56', '2.6.1.57', '2.6.1.58', '2.6.1.59', '2.6.1.6', '2.6.1.60', '2.6.1.62', '2.6.1.63', '2.6.1.64', '2.6.1.65', '2.6.1.66', '2.6.1.67', '2.6.1.7', '2.6.1.71', '2.6.1.72', '2.6.1.75', '2.6.1.76', '2.6.1.77', '2.6.1.78', '2.6.1.79', '2.6.1.8', '2.6.1.80', '2.6.1.81', '2.6.1.82', '2.6.1.83', '2.6.1.84', '2.6.1.87', '2.6.1.88', '2.6.1.89', '2.6.1.9', '2.6.1.90', '2.6.1.92', '2.6.1.96', '2.6.1.98', '2.6.1.99', '2.6.1.B16', '2.6.1.B17', '2.6.1.B3', '2.6.1.B6', '2.6.99.3', '2.8.1.4', '2.8.1.7', '2.8.1.9', '2.9.1.1', '2.9.1.2', '3.5.1.61', '3.5.1.9', '3.5.99.7', '3.7.1.3', '4.1.1.11', '4.1.1.12', '4.1.1.14', '4.1.1.15', '4.1.1.16', '4.1.1.17', '4.1.1.18', '4.1.1.19', '4.1.1.20', '4.1.1.22', '4.1.1.24', '4.1.1.25', '4.1.1.28', '4.1.1.29', '4.1.1.50', '4.1.1.53', '4.1.1.57', '4.1.1.64', '4.1.1.65', '4.1.1.81', '4.1.1.86', '4.1.1.95', '4.1.1.96', '4.1.2.26', '4.1.2.27', '4.1.2.42', '4.1.2.48', '4.1.2.49', '4.1.2.5', '4.1.3.38', '4.1.3.41', '4.1.99.1', '4.1.99.2', '4.2.1.122', '4.2.1.144', '4.2.1.145', '4.2.1.164', '4.2.1.168', '4.2.1.20', '4.2.1.22', '4.2.1.50', '4.2.3.1', '4.2.3.134', '4.2.3.2', '4.3.1.13', '4.3.1.15', '4.3.1.16', '4.3.1.17', '4.3.1.18', '4.3.1.19', '4.3.1.20', '4.3.1.27', '4.3.1.9', '4.4.1.1', '4.4.1.10', '4.4.1.11', '4.4.1.13', '4.4.1.14', '4.4.1.15', '4.4.1.16', '4.4.1.2', '4.4.1.25', '4.4.1.28', '4.4.1.35', '4.4.1.4', '4.4.1.6', '4.4.1.8', '4.4.1.9', '4.5.1.2', '4.5.1.5', '4.99.1.6', '5.1.1.1', '5.1.1.10', '5.1.1.12', '5.1.1.13', '5.1.1.15', '5.1.1.17', '5.1.1.18', '5.1.1.2', '5.1.1.21', '5.1.1.3', '5.1.1.5', '5.1.1.6', '5.1.1.9', '5.4.3.2', '5.4.3.3', '5.4.3.4', '5.4.3.5', '5.4.3.8', '5.4.3.9'], 'cofactors': ['MPL', 'NOP', 'NPL', 'PDP', 'PLP', 'PLR', 'PMP', 'PXP', 'PZP', 'UAH']}], 'Pyrroloquinoline Quinone': [{'EC': ['1.1.2.7', '1.1.5.2', '1.1.99.1', '1.1.99.20', '1.1.99.22', '1.1.99.23', '1.1.99.8', '1.2.99.3', '1.4.99.3'], 'cofactors': ['PQQ']}], 'S-adenosylmethionine': [{'EC': ['1.1.98.6', '1.16.1.8', '1.17.4.2', '1.3.99.22', '1.97.1.4', '2.1.1.1', '2.1.1.10', '2.1.1.100', '2.1.1.101', '2.1.1.102', '2.1.1.103', '2.1.1.104', '2.1.1.105', '2.1.1.106', '2.1.1.107', '2.1.1.108', '2.1.1.109', '2.1.1.11', '2.1.1.110', '2.1.1.111', '2.1.1.112', '2.1.1.113', '2.1.1.114', '2.1.1.115', '2.1.1.116', '2.1.1.117', '2.1.1.118', '2.1.1.119', '2.1.1.12', '2.1.1.120', '2.1.1.121', '2.1.1.122', '2.1.1.123', '2.1.1.124', '2.1.1.125', '2.1.1.126', '2.1.1.127', '2.1.1.128', '2.1.1.129', '2.1.1.13', '2.1.1.130', '2.1.1.131', '2.1.1.132', '2.1.1.133', '2.1.1.136', '2.1.1.137', '2.1.1.139', '2.1.1.140', '2.1.1.141', '2.1.1.142', '2.1.1.143', '2.1.1.144', '2.1.1.145', '2.1.1.146', '2.1.1.147', '2.1.1.149', '2.1.1.15', '2.1.1.150', '2.1.1.151', '2.1.1.152', '2.1.1.153', '2.1.1.154', '2.1.1.155', '2.1.1.156', '2.1.1.157', '2.1.1.158', '2.1.1.159', '2.1.1.16', '2.1.1.160', '2.1.1.161', '2.1.1.162', '2.1.1.163', '2.1.1.164', '2.1.1.165', '2.1.1.166', '2.1.1.167', '2.1.1.168', '2.1.1.169', '2.1.1.17', '2.1.1.170', '2.1.1.171', '2.1.1.173', '2.1.1.174', '2.1.1.175', '2.1.1.176', '2.1.1.178', '2.1.1.179', '2.1.1.18', '2.1.1.180', '2.1.1.182', '2.1.1.183', '2.1.1.184', '2.1.1.185', '2.1.1.186', '2.1.1.187', '2.1.1.189', '2.1.1.190', '2.1.1.192', '2.1.1.193', '2.1.1.197', '2.1.1.199', '2.1.1.2', '2.1.1.20', '2.1.1.200', '2.1.1.201', '2.1.1.202', '2.1.1.203', '2.1.1.204', '2.1.1.206', '2.1.1.207', '2.1.1.210', '2.1.1.212', '2.1.1.214', '2.1.1.215', '2.1.1.216', '2.1.1.219', '2.1.1.22', '2.1.1.220', '2.1.1.221', '2.1.1.222', '2.1.1.224', '2.1.1.226', '2.1.1.227', '2.1.1.228', '2.1.1.229', '2.1.1.230', '2.1.1.231', '2.1.1.232', '2.1.1.233', '2.1.1.234', '2.1.1.235', '2.1.1.236', '2.1.1.237', '2.1.1.240', '2.1.1.241', '2.1.1.243', '2.1.1.244', '2.1.1.25', '2.1.1.255', '2.1.1.26', '2.1.1.260', '2.1.1.261', '2.1.1.262', '2.1.1.263', '2.1.1.264', '2.1.1.265', '2.1.1.266', '2.1.1.267', '2.1.1.27', '2.1.1.271', '2.1.1.273', '2.1.1.274', '2.1.1.278', '2.1.1.28', '2.1.1.281', '2.1.1.29', '2.1.1.290', '2.1.1.293', '2.1.1.294', '2.1.1.295', '2.1.1.296', '2.1.1.297', '2.1.1.298', '2.1.1.308', '2.1.1.309', '2.1.1.31', '2.1.1.310', '2.1.1.317', '2.1.1.32', '2.1.1.325', '2.1.1.328', '2.1.1.33', '2.1.1.333', '2.1.1.337', '2.1.1.34', '2.1.1.35', '2.1.1.36', '2.1.1.37', '2.1.1.38', '2.1.1.39', '2.1.1.4', '2.1.1.40', '2.1.1.41', '2.1.1.42', '2.1.1.43', '2.1.1.44', '2.1.1.46', '2.1.1.47', '2.1.1.48', '2.1.1.49', '2.1.1.50', '2.1.1.51', '2.1.1.52', '2.1.1.53', '2.1.1.55', '2.1.1.56', '2.1.1.57', '2.1.1.59', '2.1.1.6', '2.1.1.60', '2.1.1.61', '2.1.1.62', '2.1.1.64', '2.1.1.65', '2.1.1.66', '2.1.1.67', '2.1.1.68', '2.1.1.69', '2.1.1.7', '2.1.1.70', '2.1.1.71', '2.1.1.72', '2.1.1.75', '2.1.1.76', '2.1.1.77', '2.1.1.78', '2.1.1.79', '2.1.1.8', '2.1.1.80', '2.1.1.82', '2.1.1.83', '2.1.1.84', '2.1.1.85', '2.1.1.87', '2.1.1.88', '2.1.1.89', '2.1.1.9', '2.1.1.91', '2.1.1.94', '2.1.1.95', '2.1.1.96', '2.1.1.97', '2.1.1.98', '2.1.1.99', '2.1.1.B109', '2.1.1.B74', '2.1.1.B75', '2.1.1.B76', '2.1.1.B84', '2.1.1.B85', '2.1.1.B93', '2.3.1.161', '2.3.1.184', '2.3.1.231', '2.4.99.17', '2.5.1.108', '2.5.1.120', '2.5.1.24', '2.5.1.25', '2.5.1.38', '2.5.1.4', '2.5.1.43', '2.5.1.63', '2.5.1.77', '2.6.1.62', '2.8.1.6', '2.8.1.8', '2.8.4.3', '2.8.4.4', '2.8.4.5', '3.1.21.3', '3.1.21.5', '3.3.1.2', '4.1.1.50', '4.1.99.14', '4.1.99.15', '4.1.99.17', '4.1.99.22', '4.2.1.22', '4.3.1.30', '4.3.99.3', '4.4.1.14', '4.7.1.1', '5.4.3.2', '5.4.99.58'], 'cofactors': ['0UM', '0XU', '0Y0', '0Y1', '0Y2', '36A', '37H', '4IK', '62X', '6NR', '76H', '76J', '76K', '76L', '76M', 'EEM', 'K15', 'SA8', 'SAH', 'SAM', 'SFG', 'SMM', 'SX0', 'TT8']}], 'Tetrahydrofolic acid': [{'EC': ['1.13.11.52', '1.14.13.165', '1.14.14.47', '1.14.16.1', '1.14.16.2', '1.14.16.4', '1.14.16.5', '1.14.18.1', '1.5.1.6', '1.5.3.1', '1.5.8.3', '2.1.1.13', '2.1.1.148', '2.1.1.19', '2.1.1.269', '2.1.1.341', '2.1.1.45', '2.1.2.1', '2.1.2.10', '2.1.2.11', '2.1.2.13', '2.1.2.2', '2.1.2.3', '2.1.2.5', '2.1.2.8', '2.1.2.9', '3.5.3.8', '4.1.99.3'], 'cofactors': ['1YJ', 'C2F', 'FFO', 'FON', 'FOZ', 'THF', 'THG', 'THH']}], 'Thiamine diphosphate': [{'EC': ['1.2.1.51', '1.2.1.58', '1.2.2.2', '1.2.3.13', '1.2.3.3', '1.2.4.1', '1.2.4.2', '1.2.4.4', '1.2.5.1', '1.2.7.1', '1.2.7.10', '1.2.7.11', '1.2.7.3', '1.2.7.8', '2.2.1.1', '2.2.1.12', '2.2.1.3', '2.2.1.4', '2.2.1.5', '2.2.1.6', '2.2.1.7', '2.2.1.9', '2.3.1.12', '2.3.1.190', '2.3.3.15', '2.4.1.12', '2.5.1.66', '2.7.1.78', '2.8.1.6', '3.7.1.11', '3.7.1.2', '4.1.1.1', '4.1.1.43', '4.1.1.47', '4.1.1.50', '4.1.1.7', '4.1.1.71', '4.1.1.72', '4.1.1.74', '4.1.1.75', '4.1.1.79', '4.1.1.8', '4.1.1.82', '4.1.2.22', '4.1.2.35', '4.1.2.38', '4.1.2.9', '4.2.1.113', '4.2.1.3'], 'cofactors': ['1TP', '1U0', '2TP', '5GY', '8EF', '8EL', '8EO', '8FL', '8PA', 'D7K', 'EN0', 'HTL', 'M6T', 'N1T', 'N3T', 'R1T', 'S1T', 'T5X', 'T6F', 'TD6', 'TD7', 'TD8', 'TD9', 'TDK', 'TDL', 'TDM', 'TDP', 'TDW', 'THD', 'THV', 'THW', 'THY', 'TP8', 'TPP', 'TPU', 'TPW', 'TZD', 'WWF']}], 'Topaquinone': [{'EC': ['1.4.3.13', '1.4.3.21', '1.4.3.22'], 'cofactors': ['1TY', '2TY', 'AGQ', 'G27', 'HCC', 'P2Q', 'P3Q', 'TPQ', 'TYQ', 'TYY']}], 'Ubiquinone': [{'EC': ['1.1.5.3', '1.1.5.4', '1.1.5.5', '1.1.99.1', '1.1.99.20', '1.1.99.22', '1.1.99.6', '1.10.2.2', '1.10.3.10', '1.17.5.2', '1.2.5.1', '1.2.5.3', '1.3.1.22', '1.3.5.1', '1.3.5.2', '1.3.5.5', '1.3.98.1', '1.4.5.1', '1.5.5.1', '1.6.5.3', '1.6.5.8', '1.8.1.4', '5.3.4.1'], 'cofactors': ['4YP', 'AT5', 'UQ1', 'UQ2', 'UQ5', 'UQ6']}]})

#### get_most_common_ligands()
Uses collections.counter
* **Returns**

    


#### get_pdb_entry_by_ligand(ligand_code: str)
get pdb **entry** by ligand.
Returns the first, which should be lowest e-value


#### property ligand_data()
Data for the ligands.

Note that LigandNicker has a get smiles method, that works like this, but is unrelated.


#### to_dataframe()
Converts `.data` to a pandas dataframe

## pyrosetta_help.ligands.load module


### pyrosetta_help.ligands.load.parameterised_pose_from_file(pdb_filename, wanted_ligands: Union[List[str], Dict[str, Optional[str]]] = , force_parameterisation: bool = False, neutralise_params: bool = True, save_params: bool = True, overriding_params=)
pose loading, the circutous way to not loose ligand or use PDB_component.
Assumes all mystery components are PDB residues.
Works best with ignore_unrecognized_res False
* **Parameters**

    
 * **pdb_filename** – 
 * **wanted_ligands** – a list of three letter codes or a dictionary of three letter codes to None or smiles
 * **force_parameterisation** – 
 * **neutralise_params** – protonated for pH 7
 * **save_params** – 
 * **overriding_params** – list of params filenames
* **Returns**

    


### pyrosetta_help.ligands.load.parameterised_pose_from_pdbblock(pdbblock: str, wanted_ligands: Union[List[str], Dict[str, Optional[str]]] = , force_parameterisation: bool = False, neutralise_params: bool = True, save_params: bool = True, overriding_params=)
pose loading, the circutous way to not loose ligand or use PDB_component.
Assumes all mystery components are PDB residues.
Works best with ignore_unrecognized_res False
* **Parameters**

    
 * **pdb_filename** – 
 * **wanted_ligands** – a list of three letter codes or a dictionary of three letter codes to None or smiles
 * **force_parameterisation** – 
 * **neutralise_params** – protonated for pH 7
 * **save_params** – 
 * **overriding_params** – list of params filenames
* **Returns**

    


### pyrosetta_help.ligands.load.get_smiles(ligand_code: str)
Get the smiles of a ligand.
Remember that PDBe smiles need to charged to pH 7.

## pyrosetta_help.ligands.nick module


### class pyrosetta_help.ligands.nick.LigandNicker(pdb_filename: Optional[str] = None, pdb_filehandle: Optional[io.IOBase] = None, pdb_block: Optional[str] = None, pose: Optional[pyrosetta.rosetta.core.pose.Pose] = None, chain: str = 'A', wanted_ligands: List[str] = , force_parameterisation: bool = False, neutralise_params: bool = True, save_params: bool = True, overriding_params=)
Bases: `object`

Given a pdb_file (regular initialisation) or code (`.from_pdbcode` classmethod)
and a list of 3-letter codes of wanted residues (`wanted_ligands` argument),
it loads it as Pyrosetta Pose (`donor_pose`), ready for `migrate` to loads the acceptor_pose and nick
the residues that are wanted.

If `force_parameterisation` is on it or the residue is novel it parameterises it.


#### \__init__(pdb_filename: Optional[str] = None, pdb_filehandle: Optional[io.IOBase] = None, pdb_block: Optional[str] = None, pose: Optional[pyrosetta.rosetta.core.pose.Pose] = None, chain: str = 'A', wanted_ligands: List[str] = , force_parameterisation: bool = False, neutralise_params: bool = True, save_params: bool = True, overriding_params=)
Initialisation loads the donor. `migrate` loads the acceptor.

Unfortunately, I originally wrote it to use pdb_filename only.
I should have written the init to accept a pose and
a class method the filename. And now some chucks of out there in the wild use pdb_filename.
So I could not switch it to something generic. Now it accepts four possible alternative choices:
pdb_filename, pdb_filehandle, pdb_block, pose.
* **Parameters**

    
 * **pdb_filename** – 
 * **pdb_filehandle** – 
 * **pdb_block** – 
 * **pose** – 
 * **chain** – 
 * **wanted_ligands** – 
 * **force_parameterisation** – 
 * **neutralise_params** – pH 7 protonation
 * **save_params** – 
 * **overriding_params** – overide paramaterisation and use provide params



#### constrain_migrated(wanted_vector, dex)

#### classmethod from_pdbcode(pdb_code: str, chain: str, \*args, \*\*kvargs)

#### get_mapping_between_poses(index_list: List[int])
Given a list of indices for one pose, return their aligned equivalents in the second.
Modded to be donor_pose –> acceptor_pose


#### get_surrounding_residue(pose, chain_filter, wanted_selector)

#### get_wanted_selector()

#### make_atomID_map(dex, query_pose, target_pose)

#### make_constraint_foreign_hbond(hbond, dex: Dict[int, int])

#### migrate(acceptor_pose: pyrosetta.rosetta.core.pose.Pose, acceptor_chain: str = 'A', constrained: bool = True, relaxed: bool = True, relax_radius: int = 20, relax_cycles: int = 3)
The acceptor pose is the non-empty pose.

This method aligns the sequences of the acceptor and donor pose.
It finds the mapping of the neighbourhood of the wanted residues of the donor_pose
It superimposes the poses by those residues.
It adds the residues that need nicking.
It adds constraints (optionally) based on the hydrogen bonding of the residues around the wanted residues.
onto the `.acceptor_pose`. To check:

```python
print( len(self.acceptor_pose.constraint_set().get_all_constraints())  )
```

It then optionally relaxes the neighbourhood.
* **Parameters**

    
 * **acceptor_pose** – 
 * **acceptor_chain** – 
 * **constrained** – 
 * **relaxed** – 
 * **relax_radius** – 
 * **relax_cycles** – 
* **Returns**

    


#### relax_migrated(distance: int = 20, cycles: int = 3, atom_pair_weight: int = 5)

### pyrosetta_help.ligands.nick.chain_letter_to_number(letter, pose)
# pyrosetta_help.per_atom package


### class pyrosetta_help.per_atom.AtomicInteractions(pose: pyrosetta.rosetta.core.pose.Pose, target_idx: int, threshold: int = 3, scorefxn: Optional[pyrosetta.rosetta.core.scoring.ScoreFunction] = None, weighted: bool = True, halved: bool = False)
Bases: `object`

Gets the per atom energies for the per_atom.

```python
ai = AtomicInteractions(pose, 1)  # noqa  - pose residue index 1
print(ai.describe_best())
ai.per_atom[(' N  ', 2, ' N  ')]
```

Unfortunately, bonding is not taken into account therefore the total is more favourable
as these are ignored.

```python
ai.total, ai.expected_total
```


#### \__init__(pose: pyrosetta.rosetta.core.pose.Pose, target_idx: int, threshold: int = 3, scorefxn: Optional[pyrosetta.rosetta.core.scoring.ScoreFunction] = None, weighted: bool = True, halved: bool = False)
Initialize self.  See help(type(self)) for accurate signature.


#### property best_interactions()

#### describe_atom(residue, atomname)

#### describe_best()

#### describe_interaction(target_atomname, other_resi, other_atomname)

#### property expected_total()

#### get_cc_selector()

#### get_target_selector()

#### score_types( = ['fa_atr', 'fa_rep', 'fa_sol', 'fa_elec'])

#### term_relevance_cutoff( = 0.5)

#### property total()

### class pyrosetta_help.per_atom.NeighbourInteractions(\*args, \*\*kwargs)
Bases: `object`


#### \__init__(\*args, \*\*kwargs)
Initialize self.  See help(type(self)) for accurate signature.
# pyrosetta_help.threading package


### pyrosetta_help.threading.get_alignment(target: str, template: str)
Returns alignments.


### pyrosetta_help.threading.write_grishin(target_name, target_sequence, template_name, template_sequence, outfile)

### pyrosetta_help.threading.thread(target_sequence: str, template_pose: pyrosetta.rosetta.core.pose.Pose, target_pose: Optional[pyrosetta.rosetta.core.pose.Pose] = None, target_name: str = 'target', template_name: str = 'template', fragment_sets: Optional[pyrosetta.rosetta.utility.vector1_std_shared_ptr_core_fragment_FragSet_t] = None, align: Optional[pyrosetta.rosetta.core.sequence.SequenceAlignment] = None)
Given the target sequence and the optional blank target pose (in case there’s some reason for it),
thread the sequence against the template pose — which is assumed to be a single chain.
Optionally using fragments from fragment_sets.
The three outputs are the target_pose, the threader instance and a vector of the residues threaded.

```python
print(threader.frag_libs()[1].nr_frames())
```

```python
qt = threader.get_qt_mapping(threaded)
print(f'Template residue {21} is residue {q[21]} in target')
```
* **Parameters**

    
 * **target_sequence** – 
 * **template_pose** – 
 * **target_pose** – 
 * **target_name** – 
 * **template_name** – 
 * **fragment_sets** – 
 * **align** – 
* **Returns**

    


### pyrosetta_help.threading.rangify(values)
Given a list of integers, returns a list of tuples of ranges (interger pairs).
* **Parameters**

**values** – 
* **Returns**

    


### pyrosetta_help.threading.steal_ligands(donor_pose, acceptor_pose)
Steals non-Protein residues from donor_pose and adds them to acceptor_pose

Do not use with nucleic acid polymers.
* **Parameters**

    
 * **donor_pose** – 
 * **acceptor_pose** – 
* **Returns**

    


### pyrosetta_help.threading.get_nonprotein_pose(pose)
pyrosetta.rosetta.protocols.comparative_modeling.StealLigandMover requires some weird things.
This does the same, makes a ligand only pose.
* **Parameters**

**pose** – 
* **Returns**

    


### pyrosetta_help.threading.oligomer_thread(pose, sequence)

### pyrosetta_help.threading.make_fragment_sets(\*poses: Iterable[pyrosetta.rosetta.core.pose.Pose], lengths: Iterable[int] = 3)

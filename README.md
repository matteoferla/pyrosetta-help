# pyrosetta_scripts
Some functions that I keep using over and over as a Python module.

```bash
pip3 install pyrosetta-help
```

Also in GitHub are some Markdown files, with some boilerplates.



## Starting up

A few helper functions.

### get_logger, get_log_entries

The function `configure_logger`, simply adds a stringIO handler to the log and captures the log.
The function `get_log_entries`, spits out entries of a given level.

### make_option_string

This just converts the key:value pairs to a command line string for the pyrosetta init.

* Bools are converted,
* None results in a value argument,
* Tuples are converted to xx:xx type arguments
* Dictionaries are converted to xx:xx type arguments (multiple, if multiple keys in the nested dictionary)

```python
import pyrosetta
from pyrosetta_help import make_option_string, configure_logger, get_log_entries

# capture to log
logger = configure_logger()
# give CLI attributes in a civilised way
pyrosetta.distributed.maybe_init(extra_options=make_option_string(no_optH=False,
                                                ex1=None,
                                                ex2=None,
                                                #mute='all',
                                                ignore_unrecognized_res=True,
                                                load_PDB_components=False,
                                                ignore_waters=False)
                               )
# ...
# show relevant error
print(get_log_entries('ERROR')) 
```  

## Common operations

Import a file, while dealing with the param files
```jupyterpython
from pyrosetta_help.common_ops import pose_from_file
pose = pose_from_file('combined.relaxed2.pdb', params_filenames=['35G.params','CMP.params', 'ATP.params', 'NME.params'])
```
I have somewhere one that via rdkit_to_params starts with a dict of residue:SMILES. TODO find.

Get pandas dataframe of score
```python
from pyrosetta_help import pose2pandas
scores = pose2pandas(pose)
scores.loc[scores.total_score > 10][['residue', 'total_score']]
```
Convert a selector to a list of str of NGL selector style `[resn]resi:chain` 
```jupyterpython
ligand = pyrosetta.rosetta.core.chemical.ResidueProperty.LIGAND
lig_sele = pyrosetta.rosetta.core.select.residue_selector.ResiduePropertySelector(ligand)
clarify_selector(lig_sele, pose)
```
Local relax, etc.
```jupyterpython
ed = prep_ED(pose, map_filename)
local_scorefxn = get_local_scorefxn()
local_relax = get_local_relax()
local_relax.apply(pose)
```

Note, `add_bfactor_from_score` is unstable!
    
## score_mutants

Given a list of mutants and pose, score them. scorefunction terms, interface, movement etc.

```python
from pyrosetta_help import MutantScorer, Mutation, extend_scores
model = MutantScorer(pose, modelname='test')
model.scorefxn = pyrosetta.create_score_function('ref2015')
model.strict_about_starting_residue = True
data = model.score_mutations(['p.Met1Gly', 'p.Ser29Glu'],
                            chain='V',
                            interfaces=(('isolated', 'V_ABCDEFGHIJKLMNOPQRSTWXYZ'),), #
                            preminimise=True,
                            distance=12,
                            cycles=5)
import pandas as pd
scores = pd.DataFrame(data)
extend_scores(scores)
```
    
The function `extend_scores` adds 6 columns, specifying which terms is the biggest changer.
    
## Blueprinter

A key component of using Remodel is a blueprint.
This module makes a blueprint. See doc string of class `Blueprinter` in [blueprint_maker](pyrosetta_help/blueprint_maker/__init__.py) for more.

```python
from pyrosetta_help import Blueprinter
blue = Blueprinter.from_pose(pose)
blue[20:25] = 'NATAA' # wobble
blue.wobble_span(20,25) # same as above.
del blue[15:20] # requires preceeding and suceeding residues to be NATAA though!
blue.del_span(15, 20) # same as above, but wobbles the preceeding and suceeding 1 residues
blue[22] = 'PIKAA W' # requires wobble
blue.mutate(22, 'W') # same as above, but wobbles the preceeding and suceeding residues
```

To set it:

```python
blue.set('mut.blu')
```
This equivalent to the following (handy if something needs manual correction)

```python
blue.write('mut.blu')
blue.bluprint = 'mut.blu'
# which calls `pyrosetta.rosetta.basic.options.set_file_option('remodel:blueprint', 'mut.blu')`
# so do not forget `mover.register_options()`!
```
    
This can therefore be used as normal:
    
```python
pyrosetta.rosetta.basic.options.set_boolean_option('remodel:quick_and_dirty', True)
pyrosetta.rosetta.basic.options.set_string_option('remodel:generic_aa', 'G')

rm = pyrosetta.rosetta.protocols.forge.remodel.RemodelMover()
rm.register_options()
rm.dr_cycles(5) # default 3
rm.max_linear_chainbreak(0.2) # default 0.07
rm.redesign_loop_neighborhood(False)
rm.apply(pose)
```
    
## WeightWatcher

A class to easily look at the weights for the terms of different scorefunctions.
It has various functionalities.

```python
from pyrosetta_help.weights import WeightWatcher
ww = WeightWatcher()
# List available scorefunctions:
print( ww.possible_scorefxn_names ) # dynamic attribute
# Find a scorefunction that mentions a word:
ww.find_metion('foo')
# Get the comment block of a scorefunction:
print(ww.get_scorefxn_comments('ref2015'))
# Get a scorefunction by name by calling the appropriate options first.
# NB. calling a different one will change the options.
scorefxn = ww.get_scorefxn('beta_nov16')
# Get weights dictionary (including ref) for a scorefunction or a name of one
weights = ww.get_weights('beta_nov16')
# Get a pandas table of the weights
weight_table = ww.compare(['ref2015', 'beta_nov16'], different_only=True)
```
It also has the class attribute `term_meanings`,
which is a handy dictionary to convert a score term name into a description.
E.g. converts "fa_atr" -> "Lennard-Jones attractive between atoms 
in different residues (r^6 term, London dispersion forces)." etc.
Taken from Rosetta documentations, with some edits on some terms.
    
## ChainOps

This class works around a list of dict that contain details of each chain in an model. The key keys are:

* number (pose chain number)
* chain (PDB chain letter)
* gene_name (gene name)

the instance can be subscripted with any of those three, returning the dict of that chain.

    from pyrosetta_help import ChainOps

To get a single chain pose:

    chain_pose = chain_ops.get_pose_of_chain(pose, 'B')

Transmogrifer/Murinizer deals with alignments between species.
The RosettaCM stuff is elsewhere. To do: MOVE OVER

The chain details started off by themselves, see [metadata_assembly notes](metadata_assembly.md)

## Other snippets

These could quickly be made into classes... but hey

* [phoshosite plus to pyrosetta](phospho_snippets.md)
* [distance matrix of chains](distances_snippets.md)
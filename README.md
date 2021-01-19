# pyrosetta_scripts
Some scripts that I keep using over and over.

Not quite boilerplates (as they are functions or classes that get imported as opposed to templates), but these are labelled as such in other repos.

## Init

A few helper functions.

### get_logger, get_log_entries

The function `get_logger`, simply adds a stringIO handler to the log and captures the log.
The function `get_log_entries`, spits out entries of a given level.

### make_option_string

This just converts the key:value pairs to a command line string for the pyrosetta init.

* Bools are converted,
* None results in a value argument,
* Tuples are converted to xx:xx type arguments
* Dictionaries are converted to xx:xx type arguments (multiple, if multiple keys in the nested dictionary)


    import pyrosetta
    from init_helper import make_option_string, get_logger
    
    # capture to log
    logger = get_logger()
    # give CLI attributes in a civilised way
    pyrosetta.distributed.maybe_init(extra_options=make_option_string(no_optH=False,
                                                    ex1=None,
                                                    ex2=None,
                                                    mute='all',
                                                    ignore_unrecognized_res=True,
                                                    load_PDB_components=False,
                                                    ignore_waters=False)
                                   )
    ...
    # show relevant error
    print(get_log_entries('ERROR'))   
    
## score_mutants

Given a list of mutants and pose, score them. scorefunction terms, interface, movement etc.

    from score_mutants import MutantScorer
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
    
Additionally, `MutantScorer.term_meanings` is a handy dictionary to convert a score term name into a description.
Taken from Rosetta documentations, with some edits on some terms.
    
## Blueprinter

A key component of using Remodel is a blueprint.
This module makes a blueprint. See doc string of class `Blueprinter` in [blueprint_maker](blueprint_maker/__init__.py) for more.

    from blueprint_maker import Blueprinter
    blue = Blueprinter.from_pose(pose)
    blue[10:14] = 'NATAA' # preceding loop
    del blue[15:20]
    blue[20:25] = 'NATAA' # following loop
    blue[22] = 'PIKAA W'                            
    


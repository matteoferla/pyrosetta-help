## Metadata

The metadata is a list of dictionaries. Each hold data for a chain/peptide.
Example:

    {'name': 'MOUSE MEDIATOR COMPLEX SUBUNIT 27',
     'chain': 'V',
     'module': 'UPPER TAIL',
     'gene_name': 'MED27',
     'recommended_name': 'Mediator of RNA polymerase II transcription subunit 27',
     'human_uniprot': 'Q6P2C8',
     'human_uniprot_name': 'MED27_HUMAN',
     'human_sequence': 'MADVINVSVNLEAFSQAISAIQALRSSVSRVFDCLKDGMRNKETLEGREKAFIAHFQDNLHSVNRDLNELERLSNLVGKPSENHPLHNSGLLSLDPVQDKTPLYSQLLQAYKWSNKLQYHAGLASGLLNQQSLKRSANQMGVSAKRRPKAQPTTLVLPPQYVDDVISRIDRMFPEMSIHLSRPNGTSAMLLVTLGKVLKVIVVMRSLFIDRTIVKGYNENVYTEDGKLDIWSKSNYQVFQKVTDHATTALLHYQLPQMPDVVVRSFMTWLRSYIKLFQAPCQRCGKFLQDGLPPTWRDFRTLEAFHDTCRQ',
     'mouse_uniprot': 'Q9DB40',
     'mouse_uniprot_name': 'MED27_MOUSE',
     'mouse_sequence': 'MADVLSVGVNLEAFSQAISAIQALRSSVSRVFDCLKDGMRNKETLEGREKAFIANFQDNLHSVNRDLNELERLSNLVGKPSENHPLHNSGLLSLDPVQDKTPLYSQLLQAYKWSNKLQYHAGLASGLLNQQSLKRSANQMGVSAKRRPKAQPTTLVLPPQYVDDVISRIDRMFPEMSIHLSRPNGTSAMLLVTLGKVLKVIVVMRSLFIDRTIVKGYNESVYTEDGKLDIWSKSSYQVFQKVTDHATTALLHYQLPQMPDVVVRSFMTWLRSYIKLFQAPCQRCGKFLQDGLPPTWRDFRTLEAFHDTCRQ',
     'number': 21,
     'pose_sequence': 'GVNLEAFSQAISAIQALRSSVSRVFDCLKDGMRNKETLEGREKAFIANFQDNLHSVNRDLNELERLSNLVGKPSENHPLHNSGLLSLDPLYSQLLQAYKWSNKLQYHAGLASGLLNQQSLKRSANPQYVDDVISRIDRMFPEMSIHLSRPNGTSAMLLVTLGKVLKVIVVMRSLFIDRTIVKGYNTEDGKLDIWSKSSYQVFQKVTDHATTALLHYQLPQMPDVVVRSFMTWLRSYIKLFQAPCQRCGKFLQDGLPPTWRDFRTLEA',
     'pose_gap_sequence': '-------GVNLEAFSQAISAIQALRSSVSRVFDCLKDGMRNKETLEGREKAFIANFQDNLHSVNRDLNELERLSNLVGKPSENHPLHNSGLLSLDP------LYSQLLQAYKWSNKLQYHAGLASGLLNQQSLKRSAN--------------------PQYVDDVISRIDRMFPEMSIHLSRPNGTSAMLLVTLGKVLKVIVVMRSLFIDRTIVKGYN----TEDGKLDIWSKSSYQVFQKVTDHATTALLHYQLPQMPDVVVRSFMTWLRSYIKLFQAPCQRCGKFLQDGLPPTWRDFRTLEA',
     'mouse_aln_sequence': 'MADVLSVGVNLEAFSQAISAIQALRSSVSRVFDCLKDGMRNKETLEGREKAFIANFQDNLHSVNRDLNELERLSNLVGKPSENHPLHNSGLLSLDPVQDKTPLYSQLLQAYKWSNKLQYHAGLASGLLNQQSLKRSANQMGVSAKRRPKAQPTTLVLPPQYVDDVISRIDRMFPEMSIHLSRPNGTSAMLLVTLGKVLKVIVVMRSLFIDRTIVKGYNESVYTEDGKLDIWSKSSYQVFQKVTDHATTALLHYQLPQMPDVVVRSFMTWLRSYIKLFQAPCQRCGKFLQDGLPPTWRDFRTLEAFHDTCRQ',
     'human_aln_sequence': 'MADVINVSVNLEAFSQAISAIQALRSSVSRVFDCLKDGMRNKETLEGREKAFIAHFQDNLHSVNRDLNELERLSNLVGKPSENHPLHNSGLLSLDPVQDKTPLYSQLLQAYKWSNKLQYHAGLASGLLNQQSLKRSANQMGVSAKRRPKAQPTTLVLPPQYVDDVISRIDRMFPEMSIHLSRPNGTSAMLLVTLGKVLKVIVVMRSLFIDRTIVKGYNENVYTEDGKLDIWSKSNYQVFQKVTDHATTALLHYQLPQMPDVVVRSFMTWLRSYIKLFQAPCQRCGKFLQDGLPPTWRDFRTLEAFHDTCRQ'}

The keys mean

* `name`: the name from the header
* `chain`: the chain letter (A-Z minus U)
* `module`: Head/Middle/Tail etc.
* `gene_name`: Gene symbol from Uniprot
* `recommended_name`: Protein name form Uniprot
* `human_uniprot`
* `human_uniprot_name`
* `human_sequence`
* `mouse_uniprot`
* `mouse_uniprot_name`
* `mouse_sequence`
* `number`: chain number according to Pyrosetta
* `pose_sequence`: sequence in pose
* `pose_gap_sequence`: sequence in pose with gaps based on PDB residue numbers
* `mouse_aln_sequence`: mouse aligned with human
* `human_aln_sequence`: as above

## start

Using michelanglo data

    from michelanglo_protein import ProteinCore
    
    ProteinCore.settings.verbose = True #False
    ProteinCore.settings.startup(data_folder='/users/brc/matteo/michelanglo/protein-data')
    # as seen in michelanglo_app/views/uniprot_data.py
    
## Get human
    
    import json, os
    human = json.load(open(os.path.join(ProteinCore.settings.dictionary_folder, 'taxid9606-names2uniprot.json')))
    seen = ['Q9UQD0', 'P0DP75', 'Q8WU39', 'Q5VYS4']
    seqs = []
    for gene_name in human:
        if ('MED' in gene_name or gene_name in ('CCNC', 'CDK8', 'CDK19')) and\
            'TMED' not in gene_name and\
            human[gene_name] not in seen:
            p = ProteinCore(uniprot=human[gene_name], taxid=9606).load()
            seen.append(human[gene_name])
            seqs.append({'gene_name': p.gene_name,
                         'recommended_name': p.recommended_name,
                         'human_uniprot': p.uniprot,
                         'human_uniprot_name': p.uniprot_name,
                         'human_sequence': p.sequence
                        })
             
## Get mouse
                        
    seen = []
    mouse = json.load(open(os.path.join(ProteinCore.settings.dictionary_folder, 'taxid10090-names2uniprot.json')))
    sought = [s['gene_name'] for s in seqs]
    for gene_name in mouse:
        if gene_name.replace('_MOUSE','') in sought and not mouse[gene_name] in seen:
            for entry in seqs:
                if gene_name.replace('_MOUSE','') == entry['gene_name']:
                    p = ProteinCore(uniprot=mouse[gene_name], taxid=10090).load()
                    seen.append(mouse[gene_name])
                    entry.update({#'gene_name': p.gene_name,
                                 #'recommended_name': p.recommended_name,
                                 'mouse_uniprot': p.uniprot,
                                 'mouse_uniprot_name': p.uniprot_name,
                                 'mouse_sequence': p.sequence
                                })
                    break
            else:
                raise ValueError(f"Failed to find {gene_name}!!")

## Get header

    import re
    
    chains = []
    chain = {}
    for line in open('Mediator-final-v5.pdb'):
        if 'COMPND' not in line:
            pass
        elif 'MOL_ID' in line:
            if chain:
                chains.append(chain)
                chain = {}
        elif 'MODULE' in line:
            chain['module'] = re.search('MODULE\: (.*)\;', line).group(1)
        elif 'MOLECULE' in line:
            chain['name'] = re.search('MOLECULE\: (.*)\;', line).group(1)
        elif 'CHAIN' in line:
            chain['chain'] = re.search('CHAIN\: (.*)\;', line).group(1)
        else:
            raise ValueError(f'What is {line}')
    else:
        chains.append(chain)
    chains
    
## Merge
merge seq into chains

    seqdex = {re.search('(\d+)', seq['gene_name']).group(1) if re.search('(\d+)', seq['gene_name']) else None: seq for seq in seqs}
    
    for chain in chains:
        number = re.search('(\d+)', chain['name']).group(1)
        chain.update(seqdex[number])
    
    for i, chain in enumerate(chains):
        chain['number'] = i + 1
    chains
    
## Alignment

    from Bio import pairwise2
    from typing import Tuple
    
    def align_seqs(chain: dict) -> Tuple[str, str]:
        """
        in place.
        """
        alignments = pairwise2.align.globalxs(chain['mouse_sequence'],
                                              chain['human_sequence'],
                                              -1, #open
                                              -0.1 #extend
                                             )
        al=alignments[0]
        chain['mouse_aln_sequence'] = al[0]
        chain['human_aln_sequence'] = al[1]
        return al[0], al[1]
    
    # ====================
    for chain in chains:
        align_seqs(chain)
    
## Init pose

    # init
    import pyrosetta
    from init_boilerplate import make_option_string
    pyrosetta.distributed.maybe_init(extra_options=make_option_string(no_optH=False,
                                                    ex1=None,
                                                    ex2=None,
                                                    mute='all',
                                                    ignore_unrecognized_res=True,
                                                    load_PDB_components=False,
                                                    ignore_waters=False)
                                   )
    # pose   
    pose = pyrosetta.Pose()
    filename = 'mediator.fixed.local2.per_chain.pdb'
    pyrosetta.rosetta.core.import_pose.pose_from_file(pose, filename)
    
## Pose sequence

There are two sequences I can get from a pose, the sequence of AA and the sequence based on the PDB numbering

    for chain in chains:
        chain['pose_sequence'] = pose.chain_sequence(chain['number'])
        
PDB form:
namely, YZXAAAXWV to YZXWV gap cannot be fixed my pairwise alignment, so has to be extracted the complicated way.

    # make a dict of chain: dict[resi: resn]
    previous_resi = 0
    previous_chain = None
    from collections import defaultdict
    mapping = defaultdict(dict)
    pdb_info = pose.pdb_info()
    get_resi = lambda res: int(pdb_info.pose2pdb(res).split()[0])
    # chain letter not number like pose.residue(res).chain()
    get_chain = lambda res: pdb_info.pose2pdb(res).split()[1] # str
    for res in range(1, pose.total_residue()+1):
        resi = get_resi(res)
        resn = pose.residue(res).name1()
        chain = get_chain(res) 
        mapping[chain][resi] = resn
        
    # flip
    for chain in mapping:
        seq = ['-'] * max(mapping[chain].keys())
        for resi, resn in mapping[chain].items():
            seq[resi - 1] = resn
        entry = get_entry(chain)
        entry['pose_gap_sequence'] = ''.join(seq)
        
## Write out
write

    import json
    
    with open('metadata.json', 'w') as w:
        json.dump(chains, w)
        
read
        
    import json
    
    with open('metadata.json', 'r') as r:
        chains = json.load(r)
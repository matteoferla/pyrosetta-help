from typing import *
import pyrosetta


class ResInfo:
    pdb_info = None  #: pdb_info: pyrosetta.rosetta.core.pose.PDBInfo

    def __init__(self, number: int, chain: str, icode: str = ' ', segmentID: str = '    '):
        self.number = number
        self.chain = chain
        self.icode = icode
        self.segmentID = segmentID

    @classmethod
    def get(cls, index: int):
        return cls(number=cls.pdb_info.number(index),
                   chain=cls.pdb_info.chain(index),
                   icode=cls.pdb_info.icode(index),
                   segmentID=cls.pdb_info.segmentID(index))

    def set(self,
            index: int,
            pdb_info: Optional[pyrosetta.rosetta.core.pose.PDBInfo] = None):
        if pdb_info is None:
            pdb_info = self.pdb_info
        pdb_info.set_resinfo(res=index,
                             chain_id=self.chain,
                             pdb_res=self.number,
                             ins_code=self.icode,
                             segmentID=self.segmentID)


class BlueCopier:
    def copy_pdb_info(self, original_pose: pyrosetta.Pose, final_pose: pyrosetta.Pose):
        # get original residue info
        ResInfo.pdb_info = original_pose.pdb_info()
        original = []
        previous = None
        for row in self:
            if row[0] == 0:
                # insertion code is a letter
                previous.icode
                raise NotImplementedError
            ri = ResInfo.get(row[0])
            original.append(ri)
            previous = ri
        # set info
        pdb_info = final_pose.pdb_info()
        for i, row in enumerate(original):
            row.set(i + 1, pdb_info)
        pdb_info().obsolete(False)

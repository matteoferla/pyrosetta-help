from typing import Union, List
import pyrosetta


class BlueprinterSubscripted:
    def __iter__(self):
        return filter(len, self.rows)

    def __getitem__(self, idx: Union[int, slice]) -> list:
        if isinstance(idx, int):
            for row in self:
                i = int(row[0])
                if i == idx:
                    return row
        elif isinstance(idx, slice):
            sliced = []
            for row in self:
                i = int(row[0])
                if i >= idx.start and i <= idx.stop:
                    sliced.append(row)
            return sliced
        else:
            raise TypeError(f'Expected int or slice, got a {type(idx)}')
        raise ValueError('index not found.')

    def _setrow(self, row, value):
        idx = int(row[0])
        del row[2:]
        row.append(self.get_ss(idx))
        row.append(str(value))
        return None

    def __setitem__(self, idx: Union[int, slice], value: str) -> None:
        if isinstance(idx, int):
            row = self[idx]
            self._setrow(row, value)
        elif isinstance(idx, slice):
            sliced = []
            for row in self[idx]:
                if '*' in value:
                    # star is a wildcard to mean same as original
                    self._setrow(row, value.replace('*', row[1]))
                else:
                    self._setrow(row, value)
        else:
            raise TypeError(f'Expected int or slice, got a {type(idx)}')

    def __delitem__(self, idx: Union[int, slice]) -> list:
        """
        It is not quite a deletion as the 'row' is still there but empty

        :param idx:
        :return:
        """
        deleted = []
        if isinstance(idx, int):
            row = self[idx]
            deleted.append(row.copy())
            del row[:]
        elif isinstance(idx, slice):
            for row in self[idx]:
                deleted.append(row.copy())
                del row[:]
        else:
            raise TypeError(f'Expected int or slice, got a {type(idx)}')
        return deleted

    @property
    def max(self):
        """
        :return: Highest pose index (this is different from the number of lines due to indels)
        """
        return max([int(row[0]) for row in self])

    def insert(self, idx: int, value: Union[str, List[str]], before: bool = True):
        """

        :param idx:
        :param value:
        :param before: If true, like the list method ``insert`` it inserts before
        :return:
        """
        # homogenise inputs
        if isinstance(value, str):
            values = [value]
        elif isinstance(value, list):
            values = value
        else:
            raise TypeError(f'{type(value)} is not Union[str, List[str]]')
        # insert
        for i, row in enumerate(self):
            if idx == int(row[0]):
                for j, val in enumerate(values):
                    inserendum = [0, 'X', 'D', val]
                    if before:
                        self.rows.insert(i + j, inserendum)
                    else:
                        self.rows.insert(i + j + 1, inserendum)
                break
        else:
            raise ValueError(f'Index {idx} not found.')
        # wobble sides
        if before:
            other_idx = idx - 1
        elif self.max >= idx + 1:
            other_idx = idx + 1
        else:
            other_idx = idx  # does nothing.
        for i in (idx, other_idx):
            if len(self[i]) <= 3:
                self[i] = f'PIKAA {self[1]}'

    def insert_before(self, idx: int, value: Union[str, List[str]]):
        self.insert(idx, value, before=True)

    def insert_after(self, idx: int, value: Union[str, List[str]]):
        self.insert(idx, value, before=False)

    def prepend(self, value: Union[str, List[str]]):
        self.insert(1, value, before=True)

    def append(self, value: Union[str, List[str]]):
        self.insert(self.max, value, before=False)

    # ===== output ==============================

    def __str__(self):
        return '\n'.join([' '.join(map(str, row)) for row in self])

    def to_pdb_str(self, pose):
        pose2pdb = pose.pdb_info().pose2pdb
        lines = []
        for row in self:
            if row[0] == 0:
                lines.append(' '.join(map(str, row)))
            else:
                resi, chain = pose2pdb(row[0]).split()
                lines.append(f'{resi}{chain} '+' '.join(row[1:]))
        return '\n'.join(lines)

    def write(self, bluprint_filename: str):
        with open(bluprint_filename, 'w') as w:
            w.write(str(self))

    def write_pdb_numbered(self, bluprint_filename: str, pose: pyrosetta.Pose):
        with open(bluprint_filename, 'w') as w:
            w.write(self.to_pdb_str(pose))


    def set(self, bluprint_filename: str = 'model.blu'):
        """
        Writes and Sets the blueprint in options
        """
        self.write(bluprint_filename)
        self.blueprint = bluprint_filename

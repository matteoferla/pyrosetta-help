from typing import List

class BlueprinterExpected:
    def expected_seq(self):
        return ''.join([self.get_expected_aa_from_row(row) for row in self])

    def get_expected_aa_from_row(self, row: List[str]) -> str:
        if len(row) == 0: # impossible
            return ''
        elif len(row) == 3 and row[2] == '.': # "30 K ."
            return row[1]
        elif len(row) == 3:  # "30 K H"
            return '*'  # H/L/E/D in lowercase would be ambiguous
        elif len(row) == 4 and row[3] in ('NATAA', 'NATRO'):  # "30 K H NATAA"
            return row[1]
        elif len(row) == 4 and 'PIKAA' in row[3] and len(row[3]) == len('PIKAA ')+1:  # "30 K H PIKAA R"
            return row[3][-1]
        else: # ALLAA etc..
            return '*'







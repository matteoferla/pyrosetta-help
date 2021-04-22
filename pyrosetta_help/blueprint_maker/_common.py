class BlueprinterCommon:
    def pick_native(self, index):
        """
        NATAA and PIKAA native have different results... the latter is better.
        Does nothing if index does not exist.
        """
        # this is not a typo. get item returns the row, setter sets the value (index 3)
        try:
            self[index] = f'PIKAA {self[index][1]}'
        except ValueError:
            pass

    def del_span(self, start: int, stop: int):
        self.pick_native(start - 1)
        del self[start:stop]
        self.pick_native(start + 1)

    def wobble_span(self, start: int, stop: int):
        for resi in range(start,stop+1):
            self.pick_native(resi)

    def mutate(self, resi: int, to_resn: str):
        self.pick_native(resi - 1)
        self[resi] = f'PIKAA {to_resn}'
        self.pick_native(resi + 1)

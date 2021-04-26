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

    def expand_loop_wobble(self):
        raise NotImplementedError
        # identify loops
        loop_resis = [i+1 for i, s in enumerate(self.ss) if s in ('L', 'D')]
        # escape if pointless
        if len(loop_resis) == 0:
            return
        # loops to tuple ranges
        start = loop_resis[0]
        stop = loop_resis[0]
        loop_resi_ranges = []
        for resi in loop_resis:
            if resi != start + 1:
                loop_resi_ranges.append((start, stop))
                start = resi
            stop = resi
        else:
            loop_resi_ranges.append((start, stop))
        # filter ranges that contain mutations
        # [4, 'I', 'H', 'PIKAA I']
        expandenda = []
        for row in self.rows:
            # [4, 'I', 'H', 'PIKAA I']
            if row[2] in ('L', 'D'):
                for start, stop in loop_resi_ranges:
                    if start <= row[0] <= stop:
                        expandenda.append((start, stop))
        # expand ranges
        for start, stop in expandenda:
            for resi in range(start, stop+1):
                try:
                    if self[resi][2] == '.':
                        self.pick_native(resi)
                except:
                    pass


class BlueprinterCommon:
    def del_span(self, start: int, stop: int):
        self[start - 1] = 'NATAA'
        del self[start:stop]
        self[stop + 1] = 'NATAA'

    def wobble_span(self, start: int, stop: int):
        self[start:stop] = 'NATAA'

    def mutate(self, resi: int, to_resn: str):
        self[resi - 1] = 'NATAA'
        self[resi] = f'PIKAA {to_resn}'
        self[resi + 1] = 'NATAA'

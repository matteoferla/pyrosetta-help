import pyrosetta


class Remodel:

    find_neighbors = property(
        lambda self: pyrosetta.rosetta.basic.options.get_boolean_option('remodel:design:find_neighbors'),
        lambda self, value: pyrosetta.rosetta.basic.options.set_boolean_option('remodel:design:find_neighbors', value))


    generic_aa = property(
        lambda self: pyrosetta.rosetta.basic.options.get_string_option('remodel:generic_aa'),
        lambda self, value: pyrosetta.rosetta.basic.options.set_string_option('remodel:generic_aa', value))

    quick_and_dirty = property(
        lambda self: pyrosetta.rosetta.basic.options.get_boolean_option('quick_and_dirty'),
        lambda self, value: pyrosetta.rosetta.basic.options.set_boolean_option('remodel:quick_and_dirty', value))

    blueprint = property(
        lambda self: pyrosetta.rosetta.basic.options.set_file_option('remodel:blueprint'),
        lambda self, filename: pyrosetta.rosetta.basic.options.set_file_option('remodel:blueprint', filename))

    def get_remodelmover(self,
                         dr_cycles:int=3,
                         max_linear_chainbreak:float=0.07,
                         redesign_loop_neighborhood:bool=False):
        rm = pyrosetta.rosetta.protocols.forge.remodel.RemodelMover()
        rm.register_options()
        rm.dr_cycles(dr_cycles)
        rm.max_linear_chainbreak(max_linear_chainbreak)
        rm.redesign_loop_neighborhood(redesign_loop_neighborhood)
        return rm


from enum import Enum

import numpy as np
from ase.data import covalent_radii, atomic_numbers
from pyideallatt.layered.layered import Layered

class LDH(Layered):
    '''
    LDH
    '''
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name_ = 'LayeredDoubleHydroxides'
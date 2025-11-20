from enum import Enum

import numpy as np
from ase.data import covalent_radii, atomic_numbers
from pyideallatt.layered.layered import Layered

class MolybdenumDiselenide(Layered):
    '''
    molybdenum diselenide
    '''
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name_ = 'MolybdenumDiselenide'
        raise NotImplementedError
    
class C3N4(Layered):
    '''
    C3N4
    '''
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name_ = 'C3N4'
        raise NotImplementedError
    

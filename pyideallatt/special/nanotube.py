from typing import Tuple

from pyideallatt.baselatt import BaseLatt
from pyideallatt.layered.layered import Layered

class Nanotube(BaseLatt):
    '''
    Nanotube structure generator
    '''
    def scroll(layered: Layered, chiral_axis: Tuple[int, int]):
        '''
        
        '''
        raise NotImplementedError

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name_ = 'Nanotube'
        raise NotImplementedError
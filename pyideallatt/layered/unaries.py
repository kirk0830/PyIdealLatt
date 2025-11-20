from enum import Enum

import numpy as np
from ase.data import covalent_radii, atomic_numbers
from pyideallatt.layered.layered import Layered

class Graphene(Layered):
    '''graphene'''
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name_ = 'Graphene'

        elem = kwargs['elem']
        if not isinstance(elem, list):
            elem = [elem]
        self.elem_ = elem

        a = kwargs.get('a')
        if a is None:
            # estimate the lattice parameter by the covalent radii
            a = covalent_radii[atomic_numbers[self.elem_[0]]] * np.sqrt(3)
        
        self.a_ = a
        self.b_ = a
        self.c_ = kwargs.get('vacuum', 10) * 2
        self.alpha_ = 90
        self.beta_  = 90
        self.gamma_ = 60

        self.taud_ = np.array([[1/3, 1/3, 0.5],
                               [2/3, 2/3, 0.5]])
    
class LayeredBlackPhosphorus(Layered):
    '''layered black phosphorus'''
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name_ = 'LayeredBlackPhosphorus'
        raise NotImplementedError

class LayeredUnaryType(Enum):

    GRAPHENE, = range(1)

def build(layeredtyp: LayeredUnaryType, **kwargs) -> Layered:
    '''
    build a unary layered structure

    Parameters
    ----------
    layeredtyp
        the type of layered structure
    kwargs
        the keyword arguments for the structure. `elem` must be provided, 
        and `a` the lattice parameter is optional. If `a` is not provided,
        it will be infered from the covalent radii of the element. `vacuum`
        can be provided to set the vacuum size, by default 10 on the both
        sides.

    Returns
    -------
    Layered
        the layered structure
    '''
    builder = [Graphene,]

    return builder[layeredtyp.value](**kwargs)
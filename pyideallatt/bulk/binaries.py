'''
build the ideal binary structure

Usage
-----
```python
from pyideallatt.binary import build, BulkBinaryType

v2o5 = build(BulkBinaryType.X2Y5, x='V', y='O')
print(v2o5)
'''

from pathlib import Path
from enum import Enum
from typing import TypeAlias

import numpy as np
from ase.data import covalent_radii, atomic_numbers

from pyideallatt.bulk.bulk import Bulk

class IdealX2O(Bulk):
    '''
    Ideal X2O structure, for more information, see
    Bosoni E, Beal L, Bercx M, et al. 
    How to verify the precision of density-functional-theory 
    implementations via reproducible and universal workflows[J]. 
    Nature Reviews Physics, 2024, 6(1): 45-58.
    '''
    DATAPATH_ = Path(__file__).parent / 'data' / 'x2o.yaml'

    def find_a_from_volume(v: float) -> float:
        '''find the a parameter from the volume of primitive cell'''
        return np.sqrt(2)/2 * (4*v)**(1/3)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name_ = 'IdealX2O'

        x = kwargs['x']
        assert isinstance(x, str)
        y = kwargs.get('y', 'O')
        assert isinstance(y, str)

        self.elem_ = [x, x, y]

        a = kwargs.get('a')
        if a is None:
            v = self.get_lattice_data(x)
            a = IdealX2O.find_a_from_volume(v)
            # rescale the lattice by the ratio of 
            # covalent radii of y and O
            if y != 'O':
                a *= covalent_radii[atomic_numbers[y]] / covalent_radii[8]
        self.a_ = a
        self.b_ = a
        self.c_ = a

        self.alpha_ = 60
        self.beta_  = 60
        self.gamma_ = 60

        self.taud_ = np.array([[0.25000, 0.25000, 0.25000],
                               [0.75000, 0.75000, 0.75000],
                               [0.00000, 0.00000, 0.00000]])

X2O: TypeAlias = IdealX2O

class IdealXO(Bulk):
    '''
    Ideal XO structure, for more information, see
    Bosoni E, Beal L, Bercx M, et al. 
    How to verify the precision of density-functional-theory 
    implementations via reproducible and universal workflows[J]. 
    Nature Reviews Physics, 2024, 6(1): 45-58.
    '''
    DATAPATH_ = Path(__file__).parent / 'data' / 'xo.yaml'

    def find_a_from_volume(v: float) -> float:
        '''find the a parameter from the volume of primitive cell'''
        return np.sqrt(2)/2 * (4*v)**(1/3)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name_ = 'IdealXO'

        x = kwargs['x']
        assert isinstance(x, str)
        y = kwargs.get('y', 'O')
        assert isinstance(y, str)
        self.elem_ = [x, y]

        a = kwargs.get('a')
        if a is None:
            v = self.get_lattice_data(x)
            a = IdealX2O.find_a_from_volume(v)
            # rescale the lattice by the ratio of 
            # covalent radii of y and O
            if y != 'O':
                a *= covalent_radii[atomic_numbers[y]] / covalent_radii[8]
        self.a_ = a
        self.b_ = a
        self.c_ = a

        self.alpha_ = 60
        self.beta_  = 60
        self.gamma_ = 60

        self.taud_ = np.array([[0.00000, 0.00000, 0.00000],
                               [0.50000, 0.50000, 0.50000]])

XO: TypeAlias = IdealXO

class IdealX2O3(Bulk):
    '''
    Ideal X2O3 structure, for more information, see
    Bosoni E, Beal L, Bercx M, et al. 
    How to verify the precision of density-functional-theory 
    implementations via reproducible and universal workflows[J]. 
    Nature Reviews Physics, 2024, 6(1): 45-58.
    '''
    DATAPATH_ = Path(__file__).parent / 'data' / 'x2o3.yaml'

    def find_a_from_volume(v: float) -> float:
        '''find the a parameter from the volume of primitive cell'''
        return (v)**(1/3)
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name_ = 'IdealX2O3'

        x = kwargs['x']
        assert isinstance(x, str)
        y = kwargs.get('y', 'O')
        assert isinstance(y, str)
        self.elem_ = [x, x, x, x, y, y, y, y, y, y]

        a = kwargs.get('a')
        if a is None:
            v = self.get_lattice_data(x)
            a = IdealX2O3.find_a_from_volume(v)
            # rescale the lattice by the ratio of 
            # covalent radii of y and O
            if y != 'O':
                a *= covalent_radii[atomic_numbers[y]] / covalent_radii[8]
        self.a_ = a
        self.b_ = a
        self.c_ = a

        self.alpha_ = 90
        self.beta_  = 90
        self.gamma_ = 90

        self.taud_ = np.array([[0.25000, 0.25000, 0.25000],
                               [0.75000, 0.75000, 0.25000],
                               [0.75000, 0.25000, 0.75000],
                               [0.25000, 0.75000, 0.75000],
                               [0.50000, 0.00000, 0.00000],
                               [0.00000, 0.50000, 0.00000],
                               [0.00000, 0.00000, 0.50000],
                               [0.50000, 0.50000, 0.00000],
                               [0.50000, 0.00000, 0.50000],
                               [0.00000, 0.50000, 0.50000]])

X2O3: TypeAlias = IdealX2O3

class IdealXO2(Bulk):
    '''
    Ideal XO2 structure, for more information, see
    Bosoni E, Beal L, Bercx M, et al. 
    How to verify the precision of density-functional-theory 
    implementations via reproducible and universal workflows[J]. 
    Nature Reviews Physics, 2024, 6(1): 45-58.
    '''
    DATAPATH_ = Path(__file__).parent / 'data' / 'xo2.yaml'

    def find_a_from_volume(v: float) -> float:
        '''find the a parameter from the volume of primitive cell'''
        return np.sqrt(2)/2 * (4*v)**(1/3)
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name_ = 'IdealXO2'

        x = kwargs['x']
        assert isinstance(x, str)
        y = kwargs.get('y', 'O')
        assert isinstance(y, str)
        self.elem_ = [x, y, y]

        a = kwargs.get('a')
        if a is None:
            v = self.get_lattice_data(x)
            a = IdealXO2.find_a_from_volume(v)
            # rescale the lattice by the ratio of 
            # covalent radii of y and O
            if y != 'O':
                a *= covalent_radii[atomic_numbers[y]] / covalent_radii[8]
        self.a_ = a
        self.b_ = a
        self.c_ = a

        self.alpha_ = 60
        self.beta_  = 60
        self.gamma_ = 60

        self.taud_ = np.array([[0.00000, 0.00000, 0.00000],
                               [0.25000, 0.25000, 0.25000],
                               [0.75000, 0.75000, 0.75000]])

XO2: TypeAlias = IdealXO2

class IdealX2O5(Bulk):
    '''
    Ideal X2O5 structure, for more information, see
    Bosoni E, Beal L, Bercx M, et al. 
    How to verify the precision of density-functional-theory 
    implementations via reproducible and universal workflows[J]. 
    Nature Reviews Physics, 2024, 6(1): 45-58.
    '''
    DATAPATH_ = Path(__file__).parent / 'data' / 'x2o5.yaml'

    def find_a_from_volume(v: float) -> float:
        '''find the a parameter from the volume of primitive cell'''
        return (v)**(1/3)
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name_ = 'IdealX2O5'

        x = kwargs['x']
        assert isinstance(x, str)
        y = kwargs.get('y', 'O')
        assert isinstance(y, str)
        self.elem_ = [x, x, x, x, y, y, y, y, y, y, y, y, y, y]

        a = kwargs.get('a')
        if a is None:
            v = self.get_lattice_data(x)
            a = IdealX2O5.find_a_from_volume(v)
            # rescale the lattice by the ratio of 
            # covalent radii of y and O
            if y != 'O':
                a *= covalent_radii[atomic_numbers[y]] / covalent_radii[8]
        self.a_ = a
        self.b_ = a
        self.c_ = a

        self.alpha_ = 90
        self.beta_  = 90
        self.gamma_ = 90

        self.taud_ = np.array([[0.00000, 0.00000, 0.00000],
                               [0.50000, 0.50000, 0.00000],
                               [0.50000, 0.00000, 0.50000],
                               [0.00000, 0.50000, 0.50000],
                               [0.50000, 0.00000, 0.00000],
                               [0.00000, 0.50000, 0.00000],
                               [0.00000, 0.00000, 0.50000],
                               [0.50000, 0.50000, 0.50000],
                               [0.75000, 0.75000, 0.25000],
                               [0.75000, 0.25000, 0.75000],
                               [0.25000, 0.75000, 0.75000],
                               [0.75000, 0.25000, 0.25000],
                               [0.25000, 0.75000, 0.25000],
                               [0.25000, 0.25000, 0.75000]])

X2O5: TypeAlias = IdealX2O5

class IdealXO3(Bulk):
    '''
    Ideal XO3 structure, for more information, see
    Bosoni E, Beal L, Bercx M, et al. 
    How to verify the precision of density-functional-theory 
    implementations via reproducible and universal workflows[J]. 
    Nature Reviews Physics, 2024, 6(1): 45-58.
    '''
    DATAPATH_ = Path(__file__).parent / 'data' / 'xo3.yaml'

    def find_a_from_volume(v: float) -> float:
        '''find the a parameter from the volume of primitive cell'''
        return (v)**(1/3)
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name_ = 'IdealXO3'

        x = kwargs['x']
        assert isinstance(x, str)
        y = kwargs.get('y', 'O')
        assert isinstance(y, str)
        self.elem_ = [x, y, y, y]

        a = kwargs.get('a')
        if a is None:
            v = self.get_lattice_data(x)
            a = IdealXO3.find_a_from_volume(v)
            # rescale the lattice by the ratio of 
            # covalent radii of y and O
            if y != 'O':
                a *= covalent_radii[atomic_numbers[y]] / covalent_radii[8]
        self.a_ = a
        self.b_ = a
        self.c_ = a

        self.alpha_ = 90
        self.beta_  = 90
        self.gamma_ = 90

        self.taud_ = np.array([[0.00000, 0.00000, 0.00000],
                               [0.50000, 0.00000, 0.00000],
                               [0.00000, 0.50000, 0.00000],
                               [0.00000, 0.00000, 0.50000]])

XO3: TypeAlias = IdealXO3

class BulkBinaryType(Enum):

    X2Y, XY, X2Y3, XY2, X2Y5, XY3 = range(6)
    
def build(bulktyp: BulkBinaryType, **kwargs) -> Bulk:
    '''
    build a binary bulk structure

    Parameters
    ----------
    bulktyp: BulkBinaryType
        the type of binary bulk structure
    kwargs: dict
        the keyword arguments for the bulk structure. For `X2Y`, `XY`, 
        `X2Y3`, `X2O3`, `XO2`, `X2O5`, `XO3`, `x` must be defined for
        specifying the element at the x-site, `y` is O by default. 
        `a` is the lattice parameter of the primitive cell, if not defined,
        an estimated value will be used.
    '''
    builder = [IdealX2O,
               IdealXO,
               IdealX2O3,
               IdealXO2,
               IdealX2O5,
               IdealXO3]
    
    return builder[bulktyp.value](**kwargs)
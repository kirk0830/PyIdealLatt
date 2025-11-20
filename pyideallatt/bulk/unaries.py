'''
build the ideal unary bulk structures

Usage
-----
```python
from pyideallatt.unaries import build, BulkUnaryType

copper = build(BulkUnaryType.FACECENTEREDCUBIC, elem='Cu', a=3.6)
print(copper)

# if the `a` is not provided, will infered from the covalent radii
# of the element
copper = build(BulkUnaryType.FACECENTEREDCUBIC, elem='Cu')
print(copper)
```
'''
from pathlib import Path
from enum import Enum
from typing import TypeAlias

import numpy as np
from ase.data import covalent_radii, atomic_numbers
from pyideallatt.bulk.bulk import Bulk

class SimpleCubic(Bulk):
    '''simple-cubic'''
    DATAPATH_ = Path(__file__).parent / 'data' / 'sc.yaml'

    def find_a_from_volume(v: float) -> float:
        '''find the a parameter from the volume of primitive cell'''
        return (v)**(1/3)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name_ = 'SimpleCubic'

        elem = kwargs['elem']
        if not isinstance(elem, list):
            elem = [elem]

        assert len(elem) == 1
        assert all(isinstance(e, str) for e in elem)
        self.elem_ = elem

        a = kwargs.get('a')
        if a is None:
            v = self.get_lattice_data(elem)
            a = SimpleCubic.find_a_from_volume(v)
        self.a_ = a
        self.b_ = a
        self.c_ = a
        self.alpha_ = 90
        self.beta_  = 90
        self.gamma_ = 90
        
        self.taud_ = np.array([[0.5, 0.5, 0.5]])

SC: TypeAlias = SimpleCubic 

class BodyCenteredCubic(Bulk):
    '''body-centered cubic'''
    DATAPATH_ = Path(__file__).parent / 'data' / 'bcc.yaml'

    def find_a_from_volume(v: float) -> float:
        '''find the a parameter from the volume of primitive cell'''
        return np.sqrt(3)/2 * (2*v)**(1/3)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name_ = 'BodyCenteredCubic'

        elem = kwargs['elem']
        if not isinstance(elem, list):
            elem = [elem]

        assert len(elem) == 1
        assert all(isinstance(e, str) for e in elem)
        self.elem_ = elem

        a = kwargs.get('a')
        if a is None:
            v = self.get_lattice_data(elem)
            a = BodyCenteredCubic.find_a_from_volume(v)
        self.a_ = a
        self.b_ = a
        self.c_ = a

        self.alpha_ = 109.4712
        self.beta_  = 109.4712
        self.gamma_ = 109.4712

        self.taud_ = np.array([[0, 0, 0]])

BCC: TypeAlias = BodyCenteredCubic

class FaceCenteredCubic(Bulk):
    '''face-centered cubic'''
    DATAPATH_ = Path(__file__).parent / 'data' / 'fcc.yaml'

    def find_a_from_volume(v: float) -> float:
        '''find the a parameter from the volume of primitive cell'''
        return np.sqrt(2)/2 * (4*v)**(1/3)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name_ = 'FaceCenteredCubic'

        elem = kwargs['elem']
        if not isinstance(elem, list):
            elem = [elem]

        assert len(elem) == 1
        assert all(isinstance(e, str) for e in elem)
        self.elem_ = elem

        a = kwargs.get('a')
        if a is None:
            a = 4 * covalent_radii[atomic_numbers[elem]] / np.sqrt(2) / 2
        self.a_ = a
        self.b_ = a
        self.c_ = a

        self.alpha_ = 60
        self.beta_  = 60
        self.gamma_ = 60

        self.taud_ = np.array([[0, 0, 0]])

FCC: TypeAlias = FaceCenteredCubic

def calculate_equilateral_triangle_basedge(a, theta, rad=False):
    '''calculate the edge length of a equal-lateral-triangle'''
    theta = np.deg2rad(theta) if not rad else theta
    return np.sqrt(2 * a**2 * (1 - np.cos(theta)))

class Diamond(Bulk):
    '''diamond'''
    DATAPATH_ = Path(__file__).parent / 'data' / 'diamond.yaml'

    def find_a_from_volume(v: float) -> float:
        '''find the a parameter from the volume of primitive cell'''
        return np.sqrt(2)/2 * (4*v)**(1/3)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name_ = 'Diamond'

        elem = kwargs['elem']
        if not isinstance(elem, list):
            elem = [elem for _ in range(2)]

        assert len(elem) == 2
        assert all(isinstance(e, str) for e in elem)
        self.elem_ = elem

        a = kwargs.get('a')
        if a is None:
            tmp = covalent_radii[atomic_numbers[elem]] * 2
            a = calculate_equilateral_triangle_basedge(tmp, 109.4712)
        self.a_ = a
        self.b_ = a
        self.c_ = a

        self.alpha_ = 60
        self.beta_  = 60
        self.gamma_ = 60

        self.taud_ = np.array([[0.25, 0.25, 0.25], 
                               [0.75, 0.75, 0.75]])

class BetaTin(Bulk):
    '''beta-tin'''
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name_ = 'BetaTin'
        
        elem = kwargs['elem']
        if not isinstance(elem, list):
            elem = [elem for _ in range(2)]

        assert len(elem) == 2
        assert all(isinstance(e, str) for e in elem)
        self.elem_ = elem

        a = kwargs.get('a')
        if a is None:
            a = calculate_equilateral_triangle_basedge(
                a=2*covalent_radii[atomic_numbers[elem]], 
                theta=94.0229)
        
        self.a_ = a
        self.b_ = a
        self.c_ = a
        self.alpha_ = 97.53393241
        self.beta_  = 97.53393241
        self.gamma_ = 137.54208131

        self.taud_ = np.array([[0.  , 0.  , 0.  ], 
                               [0.25, 0.75, 0.50]])
        
class BulkUnaryType(Enum):
    '''
    The enumerator for the unary bulk structures.
    '''
    SIMPLECUBIC, \
    BODYCENTEREDCUBIC, \
    FACECENTEREDCUBIC, \
    DIAMOND, \
    BETATIN = range(5)

def build(bulktyp: BulkUnaryType, **kwargs) -> Bulk:
    '''
    build a unary bulk structure

    Parameters
    ----------
    bulktyp
        the type of the bulk structure
    kwargs
        the arguments for the bulk structure. For `SIMPLECUBIC`,
        `BODYCENTEREDCUBIC`, `FACECENTEREDCUBIC`, `DIAMOND` and `BETATIN`,
        `elem` must be provided. `a` is the lattice parameter, if not provided,
        an estimated value will be used.

    Returns
    -------
    bulk
        the bulk structure
    '''
    builder = [SimpleCubic, 
               BodyCenteredCubic, 
               FaceCenteredCubic, 
               Diamond, 
               BetaTin]
    
    return builder[bulktyp.value](**kwargs)


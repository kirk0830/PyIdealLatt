import yaml
import numpy as np
from pathlib import Path
from ase.data import covalent_radii, atomic_numbers
from pyideallatt.molecule.molecule import Molecule

class Dimer(Molecule):
    '''
    0D molecule: the homonuclear diatomic molecule
    '''
    DATA_ = None
    DATAPATH_ = Path(__file__).parent / 'data' / 'dimer.yaml'

    @classmethod
    def load_lattice_data(cls):
        '''load the lattice data from the data file'''
        if cls.DATAPATH_ is None:
            raise RuntimeError(f'{cls.__name__} does not have a '
                'configuration on the lattice data file. '
                'Please either change the implementation of __init__ '
                'or configure a correct DATAPATH_ attribute.')
        with open(cls.DATAPATH_, 'r') as f:
            cls.DATA_ = yaml.safe_load(f)

    @classmethod
    def get_lattice_data(cls, elem: str) -> float:
        '''look up the lattice data from the data file'''
        if cls.DATA_ is None:
            cls.load_lattice_data()
        return cls.DATA_[elem]
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name_ = 'Dimer'

        elem = kwargs['elem']
        if not isinstance(elem, list):
            elem = [elem, elem]

        assert len(elem) == 2
        assert all(isinstance(e, str) for e in elem)
        self.elem_ = elem

        a = kwargs.get('a')
        if a is not None:
            self.a_ = a
            self.b_ = a
            self.c_ = a

        bl = kwargs.get('bl')
        if bl is None:
            bl = self.get_lattice_data(elem[0])
            if len(set(elem)) > 1:
                bl = bl / 2 * (1 + covalent_radii[atomic_numbers[elem[0]]] / covalent_radii[atomic_numbers[elem[1]]])
        self.bl_ = bl

        tauc_ = np.array([[- bl/2., 0., 0.],
                          [+ bl/2., 0., 0.]]) \
              + np.array([[self.a_/2, self.b_/2, self.c_/2]])
        cell_ = np.array([[self.a_, 0., 0.],
                          [0., self.b_, 0.],
                          [0., 0., self.c_]])
        
        # cell @ taud_.T = tauc_.T
        self.taud_ = np.linalg.solve(cell_, tauc_.T).T

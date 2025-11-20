import yaml

import numpy as np
from scipy.optimize import minimize
from ase.filters import StrainFilter, FrechetCellFilter
from ase.optimize import QuasiNewton

from pyideallatt.baselatt import BaseLatt

class Bulk(BaseLatt):
    '''
    3D periodic structure
    '''
    DATA_ = None
    DATAPATH_ = None

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
        self.name_ = 'Bulk'

    def relax(self, calculator, **kwargs) -> None:
        '''
        relax the lattice and atoms' positions

        Parameters
        ----------
        calculator
            the calculator that can evaluate the energy/forces/stress. It 
            should be the derived class of 
            :class:`~ase.calculators.calculator.Calculator`.
        kwargs
            additional arguments for the optimizer, see 
            :class:`~ase.optimize.optimize.Optimizer` for details. There
            are mainly 4 usefull arguments: `fmax`, `logfile`, `trajectory`
            and `steps`.
        '''
        atoms = self.toase()
        atoms.calc = calculator
        
        ucf = FrechetCellFilter(atoms)
        dyn = QuasiNewton(ucf, **kwargs)
        
        options = {k: kwargs.get(k) for k in ['fmax', 'steps']}
        options = {k: v for k, v in options.items() if v is not None}
        dyn.run(**options)

        # update the structure
        cellpar = atoms.cell.cellpar()
        
        self.a_     = cellpar[0]
        self.b_     = cellpar[1]
        self.c_     = cellpar[2]
        self.alpha_ = cellpar[3]
        self.beta_  = cellpar[4]
        self.gamma_ = cellpar[5]

        self.taud_ = atoms.get_scaled_positions()

    def relax_cell_edge(self, calculator, **kwargs) -> None:
        '''
        relax solely the length of lattice vectors, keeping the scaled
        position of atoms and lattice angles unchanged.

        Parameters
        ----------
        calculator
            the calculator that can evaluate the energy/forces/stress. It 
            should be the derived class of 
            :class:`~ase.calculators.calculator.Calculator`.
        kwargs
            additional arguments for the optimizer, see 
            :class:`~ase.optimize.optimize.Optimizer` for details. There
            are mainly 4 usefull arguments: `fmax`, `logfile`, `trajectory`
            and `steps`.
        '''
        atoms = self.toase()
        atoms.calc = calculator

        sf = StrainFilter(atoms, mask=[1, 1, 1, 0, 0, 0])
        dyn = QuasiNewton(sf, **kwargs)

        options = {k: kwargs.get(k) for k in ['fmax', 'steps']}
        options = {k: v for k, v in options.items() if v is not None}
        dyn.run(**options)

        # update the structure
        cellpar = atoms.cell.cellpar()
        
        self.a_ = cellpar[0]
        self.b_ = cellpar[1]
        self.c_ = cellpar[2]

    def relax_cell_scale(self, calculator, **kwargs):
        '''
        relax the three lengths of lattice vectors with the same ratio,
        keeping the scaled position of atoms and lattice angles unchanged.

        Parameters
        ----------
        calculator
            the calculator that can evaluate the energy. It 
            should be the derived class of 
            :class:`~ase.calculators.calculator.Calculator`.
        kwargs
            additional arguments for the optimizer, see 
            :class:`~scipy.optimize.minimize` for details.
        '''
        def f(eps):
            atoms = self.toase()
            cell = atoms.get_cell()
            atoms.set_cell(cell * (1 + eps), scale_atoms=True)
            atoms.calc = calculator
            return atoms.get_potential_energy()
        
        res = minimize(f, 
                       x0=0.0, 
                       method='L-BFGS-B', 
                       bounds=[(-1, np.inf)],
                       tol=1e-6)
        eps = res.x[0]

        self.a_ *= (1 + eps)
        self.b_ *= (1 + eps)
        self.c_ *= (1 + eps)

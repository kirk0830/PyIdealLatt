import numpy as np
from scipy.optimize import minimize
from ase.filters import StrainFilter, FrechetCellFilter
from ase.optimize import QuasiNewton
from ase.geometry import cellpar_to_cell

from pyideallatt.baselatt import BaseLatt

class Layered(BaseLatt):
    '''
    2D periodic structure (XY-plane)
    '''
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name_ = 'Layered'

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
        
        ucf = FrechetCellFilter(atoms, mask=[1, 1, 0, 1, 1, 1])
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

        sf = StrainFilter(atoms, mask=[1, 1, 0, 0, 0, 0])
        dyn = QuasiNewton(sf, **kwargs)

        options = {k: kwargs.get(k) for k in ['fmax', 'steps']}
        options = {k: v for k, v in options.items() if v is not None}
        dyn.run(**options)

        # update the structure
        cellpar = atoms.cell.cellpar()
        
        self.a_ = cellpar[0]
        self.b_ = cellpar[1]

    def relax_cell_scale(self, calculator, **kwargs) -> None:
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
            cellpar = atoms.cell.cellpar()
            cellpar[0] *= (1 + eps)
            cellpar[1] *= (1 + eps)
            atoms.set_cell(cellpar_to_cell(cellpar), scale_atoms=True)
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
from abc import ABC, abstractmethod
from typing import List, Dict

from ase import Atoms
from ase.optimize import BFGS
from ase.geometry import cellpar_to_cell

class BaseLatt(ABC):
    '''
    The base class for all structures. This class provides a template
    for implementing all structures
    '''
    def __init__(self, **kwargs):
        self.name_ = 'BaseLatt'
        
        # cell
        self.a_ = None
        self.b_ = None
        self.c_ = None
        self.alpha_ = None
        self.beta_ = None
        self.gamma_ = None
        
        # atoms (direct coordinates)
        self.elem_ = None
        self.taud_ = None

    def todict(self) -> Dict:
        return {
            'name': self.name_,
            'a': self.a_,
            'b': self.b_,
            'c': self.c_,
            'alpha': self.alpha_,
            'beta': self.beta_,
            'gamma': self.gamma_,
            'elem': self.elem_,
            'taud': self.taud_
        }

    @classmethod
    def fromdict(cls, mydict: Dict) -> 'BaseLatt':
        return cls(**mydict)

    def copy(self) -> 'BaseLatt':
        return self.fromdict(self.todict())

    def __repr__(self) -> str:
        mydict = self.todict()
        mydict['elem'] = ''.join(mydict['elem'])
        mydict['taud'] = '...'
        return f'<{",".join([f"{k}={v}" for k, v in mydict.items()])}>'
    
    def __str__(self) -> str:
        return self.__repr__()

    def valid(self) -> bool:
        ''' check if the structure is valid '''
        return all(x is not None for x in self.todict().values())

    def cellpar(self) -> List[float]:
        ''' get the cell parameters '''
        return [self.a_, self.b_, self.c_, self.alpha_, self.beta_, self.gamma_]

    def toase(self) -> Atoms:
        ''' export structures to the ase.atoms.Atoms '''
        return Atoms(symbols=self.elem_,
                     scaled_positions=self.taud_,
                     cell=cellpar_to_cell(self.cellpar()))
    
    @classmethod
    def fromase(cls, atoms: Atoms, **kwargs) -> 'BaseLatt':
        ''' import structures from the ase.atoms.Atoms '''
        cellpar = atoms.cell.cellpar()
        elem = atoms.get_chemical_symbols()
        taud = atoms.get_scaled_positions()
        return cls(a=cellpar[0], b=cellpar[1], c=cellpar[2],
                   alpha=cellpar[3], beta=cellpar[4], gamma=cellpar[5],
                   elem=elem, taud=taud, **kwargs)

    @abstractmethod
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
        # for the implementation of this method would be different for
        # 3D, 2D, 1D and 0D structures, this method is not implemented
        # here.
        raise NotImplementedError

    @abstractmethod
    def relax_cell_edge(self, calculator, **kwargs) -> None:
        '''
        relax solely the lattice constants

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
        raise NotImplementedError

    def relax_atom(self, calculator, **kwargs) -> None:
        '''
        relax solely the atoms' positions 

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
        dyn = BFGS(atoms, **kwargs)
        dyn.run(fmax=kwargs.get('fmax', 0.05), 
                steps=kwargs.get('steps', 1000))
        
        # update the atoms
        self.taud_ = atoms.get_scaled_positions()

    @abstractmethod
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
        raise NotImplementedError
    
    @abstractmethod
    def scale(self, eps, bycopy=False):
        '''
        scale the lattice vectors with the same ratio. Because the `scale`
        behaves differently for different structures (especially for 
        different dimensions), this method is not implemented here.

        Parameters
        ----------
        eps
            the ratio of scaling.
        bycopy
            if True, return a new structure, otherwise, modify the current
            structure.
        
        Returns
        -------
        None | BaseLatt
            if `bycopy` is True, return a new structure, otherwise, return
            None, the current structure is modified.
        '''
        raise NotImplementedError
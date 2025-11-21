from pyideallatt.baselatt import BaseLatt

class Molecule(BaseLatt):
    '''
    0D molecule
    '''
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name_ = 'Molecule'

        a = kwargs.get('a')
        if a is None:
            a = 20.
        self.a_ = a
        self.b_ = a
        self.c_ = a

        self.alpha_ = 90.
        self.beta_  = 90.
        self.gamma_ = 90.

        # self.taud_ would be defined in the derived class

    def relax(self, calculator, **kwargs) -> None:
        '''
        relax the atoms' positions

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
        return self.relax_atom(calculator, **kwargs)

    def relax_cell_edge(self, calculator, **kwargs):
        ''' Any cell parameters should not be changed '''
        raise RuntimeError('Any cell parameters should not be changed for molecule')
    
    def relax_cell_scale(self, calculator, **kwargs):
        ''' Any cell parameters should not be changed '''
        raise RuntimeError('Any cell parameters should not be changed for molecule')
    
    def scale(self, eps: float, bycopy: bool = False) -> 'Molecule':
        '''
        scale the bond length of the molecule

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
        toscale = self.copy() if bycopy else self
        
        # we rescale the bonds by the ratio of eps
        center = toscale.taud_.mean(axis=0)
        toscale.taud_ = (toscale.taud_ - center) * eps + center

        return toscale if bycopy else None
    

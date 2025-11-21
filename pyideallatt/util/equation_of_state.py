'''
Usage
-----
```python
from pathlib import Path
here = Path(__file__).parent
from ase.calculators.espresso import Espresso, EspressoProfile
from mattersim.forcefield import MatterSimCalculator
from pyideallatt.util.equation_of_state import calculate, delta
from pyideallatt.bulk.unaries import build, BulkUnaryType

# calculate the delta value between the DFT code and Machine-learning
# forcefield...

copper = build(BulkUnaryType.FACECENTEREDCUBIC, elem='Cu')

profile = EspressoProfile(
        command='mpirun -np 16 pw.x',
        pseudo_dir=here
)
dft_calculator = Espresso(
    profile=profile,
    pseudopotentials={'Cu': 'Cu.paw.z_11.ld1.psl.v1.0.0-low.upf'},
    kpts=(11, 11, 11),
    input_data={
        'system': {
            'ecutwfc': 50,
            'ecutrho': 200
        },
        'tstress': True,
        'tprnfor': True,
        'occupations': 'smearing',
        'degauss': 0.01,
        'smearing': 'cold',
        'conv_thr': 1e-8,
        'mixing_mode': 'local-TF',
        'mixing_beta': 0.35,
        'diagonalization': 'david',
        'startingwfc': 'random'
        }
)

mlp_calculator = MatterSimCalculator(
    # default is the mattersim-v1.0.0-1M.pth
)

# calculate the EOS predicted by DFT
copper_dft = copper.copy()
copper_dft.relax(calculator=dft_calculator)
eos_dft = calculate(copper_dft,
                    calculator=dft_calculator)

# calculate the EOS predicted by MLP
copper_mlp = copper.copy()
copper_mlp.relax(calculator=mlp_calculator)
eos_mlp = calculate(copper_mlp,
                    calculator=mlp_calculator)

# calculate the delta value
v1, b1, bp1 = eos_dft[3], eos_dft[1], eos_dft[2]
v2, b2, bp2 = eos_mlp[3], eos_mlp[1], eos_mlp[2]
d = delta(v1, b1, bp1, v2, b2, bp2) / len(copper.elem_) * 1e3
print(f'the delta value is {d:.2f} meV/atom')
```
'''

from typing import Tuple

import numpy as np
from scipy.optimize import curve_fit
from ase.utils.deltacodesdft import delta
from ase.calculators.calculator import Calculator

from pyideallatt.baselatt import BaseLatt

def birchmurnaghan(v, e0, b0, b0p, v0):
    return e0 + 9 * v0 * b0 / 16 * (b0p * ((v0 / v)**(2/3) - 1)**3 + ((v0 / v)**(2/3) - 1)**2 * (6 - 4 * (v0 / v)**(2/3)))

def calculate(latt: BaseLatt, 
              lo: float = 0.96,
              hi: float = 1.04, 
              n: int = 7,
              calculator: Calculator = None) -> Tuple[float, float, float, float]:
    '''
    calculate the equation of state on the given lattice

    Parameters
    ----------
    latt
        the lattice to calculate the equation of state
    lo
        the lower bound of the volume
    hi
        the upper bound of the volume
    n
        the number of volumes to calculate
    calculator
        the calculator to calculate the energy
        
    Returns
    -------
    Tuple[float, float, float, float]
        e0, b0, b0p and v0
    '''
    assert calculator is not None
    def f(eps: float) -> Tuple[float, float]:
        '''get the volume and energy under present scaling'''
        scaled = latt.scale(eps, bycopy=True).toase()
        scaled.calc = calculator
        return scaled.get_volume(), scaled.get_potential_energy()
    
    e, v = zip(*[f(eps**(1/3)) for eps in np.linspace(lo, hi, n)])
    
    imin = np.argmin(e) # for the initial guess of the e0 and v0
    popt, _ = curve_fit(birchmurnaghan, v, e, p0=(e[imin], 1, 1, v[imin]))

    e0, b0, b0p, v0 = popt
    return e0, b0, b0p, v0


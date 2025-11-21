from pathlib import Path
here = Path(__file__).parent
from ase.calculators.espresso import Espresso, EspressoProfile
from mattersim.forcefield import MatterSimCalculator
from pyideallatt.util.equation_of_state import calculate, delta
from pyideallatt.bulk.unaries import build, BulkUnaryType

# calculate the delta value between the DFT code and Machine-learning
# forcefield...

copper = build(BulkUnaryType.FACECENTEREDCUBIC, elem='Cu')

qeparam = {
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
profile = EspressoProfile(command='mpirun -np 16 pw.x',
                          pseudo_dir=here)
dft_calculator = Espresso(profile=profile,
                          pseudopotentials={
                              'Cu': 'Cu.paw.z_11.ld1.psl.v1.0.0-low.upf'
                              },
                          kpts=(11, 11, 11),
                          input_data=qeparam)

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
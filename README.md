# PyIdealLatt

![PyIdealLatt Logo](docs/assets/logo.bmp)

## Overview

PyIdealLatt is a lightweight Python package for generating ideal crystal structures, specifically designed to support large-scale model training in condensed matter physics and materials science. With minimal dependencies and a clean API, it provides fast and reliable crystal structure generation for AI/ML applications.

## Key Features

- 🚀 **Lightweight Design**: Minimal dependencies for easy deployment and integration
- 🔬 **Comprehensive Structure Types**: Support for cubic, FCC, BCC, diamond, and binary compound structures
- ⚡ **Smart Parameter Inference**: Automatic lattice parameter estimation from covalent radii
- 📊 **Experimentally Grounded**: Uses experimental lattice volume data for accuracy
- 🛠️ **Optimization Ready**: Built-in structure optimization capabilities via ASE integration
- 🧪 **ML-Friendly**: Optimized for generating training data in materials science AI applications

## Installation

```bash
pip install pyideallatt
```

## Quick Start

### Unary Structures

```python
from pyideallatt.bulk.unaries import build, BulkUnaryType

# Generate FCC copper structure
copper = build(BulkUnaryType.FACECENTEREDCUBIC, elem='Cu', a=3.6)
print(copper)

# Auto-infer lattice parameters
copper = build(BulkUnaryType.FACECENTEREDCUBIC, elem='Cu')
print(copper)

# Generate BCC iron structure
iron = build(BulkUnaryType.BODYCENTEREDCUBIC, elem='Fe')
print(iron)

# Generate diamond structure
diamond = build(BulkUnaryType.DIAMOND, elem='C')
print(diamond)
```

### Binary Compound Structures

```python
from pyideallatt.bulk.binaries import build, BulkBinaryType

# Generate X₂O structure (e.g., Ca₂O)
calcium_oxide = build(BulkBinaryType.X2Y, x='Ca')
print(calcium_oxide)

# Generate XO structure (e.g., MgO)
magnesium_oxide = build(BulkBinaryType.XY, x='Mg')
print(magnesium_oxide)

# Generate XO₂ structure (e.g., SiO₂)
silicon_dioxide = build(BulkBinaryType.XY2, x='Si')
print(silicon_dioxide)

# Support for different anions
titanium_nitride = build(BulkBinaryType.XY, x='Ti', y='N')
print(titanium_nitride)
```

### Structure Optimization

```python
from ase.calculators.emt import EMT

# Create structure
structure = build(BulkUnaryType.FACECENTEREDCUBIC, elem='Cu')

# Lattice and atomic position optimization
calculator = EMT()
structure.relax(calculator, fmax=0.01)

# Optimize only lattice edges
structure.relax_cell_edge(calculator, fmax=0.01)

# Uniform scaling optimization
structure.relax_cell_scale(calculator)
```

## Supported Crystal Structures

### Unary Structures

| Structure Type | Class Name | Description |
|----------------|------------|-------------|
| Simple Cubic | `SimpleCubic` | Simple Cubic (SC) |
| Body-Centered Cubic | `BodyCenteredCubic` | Body-Centered Cubic (BCC) |
| Face-Centered Cubic | `FaceCenteredCubic` | Face-Centered Cubic (FCC) |
| Diamond | `Diamond` | Diamond Structure |
| β-Tin | `BetaTin` | β-Tin Structure |

### Binary Compound Structures

| Formula | Class Name | Description |
|---------|------------|-------------|
| XO | `IdealXO` | 1:1 compounds (e.g., NaCl-type) |
| X₂O | `IdealX2O` | 2:1 compounds |
| XO₂ | `IdealXO2` | 1:2 compounds (e.g., TiO₂-type) |
| X₂O₃ | `IdealX2O3` | 2:3 compounds (e.g., Al₂O₃-type) |
| X₂O₅ | `IdealX2O5` | 2:5 compounds |

## Minimal Dependencies

PyIdealLatt is designed with a lightweight philosophy:

- **Core Dependencies**: `numpy`, `ase`, `scipy`, `pyyaml`
- **Optional Dependencies**: None required for basic functionality
- **Installation Size**: < 50MB total footprint
- **Startup Time**: < 1 second for typical operations

## Use Cases

### Machine Learning Training Data
```python
# Generate training dataset
from ase.io import write
from ase.calculators.cp2k import CP2K

# the calculator you want to use to make labels
calculator = CP2K(
    # your CP2K input parameters
)

elements = ['Cu', 'Ag', 'Au', 'Al', 'Fe', 'Ni']
structures = []

for elem in elements:
    fcc = build(BulkUnaryType.FACECENTEREDCUBIC, elem=elem)
    bcc = build(BulkUnaryType.BODYCENTEREDCUBIC, elem=elem)
    structures.extend([fcc, bcc])

# Export for ML training
for i, struct in enumerate(structures):
    s = struct.toase()
    s.calc = calculator
    s.get_potential_energy()
    write(f'training_data/struct_{i}.extxyz', struct.toase())
```

### High-Throughput Screening
```python
# Rapid structure generation for screening
compounds = [('Ti', 'O'), ('Al', 'O'), ('Si', 'O')]

for xsite, ysite in compounds:
    for bintyp in [BulkBinaryType.XY, BulkBinaryType.X2O, BulkBinaryType.XY2]::
        try:
            structure = build(bintyp, x=xsite, y=ysite)
            print(f"Generated: {structure.name_}")
        except:
            continue
```

## Contributing

We welcome contributions! Please follow these steps:

1. Fork this repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

- Project Homepage: https://github.com/kirk0830/PyIdealLatt
- Issues: https://github.com/kirk0830/PyIdealLatt/issues

---

**Note**: This package is primarily intended for academic research and educational purposes. Please verify generated structure parameters before use in production environments.
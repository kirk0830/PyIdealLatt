from pathlib import Path
from enum import Enum
from typing import TypeAlias

import numpy as np
from ase.data import covalent_radii, atomic_numbers

from pyideallatt.bulk.bulk import Bulk

class Spinel(Bulk):
    '''
    Spinel
    '''
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        raise NotImplementedError
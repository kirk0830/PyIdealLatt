from pyideallatt.baselatt import BaseLatt

class Nanowire(BaseLatt):
    '''
    Nanowire structure generator
    '''
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name_ = 'Nanowire'
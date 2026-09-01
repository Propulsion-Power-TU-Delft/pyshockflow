from pyshockflow import Driver
from pyshockflow import Config

inputFiles = [
    # 'input_ideal_standard.ini',
    # 'input_ideal_vinokur.ini',
    # 'input_real_arabi.ini',
    # 'input_real_vinokur.ini',
    # 'input_real_hllc.ini',
    'input_real_ausm+up.ini',
]

for inputFile in inputFiles:
    config = Config(inputFile)
    tube = Driver(config)
    tube.solve()

from pyshockflow import Driver
from pyshockflow import Config

inputFiles = ['input_Test%i.ini' % i for i in [1, 3, 5]]
for inputFile in inputFiles:
    config = Config(inputFile)
    tube = Driver(config)
    tube.solve()

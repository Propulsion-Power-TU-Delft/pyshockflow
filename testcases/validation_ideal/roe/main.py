from pyshockflow import Driver
from pyshockflow import Config

testNumbers = [1,3,5]
inputFiles = ['input_Test%i.ini' % i for i in testNumbers]
for inputFile in inputFiles:
    config = Config(inputFile)
    tube = Driver(config)
    tube.solve()

from pyshockflow import Driver
from pyshockflow import Config

configList = ['input.ini']

for configFile in configList:
    config = Config(configFile)
    tube = Driver(config)
    tube.solve()

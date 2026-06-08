from pyshockflow import Driver
from pyshockflow import Config

pressureList = [45]
configList = ['input_%ikPa.ini' %(pressure) for pressure in pressureList]

for configFile in configList:
    config = Config(configFile)
    tube = Driver(config)
    tube.solve()

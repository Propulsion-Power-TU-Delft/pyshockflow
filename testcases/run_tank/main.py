from pyshockflow import Driver
from pyshockflow import Config


config = Config('input_0.450.ini')
tube = Driver(config)
tube.solve()


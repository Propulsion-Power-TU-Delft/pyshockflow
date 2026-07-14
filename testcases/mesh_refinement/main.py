from pyshockflow import Driver
from pyshockflow import Config

ins = [
    'input_normal.ini', 
    'input_refined.ini',
    'input_refined_adaptive.ini',
    ]

for input in ins:
    config = Config(input)
    tube = Driver(config)
    tube.solve()

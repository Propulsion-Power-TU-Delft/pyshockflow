from pyshockflow import Driver
from pyshockflow import Config

inputs = [
    'input_100.ini', 
    'input_500.ini', 
    'input_1000.ini',
    ]

for input_file in inputs:
    print(f"Running simulation for {input_file}...")
    config = Config(input_file)
    tube = Driver(config)
    tube.solve()


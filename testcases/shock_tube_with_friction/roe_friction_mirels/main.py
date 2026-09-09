import sys
from pathlib import Path
from pyshockflow import Driver, Config

config = Config('input_Test1.ini')
tube = Driver(config)
tube.solve()
print("Mirels turbulent shock tube run finished successfully!")

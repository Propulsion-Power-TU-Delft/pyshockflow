from pyshockflow import Driver, Config

config = Config('input_Test1.ini')
tube = Driver(config)
tube.solve()
print("Reflected shock simulation with Mirels friction completed successfully!")

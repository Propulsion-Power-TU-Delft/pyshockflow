# pyshockflow #

### What is this repository for? ###

* Resolution of quasi-1D flow problems for ideal and real gases.
* Learning of numerical flux schemes, effects and implementation.
* Implementation and testing of new numerical schemes to solve the 1D Euler Equations.
* Study of non-ideal and non-classical gas dynamics effects.


### How do I get set up? ###

* Git clone the present folder in your system

* Go to the root folder and create a new python environment (pyshockflow) through yml file:
```bash
conda env create -f environment.yml
```

* Activate the environment:
```bash
conda activate pyshockflow
```

* Install the packages:
```bash
python -m pip install -e .
```

* Navigate to the regression tests folder and run the automated tests:
```bash
cd test/regression_tests/
pytest test_cases.py
```
this should print completed tests if everyting works correctly

* The testcases folder has several other cases, where input parameters are specified inside every `input.ini` file. To run a case:
```bash
python main.py
```



### Notes ###
* The code should be compatible with Mac OS, Windows and Linux systems.

* The technical details can be found in the related article ... (soon to be submitted)


### Contribution guidelines ###

* Validate the modifications by means of detailed test cases
* Ensure regression tests still work and extend them if needed

### Authors and contacts ###

- **Francesco Neri**, TU Delft, `f.neri@tudelft.nl`
- **Matteo Pini**, TU Delft, `m.pini@tudelft.nl`

### References ###

[1] Toro, Eleuterio F. Riemann solvers and numerical methods for fluid dynamics: a practical introduction. Springer Science & Business Media, 2013.

[2] Arabi, Sina, Jean-Yves Trépanier, and Ricardo Camarero. "A simple extension of Roe's scheme for real gases." Journal of Computational Physics 329 (2017): 16-28.

[3] Blazek, Jiri. Computational fluid dynamics: principles and applications. Butterworth-Heinemann, 2015.

[4] D'Alessandro, Simone, Marco Pizzarelli, and Francesco Nasuti. "A hybrid real/ideal gas mixture computational framework to capture wave propagation in liquid rocket combustion chamber conditions." Aerospace 8.9 (2021): 250.

[5] Hirsch, Charles. Numerical computation of internal and external flows: The fundamentals of computational fluid dynamics. Elsevier, 2007.

[6] Tom-Robin Teschner. [CFD University](https://cfd.university/).
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi as PI
from PyShockflow.styles import *

xStart = 0
xEnd = 1.0
D_tube = 30.3*1E-03
D_throat = 22.0*1E-03
A_tube = PI*(D_tube**2)/4
A_throat = PI*(D_throat**2)/4
print(f"The area of the tube is {A_tube:.3e}.")
print(f"The area of the throat is {A_throat:.3e}.")
print(f"The ratio between throat and tube section is {A_throat/A_tube:.6f}.")
nPoints = 100

x = np.linspace(xStart, xEnd, nPoints)

# parabolic function with throat in the middle
Acoeff = -4*(A_throat-A_tube)
Bcoeff = -Acoeff
Ccoeff = A_tube
print(f"The equation of the nozzle is: A = {Acoeff:.3e}*z^2 + {Bcoeff:.3e}*z + {Ccoeff:.3e}")
z = (x-x[0])/(x[-1]-x[0])
Area = (Acoeff*z**2 + Bcoeff*z + Ccoeff)

plt.figure()
plt.plot([xStart, (xStart+xEnd)/2, xEnd], [A_tube, A_throat, A_tube], 'o', label='Control Points')
plt.plot(x, Area, label='Quadratic nozzle')
plt.xlabel('x')
plt.ylabel('A')
plt.legend()
plt.grid(alpha=.3)

with open('nozzle.csv', 'w') as file:
    file.write('x,A\n')
    for i in range(len(x)):
        file.write('%.6e,%.6e\n' %(x[i], Area[i]))

plt.show()




import matplotlib.pyplot as plt
import numpy as np
import cmath as cm

# dy/dx = y
def f(x, y):
	return y


xmin = 0.0 	# int(input("Enter the start time value: "))
xmax = 2.0 	# int(input("Enter the final time value: "))
step = 0.1 	# int(input("Enter the step size: "))
y0 = 1.0 	# int(input("Enter the y starting value: "))

nsteps = int((xmax - xmin)/step)
print "Number of itterations: ", nsteps

x = np.linspace(xmin, xmax, nsteps)
y = np.linspace(xmin, xmax, nsteps)
y[0] = y0

for n in range(0, nsteps - 1):
	k1 = step*f(x[n], y[n])
	k2 = step * f(x[n] + step/2, y[n] + k1/2)
	y[n+1] = y[n] + k2

plt.figure()
plt.plot(x, y)
plt.show()
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import math

def dirichletInitial(x):
	return 0
	
def dirichletInitialDt(x):
	return 1
	
def kinkInitial(x):
	l = 1
	v = 2
	x0 = 0.5
	m = math.sqrt(l) / (v * math.sqrt(2))
	return v * math.tanh(m*(x - x0))
	
def initialdt(x):
	return 0
	
def zeroPotential(phi):
	return 0
	
def massPotential(phi):
	m = 3
	return m**2 * phi
	
def potential(phi):
	l = 1
	v = 2
	return l * (phi**2 - v**2) * phi
	
	
	
def dirichletBoundary(u, du):
	u[0] = 0
	u[1] = 0
	u[-1] = 0
	u[-2] = 0
	
def periodicBoundary(u, du):
	average = (u[0] + u[-1])/2
	u[0] = average
	u[-1] = average
	u[1] = average
	u[-2] = average
	average = (du[0] + du[-1])/2
	du[0] = average
	du[-1] = average
	
def kinkBoundary(u,du):
	u[0] = -2
	u[1] = -2
	u[-1] = 2
	u[-2] = 2
	

pictureEvery = 8


def getSecondDerivative(f,h):
	d2 = []
	
	#left side
	top = f[2] - 2*f[1] + f[0]
	d2.append(top / (h * h))
	
	top = f[3] - 2*f[2] + f[1]
	d2.append(top / (h * h))
	
	#middle
	N = len(f)
	for n in range(2,N-2):
		top = f[n+2] - 2*f[n] + f[n-2]
		d2.append(top / (4 * h * h))
		
		
	#right side
	top = f[-2] - 2*f[-3] + f[-4]
	d2.append(top / (h * h))
	
	top = f[-1] - 2*f[-2] + f[-3]
	d2.append(top / (h * h))

	return d2

def run(x0,x1,dx,finishTime,dt,initial,initialdt,potential,boundaryconditions):
	xs = np.arange(x0,x1,dx)
	u = [initial(x) for x in xs]
	du = [initialdt(x) for x in xs]
	
	fig, ax = plt.subplots(figsize=(5, 4))
	ax.set(xlim=(x0, x1), ylim=(-5, 5))
	line = ax.plot(xs, u, color='k', lw=2)[0]

	N = len(xs)
	def animate(i):
		d2 = getSecondDerivative(u,dx)
					
		#time step CHANGE TO RK2
		for n in range(0,N):
			c = d2[n] - potential(u[n])
			du[n] += c * dt
			u[n] += du[n] * dt
					
					
		#apply boundary conditions
		boundaryconditions(u,du)
			
		line.set_ydata(u)
		
		
	anim = FuncAnimation(fig, animate, interval= dt*1000, frames=int(finishTime/dt))
	anim.save('output.mp4')
	#plt.draw()
	#plt.show()
		
		
		
run(0,6,0.02, 15, 1/60, dirichletInitial, dirichletInitialDt, zeroPotential, dirichletBoundary)

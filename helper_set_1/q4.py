import matplotlib.pyplot as plt
import math
import os

os.makedirs("q4pics", exist_ok=True)

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
	u = []
	du = []
	xs = []
	
	x=x0
	while (x < x1):
		u.append(initial(x))
		du.append(initialdt(x))
		xs.append(x)
		x += dx

	N = len(xs)
	time = 0
	pictureCounter = 0
	while(time < finishTime):
		d2 = getSecondDerivative(u,dx)
			
		#time step CHANGE TO RK2
		for n in range(0,N):
			c = d2[n] - potential(u[n])
			du[n] += c * dt
			u[n] += du[n] * dt
			
			
		#apply boundary conditions
		boundaryconditions(u,du)
		
		
		if (pictureCounter % pictureEvery == 0):
			#draw picture
			plt.axis((x0,x1,-4.5,4.5))
			plt.plot(xs,u,'r')
			plt.xlabel("x")
			plt.ylabel("y")
			plt.grid(True)
			plt.savefig("q4pics/a" + format(int(pictureCounter/pictureEvery), '05d') + ".png")
			plt.clf()
		


		pictureCounter += 1
		time += dt
		
		
		
run(0,8,0.02, 15, 1/240, dirichletInitial, dirichletInitialDt, zeroPotential, dirichletBoundary)

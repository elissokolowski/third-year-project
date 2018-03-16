import math
import numpy as np

deltaphi = 2 * math.pi

def potential(phi, alpha):
	return (1 - math.cos(phi)) * (1 - alpha * math.sin(phi)**2)

def f1(x, alpha):
	return math.sqrt(2 * potential(x, alpha))
	
def integrate(f,x0,x1,dx):
	xs = np.arange(x0,x1,dx)
	
	answer = 0
	for x in xs:
		answer += f(x)*dx
	return answer
	
def integrate2d(f, x0, x1, y0, y1, dx, dy):
	xs = np.arange(x0,x1,dx)
	ys = np.arange(y0,y1,dy)
	
	answer = 0
	for x in xs:
		sumAtX = 0
		for y in ys:
			sumAtX += f(x, y)*dy
		answer += sumAtX * dx

	return answer

def f2(x,y, alpha):
	return (potential(x, alpha) + potential(y, alpha) - potential(x+y, alpha)) / math.sqrt(potential(x, alpha) * potential(y, alpha))

def getExactPhaseShift(alpha, gamma0, gamma1):
	M = integrate(lambda phi: f1(phi, alpha), 0, deltaphi, 0.01)
	deltaXgammaBeta = 1/(2 * M) * integrate2d(lambda x, y : f2(x,y, alpha), 0.0001, deltaphi, 0.0001, deltaphi, 0.01, 0.01)

	gammaBetas = np.arange(gamma0, gamma1, 0.2)
	deltaXs = [deltaXgammaBeta / gammaBeta for gammaBeta in gammaBetas]

	return gammaBetas, deltaXs, deltaXgammaBeta


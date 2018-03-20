import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import math
import sys
import os
from array import array

import phaseShift
import plotPhaseShifts



import twokinks as initial
#import kinkantikink as initial
alpha = 0
writeStep = 50                   # no. compute steps between each write
betas = [0.9, 0.96, 0.98, 0.99]  #set which values of beta to do






def potential(phi):
	s = math.sin
	c = math.cos
	return s(phi) * ( 1 - alpha * (s(phi)**2 + 2*c(phi) - 2*c(phi)**2))

x0 = - initial.startSpacing - 5                        # left simulation boundary
x1 = initial.startSpacing + 5                      # right simulation boundary
initial_phi = initial.doubleKinkInitial     # initial phi(x)
initial_pi = initial.doubleKinkInitialDot	# initial pi(x)

class List:
	def __init__(self, l):
		self.li = l
	def __iadd__(self, other):
		self.li = [x+y for x,y in zip(self.li, other.li)]
		return self
	def __add__(self, other):
		return List([x+y for x,y in zip(self.li, other.li)])
	def __mul__(self, other):
		return List([x * other for x in self.li])
	def __getitem__(self, key):
		return self.li[key]
	def __setitem__(self, key, value):
		self.li[key] = value
	def length(self):
		return len(self.li)
	def __len__(self):
		return len(self.li)

def getSecondDerivative(f,h):
	d2 = []

	#left side
	top = f[2] - 2*f[1] + f[0]
	d2.append(top / (h * h))

	#middle
	N = f.length()
	for n in range(1,N-1):
		top = f[n+1] - 2*f[n] + f[n-1]
		d2.append(top / (h * h))

	#right side
	top = f[-1] - 2*f[-2] + f[-3]
	d2.append(top / (h * h))

	return d2


def f(phi, pi, dx):
	d2 = getSecondDerivative(phi, dx)
	dVdp = [potential(p) for p in phi]
	return List([(d-v) for d,v in zip(d2, dVdp)])

def g(phi, pi, dx):
	return pi


def rk4(phi, pi, dt, dx):
	k1 = f(phi, pi, dx) * dt
	l1 = g(phi, pi, dx) * dt

	k2 = f(phi + l1 * 0.5, pi + k1 * 0.5, dx) * dt
	l2 = g(phi + l1 * 0.5, pi + k1 * 0.5, dx) * dt

	k3 = f(phi + l2 * 0.5, pi + k2 * 0.5, dx) * dt
	l3 = g(phi + l2 * 0.5, pi + k2 * 0.5, dx) * dt

	k4 = f(phi + l3, pi + k3, dx) * dt
	l4 = g(phi + l3, pi + k3, dx) * dt

	pi += (k1 + k2 * 2 + k3 * 2 + k4) * (float(1)/6)
	phi += (l1 + l2 * 2 + l3 * 2 + l4) * (float(1)/6)

	
def getdx(gamma, howManyPointsAcross):
	def inverse(y):
		return math.log(abs(math.tan(y/4))) / gamma
	width = abs(inverse(0.2) - inverse(2*math.pi - 0.2))
	return width / howManyPointsAcross


def run(beta, timeStepMethod, dx = 0):
	gamma = 1 / math.sqrt(1 - beta**2)
	print("\nRunning with beta = " + str(beta) + "    gamma = " + str(gamma))
	if dx == 0:
		dx = getdx(gamma, 80)
	dt = dx / 2
	print("dx = " + str(dx))
	xs = np.arange(x0,x1,dx)
	phi = List([initial_phi(x,beta) for x in xs])
	pi = List([initial_pi(x,beta) for x in xs])

	N = len(xs)
	
	folder = "images/" + str(beta)
	if not os.path.isdir(folder):
		os.mkdir(folder)
	if not os.path.isdir("data/"):
		os.mkdir("data/")
	dataFile = open("data/" + str(beta), 'wb')
	#writing header for binary data file
	#length, x0, x1, dx, timeBetween
	array('i' , [len(phi)]).tofile(dataFile)
	array('d', [x0,x1, dx, dt*writeStep]).tofile(dataFile)

	time = 0
	pictureCounter = 1
	finishTime = 2 * initial.startSpacing / beta #distance / speed
	print("Finish time is " + str(finishTime))

	while time < finishTime:
		timeStepMethod(phi, pi, dt, dx)
		initial.boundary(phi, pi)
		
		if pictureCounter % writeStep == 0:
			#array('d', phi.li).tofile(dataFile)
			
			plt.clf()
			plt.plot(xs,phi,'g')
			plt.plot(xs,[initial.doubleKink(x,time,beta) for x in xs],'b--')
			plt.xlabel("x")
			plt.ylabel(r"$\phi$")
			plt.axis([x0,x1, initial.plotHeight0, initial.plotHeight1])
			plt.grid(True)
			plt.savefig(folder + "/phi" + str(int(pictureCounter / writeStep)) + ".png")
		
		time += dt
		pictureCounter += 1
		
	dataFile.close()
		
	#get phase difference
	phaseS = phaseShift.getPhaseShift(xs, time, lambda x,t : initial.doubleKink(x,t,beta), phi, initial.pointToGetPhaseShiftAt)
	print("phase shift = " + str(phaseS))
	return phaseS

phaseShifts = [run(beta, rk4) for beta in betas]
betaGammas = [beta / math.sqrt(1 - beta**2) for beta in betas]

phaseShift = open("phaseshifts.txt", 'w')
for b, p in zip(betaGammas, phaseShifts):
	phaseShift.write(str(b) + " " + str(p) + "\n")
phaseShift.close()

plotPhaseShifts.plot(alpha, 1.5, 10)






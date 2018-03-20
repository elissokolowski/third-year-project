import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import math
import datetime as d
import time as tm
import fnmatch as fn
import sys
import os
from array import array

import phaseShift
import plotPhaseShifts



import twokinks as initial
#import kinkantikink as initial
alpha = 0.5
writeStep = 20 # no. compute steps between each write
betas = np.arange(0.96, 0.999, 0.003)  #set which values of beta to do

experiment = "test_2"
dataFolder = experiment + "/"
rawFolder = dataFolder + "raw/"
imagesFolder = dataFolder + "images/"
anim_file = dataFolder + "anim.mp4"


def gamma(v):
	return 1 / math.sqrt(1 - v**2)

def potential(phi):
	s = math.sin
	c = math.cos
	return s(phi) * ( 1 - alpha * (s(phi)**2 + 2*c(phi) - 2*c(phi)**2))

x0 = - initial.startSpacing - 5                        # left simulation boundary
x1 = initial.startSpacing + 5                      # right simulation boundary
initial_phi = initial.doubleKinkInitial     # initial phi(x)
initial_pi = initial.doubleKinkInitialDot	# initial pi(x)

def getSecondDerivative(f,h):
	d2 = []

	#left side
	top = f[2] - 2*f[1] + f[0]
	d2.append(top / (h * h))

	#middle
	N = len(f)
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
	return np.array([(d-v) for d,v in zip(d2, dVdp)])

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

def run(experiment, beta, timeStepMethod, dx = 0):
	gam = gamma(beta)
	dataFolder = experiment + "/"
	rawFolder = dataFolder + "raw/"
	imagesFolder = dataFolder + "images/"
	anim_file = dataFolder + "anim.mp4"

	print("\nRunning with beta = " + str(beta) + "    gamma = " + str(gam))

	dx = getdx(gam,80)
	dt = dx / float(1.2)

	print("dx = " + str(dx))

	xs = np.arange(x0,x1,dx)
	phi = np.array([initial_phi(x,beta) for x in xs])
	pi = np.array([initial_pi(x,beta) for x in xs])

	fig, ax = plt.subplots(figsize=(5, 4))
	ax.set(xlim=(x0, x1), ylim=(initial.plotHeight0, initial.plotHeight1))
	line = ax.plot(xs, phi, color='k', lw=1)[0]

	if not os.path.exists(dataFolder):
		os.mkdir(dataFolder)
	if not os.path.exists(imagesFolder):
		os.mkdir(imagesFolder)
	if not os.path.exists(rawFolder):
		os.mkdir(rawFolder)


	time = 0
	pictureCounter = 0
	finishTime = int(2 * initial.startSpacing / beta) #distance / speed
	steps = finishTime / dt
	print("Finish time is " + str(finishTime))
	print("Number of time steps: " + str(steps))

	while time < finishTime:
		timeStepMethod(phi, pi, dt, dx)
		initial.boundary(phi, pi)
		
		if pictureCounter % writeStep == 0:
			np.savetxt(rawFolder + str(pictureCounter) + ".txt", phi, delimiter="\n")
			sys.stdout.write("\r" + str(int(pictureCounter * 100 / steps)) + "%")
			sys.stdout.flush()

		time += dt
		pictureCounter += 1

	n_files = len(fn.filter(os.listdir(rawFolder), '*[0-9].txt'))
	def animate(i):
		fname = rawFolder + str(i * writeStep) + ".txt"
		with open(fname, "r") as f:
			content = f.readlines()

		content = [float(x.strip()) for x in content]

		line.set_ydata(content)

	anim = FuncAnimation(fig, animate, interval=dt * 1000, frames=n_files)
	anim.save(anim_file)
		
	#get phase difference
	phaseS = phaseShift.getPhaseShift(xs, time, lambda x,t : initial.doubleKink(x,t,beta), phi, initial.pointToGetPhaseShiftAt)
	print("phase shift = " + str(phaseS))
	return phaseS

phaseShifts = [run(experiment, beta, rk4) for beta in betas]
betaGammas = [beta / gamma(beta) for beta in betas]

phaseShift = open(dataFolder + "phaseshifts.txt", 'w')
for b, p in zip(betaGammas, phaseShifts):
	phaseShift.write(str(b) + " " + str(p) + "\n")
phaseShift.close()

plotPhaseShifts.plot(alpha, 1.5, 10, dataFolder)






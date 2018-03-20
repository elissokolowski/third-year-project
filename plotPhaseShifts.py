import matplotlib.pyplot as plt
import numpy as np

import exactPhaseShift


def plot(alpha, a, b):
	betas, phases = exactPhaseShift.getExactPhaseShift(alpha, a, b)
	
	f = open("phaseshifts.txt", 'r')
	betaGammas = []
	phaseShifts = []
	
	for line in f:
		nums = line.split(" ")
		betaGammas.append(float(nums[0].strip()))
		phaseShifts.append(float(nums[1].strip()))
	
	plt.clf()
	plt.scatter(betaGammas, phaseShifts)
	plt.plot(betas ,phases, 'g')
	plt.xlabel(r"$\beta \gamma$")
	plt.ylabel("phase shift")
	plt.grid(True)
	plt.savefig("phaseshift.png")

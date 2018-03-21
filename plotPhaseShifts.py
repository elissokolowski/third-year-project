import matplotlib.pyplot as plt
import numpy as np
import math
import os

import exactPhaseShift

def getdx(gamma, howManyPointsAcross):
	def inverse(y):
		return math.log(abs(math.tan(y/4))) / gamma
	width = abs(inverse(0.2) - inverse(2*math.pi - 0.2))
	return width / howManyPointsAcross


def gamma(v):
	return 1 / math.sqrt(1 - v**2)

class Settings:
	def __init__(self, experiment):
		settingsFile = open("experiments/" + experiment + "/settings.txt", 'r')
		self.kinkType = settingsFile.readline().rstrip()
		if self.kinkType != "kinkkink" and self.kinkType != "kinkantikink":
			print("Error: not valid kinks type")
		self.beta = float(settingsFile.readline())
		self.gam = gamma(self.beta)
		self.alpha = float(settingsFile.readline())
		dx = settingsFile.readline()
		if dx == "" or dx == "\n":
			self.dx = getdx(self.gam,80)
		else:
			self.dx = float(dx)


def plot(alpha, a, b):
	betas, phases = exactPhaseShift.getExactPhaseShift(alpha, a, b)
	
	possibleExperiments = os.listdir("experiments/")
	experiementsToPlot = []
	for possible in possibleExperiments:
		s = Settings(possible)
		if alpha == s.alpha:
			experiementsToPlot.append((possible, s))
	

	betaGammas = [s[1].beta * s[1].gam for s in experiementsToPlot]
	phaseShifts = []
	
	for experiment in experiementsToPlot:
		f = open("experiments/" + experiment[0] + "/phaseshift.txt")
		phaseShifts.append(float(f.readline().strip()))
	
	plt.clf()
	plt.scatter(betaGammas, phaseShifts)
	plt.plot(betas ,phases, 'g')
	plt.xlabel(r"$\beta \gamma$")
	plt.ylabel("phase shift")
	plt.grid(True)
	plt.savefig("phaseshift.png")

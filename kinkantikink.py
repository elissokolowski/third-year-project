import math

startSpacing = 10
pointToGetPhaseShiftAt = -math.pi
plotHeight0 = -7
plotHeight1 = 7

def kink1(x,t, beta):
	gamma = 1 / math.sqrt(1 - beta**2)
	return 4 * math.atan(math.exp(gamma*(x+startSpacing - beta * t)))
def kink1dot(x, beta):
	gamma = 1 / math.sqrt(1 - beta**2)
	return -4*beta*gamma*math.exp(gamma*(x+startSpacing))/ (math.exp(gamma*(x+startSpacing))**2 +1)

def kink2(x,t,beta):
	beta *= -1
	gamma = 1 / math.sqrt(1 - beta**2)
	return 4 * math.atan(math.exp(-gamma*(x-startSpacing - beta * t)))
def kink2dot(x,beta):
	beta *= -1
	gamma = 1 / math.sqrt(1 - beta**2)
	return 4*beta*gamma*math.exp(-gamma*(x-startSpacing))/ (math.exp(-gamma*(x-startSpacing))**2 +1)



def doubleKink(x,t,beta):
	return kink1(x,t,beta) + kink2(x,t,beta) - 2*math.pi
doubleKinkInitial = lambda x, beta : doubleKink(x,0,beta)
def doubleKinkInitialDot(x,beta):
	return kink1dot(x,beta) + kink2dot(x,beta)
def boundary(phi, pi):
	phi[0] = 0
	phi[-1] = 0
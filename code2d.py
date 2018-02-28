import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import math

lam = 10
v = 2
m = math.sqrt(lam) * v

def kink1(x):
	beta = -0.5
	gamma = 1 / math.sqrt(1 - beta**2)
	return v * math.tanh(m/math.sqrt(2) * gamma * x)
def kink1dot(x):
	beta = -0.5
	gamma = 1 / math.sqrt(1 - beta**2)
	return v * (1 - math.tanh(m/math.sqrt(2) * gamma * x)**2) * m * gamma * beta / math.sqrt(2)
	
def kink2(x):
	beta = 0.5
	gamma = 1 / math.sqrt(1 - beta**2)
	x0 = 12
	return v * math.tanh(-m/math.sqrt(2) * (gamma * x - x0))
def kink2dot(x):
	beta = 0.5
	gamma = 1 / math.sqrt(1 - beta**2)
	x0 = 12
	return -1 * v * (1 - math.tanh(-m/math.sqrt(2) * (gamma * x - x0))**2) * m * gamma * beta / math.sqrt(2)
	
def td_kink1(x,y):
	beta = -0.5
	gamma = 1 / math.sqrt(1 - beta**2)
	return v * math.tanh(m/math.sqrt(2) * gamma * x)

def td_kink1dot(x,y):
	beta = -0.5
	gamma = 1 / math.sqrt(1 - beta**2)
	return v * (1 - math.tanh(m/math.sqrt(2) * gamma * x)**2) * m * gamma * beta / math.sqrt(2)

def td_kinksin(x,y):
	return np.sin(x**2 + y**2) / (x**2 + y**2 + 0.0001)

def set_zero(x,y):
	return 0

def doubleKinkInitial(x):
	if (x < 0):
		return kink1(x)
	elif (x > 0):
		return kink2(x)
	else:
		return 2
def doubleKinkInitialDot(x):
	if (x < 0):
		return kink1dot(x)
	elif (x > 0):
		return kink2dot(x)
	else:
		return 0

def zeroPotential(phi):
	return 0

def doubleWellPotential(phi):
	return phi * lam * (phi**2 - v**2)

def td_doubleWellPotential(phi):
	return phi * lam * (phi**2 - v**2)

def td_boundary(phi,pi):
	x_pts, y_pts = phi.shape

	for j in range(0, y_pts):
		phi[0,j] = 0
		phi[x_pts - 1, j] = 0

	for i in range(0, x_pts):
		phi[i,0] = 0
		phi[i, y_pts - 1] = 0

def boundary(phi, pi):
	phi[0] = -2
	phi[-1] = 2
	

#1 = compare single point to exactResult
#2 = compare global error with exact
#3 = compare 3 different resolutions
errorMode = 2
#used for number 1 and 3
pointToCompareTo = 1
#used for 1 and 2
def exactResult(x,t):
	beta = -0.5
	gamma = 1 / math.sqrt(1 - beta**2)
	return v * math.tanh(m/math.sqrt(2) * gamma * (x - beta * t))
	
x0 = -3
y0 = -3
x1 = 3
y1 = 3
finishTime = 10
initial = td_kinksin
phiDotInitial = td_kink1dot
potential = zeroPotential
	
	
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

def calculateTotalError(phi, exact, xs, t):
	dx = xs[1] - xs[0]
	ErrorAtPoints = [(phi[n] - exact(xs[n], t))**2 for n in range(len(xs))]
	return math.sqrt(0.5 * dx * (2 * sum(ErrorAtPoints) - ErrorAtPoints[0] - ErrorAtPoints[-1]))

def getSecondDerivative(phi,h):
	
	d2 = np.zeros(phi.shape)
	x_pts , y_pts = phi.shape
	
	for j in range(0, y_pts):
		
		for i in range(0,x_pts):
			d2_x = 0
			d2_y = 0
			
			if i == 0:
				d2_x = phi[i+2,j] - 2*phi[i + 1,j] + phi[i,j]
			elif i == x_pts - 1:
				d2_x = phi[i-2,j] - 2*phi[i - 1,j] + phi[i,j]
			else:
				d2_x = phi[i+1,j] - 2*phi[i,j] + phi[i-1,j]
			
			if j == 0:
				d2_y = phi[i,j+2] - 2*phi[i,j+1] + phi[i,j]
			elif j == y_pts - 1:
				d2_y = phi[i,j-2] - 2*phi[i,j-1] + phi[i,j]
			else:
				d2_y = phi[i, j+1] - 2*phi[i,j] + phi[i,j-1]
			
			d2[i,j] = (d2_x + d2_y) / (h * h)
		
	return d2	
	
def f(phi, pi, dx):
	d2 = getSecondDerivative(phi, dx)
	pot = np.zeros(phi.shape)
	x_pts, y_pts = phi.shape

	for j in range(0, y_pts):
		for i in range(0, x_pts):
			pot[i,j] = potential(phi[i,j])

	return d2 - pot
		
def g(phi, pi, dx):
	return pi
		
		
def euler(phi, pi, dt, dx):
	pi += f(phi, pi, dx) * dt
	phi += g(phi, pi, dx) * dt

def rk2(phi, pi, dt, dx):
	k1 = f(phi, pi, dx) * dt
	l1 = g(phi, pi, dx) * dt
	
	k2 = f(phi + l1 * 0.5, pi + k1 * 0.5, dx) * dt
	l2 = g(phi + l1 * 0.5, pi + k1 * 0.5, dx) * dt
	
	pi += k2
	phi += l2
	
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
	

def run(dx, dt, timeStepMethod, outputDir = "output"):
	xs = np.arange(x0,x1,dx)
	ys = np.arange(y0,y1,dx)

	phi = np.zeros((len(xs),len(ys)))
	pi = np.zeros((len(xs),len(ys)))

	x_pts, y_pts = phi.shape
	for j in range(0, y_pts):
		for i in range(0, x_pts):
			phi[i,j] = initial(xs[i],ys[j])
			pi[i,j] = phiDotInitial(xs[i],ys[j])

	outputFile = outputDir + "/output.mp4"
	dataDir = outputDir + "/data"

	
	im = plt.imshow(phi, animated=True)

	def animate(i):
		timeStepMethod(phi, pi, dt, dx)
		
		#apply boundary conditions
		td_boundary(phi, pi)
		
		#add error to Error
		#redo errors in 2d
		#if errorMode == 1:
		#	Error[1].append(abs(phi[ErrorPointX, ErrorPointY] - exactResult(xs[ErrorPointX, ErrorPointY], dt * i)))
		#elif errorMode == 2:
		#	Error[1].append(calculateTotalError(phi, exactResult, xs, dt * i))
		#elif errorMode == 3:
		#	Error[1].append(phi)
			
		im.set_array(phi)
		return im,
	
		
	anim = FuncAnimation(fig, animate, interval= dt*1000, frames=int(finishTime/dt))
	anim.save(outputFile, fps=30)
	print("done")

fig = plt.figure()
run(float(1)/20, float(1)/30, rk4, "output125")
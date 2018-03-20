def getNforPointIncreasing(list, point):
	n = -1
	while list[n] > point:
		n-=1
	n+=1
	return n + (point - list[n]) / (list[n-1] - list[n])

def getXatPointIncreasing(xs, list, point):
	n = getNforPointIncreasing(list, point)
	between = n - int(n)
	return (1 - between) * xs[int(n)] + between * xs[int(n)-1]
	
def getNforPointDecreasing(list, point):
	n = 0
	while list[n] > point:
		n+=1
	n-=1
	return n + (point - list[n]) / (list[n+1] - list[n])

def getXatPointDecreasing(xs, list, point):
	n = getNforPointDecreasing(list, point)
	between = n - int(n)
	return (1 - between) * xs[int(n)] + between * xs[int(n)+1]
	

#assuming that both 'exact' and 'calculated' start above 'point' and decrease to below
def getPhaseShift(xs, t, exact, calculated, point):
	exactList = [exact(x, t) for x in xs]
	
	Xe = 0
	Xc = 0
	
	if exactList[0] < point:
		Xe = getXatPointIncreasing(xs, exactList, point)
		Xc = getXatPointIncreasing(xs, calculated, point)
	else:
		Xe = getXatPointDecreasing(xs, exactList, point)
		Xc = getXatPointDecreasing(xs, calculated, point)
	
	return abs(Xe - Xc)
	



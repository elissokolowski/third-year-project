def getNforPoint(list, point):
	n = 0
	while list[n] > point:
		n+=1
	n-=1
	return n + (point - list[n]) / (list[n+1] - list[n])

def getXatPoint(xs, list, point):
	n = getNforPoint(list, point)
	between = n - int(n)
	return (1 - between) * xs[int(n)] + between * xs[int(n)+1]

#assuming that both 'exact' and 'calculated' start above 'point' and decrease to below
def getPhaseShift(xs, t, exact, calculated, point):
	exactList = [exact(t, x) for x in xs]
	
	Xe = getXatPoint(xs, exactList, point)
	Xc = getXatPoint(xs, calculated, point)
	
	return abs(Xe - Xc)

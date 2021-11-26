import math
a = 6378
angle = 0.4*2*math.pi/360
beta = -(349.01-288.15)/0.55/10**6           #-3*10**-4
print(-1/4/6378, beta)
ke = 1/(1+a*beta)
print('ke', ke)

ans = ke*a*(math.cos(angle)/math.cos(angle+60/(ke*a)) -1)
print(ans)
ans = ((60*(1+beta*a)/math.cos(angle) + a*math.sin(angle))**2 - a**2*(math.sin(angle))**2)/(2*a*(1+beta*a))
print(ans)
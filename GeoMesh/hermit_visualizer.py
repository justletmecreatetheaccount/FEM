import numpy as np
import matplotlib.pyplot as plt
from numpy import pow, sqrt

# ==================================================================
h,x0,y0,r0,h0,d0,x1,y1,r1,h1,d1 = 0.100000,-0.500000,-1.000000,0.500000,0.005000,0.400000,0.250000,0.500000,0.125000,0.020000,0.100000
d = np.linspace(0, 1.052, 1000) #https://www.desmos.com/calculator/6wbg63f8dy
y = np.zeros(d.shape)
d_1 =  (1.052 - d1) # distance to the d1 border
d_0 = d0
# ==================================================================


def compute_hermit(d):
    y = h
    if (d <= d_0):
        a0 = h0
        a2 = (3*(h - h0))/pow(d0,2)
        a3 = (2*(h0 - h))/pow(d0,3)
        y = a0 + a2*pow(d,2) + a3*pow(d,3)
    if (d >= d_1):
        print(d)
        a0 = h1
        a2 = (3*(h - h1))/pow(d_1,2)
        a3 = (2*(h1 - h))/pow(d_1,3)
        y = min(h,a0 + a2*pow(d,2) + a3*pow(d,3))
        
    return y




# ==================================================================
for i,_d in enumerate(d):
    y[i] = compute_hermit(_d)
plt.plot(d, y, label='Hermite Interpolation')
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.savefig("Hermit.png")
# ==================================================================
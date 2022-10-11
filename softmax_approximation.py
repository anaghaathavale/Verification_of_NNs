import matplotlib.pyplot as plt
import numpy as np
import math
from sklearn.metrics import max_error
from scipy.optimize import minimize_scalar
import sympy as sp
from sympy import exp

arr=[]
def fun(t,u): #(t,u are the x co-ordinates of the end-points for the line)

	#x = np.array(np.linspace(-6.906, 0.0, num=100))
	x = np.array(np.linspace(t, u, num=100))

	y_array=[]
	for i in x:
		y_array.append((1/(1 + math.exp(-(i)))))
	y = np.array(y_array)

	#create basic scatterplot
	#plt.plot(x, y, 'o')

	#obtain m (slope) and b(intercept) of linear regression line
	m, b = np.polyfit(x, y, 1) #'1' defines the degree of the polynomial we want to fit.

	model = np.array([m,b])

	predict = np.poly1d(model)

	x = sp.Symbol('x', real=True)
	interval = [t, u]
	f = (1 / (1 + exp(-(x)))) - m * x - b
	fd1 = f.diff(x) #df/dx
	
	#print(fd1)
	minf = min(f.subs(x, interval[0]),f.subs(x, interval[1])) #substitute x with interval[0] and interval[1] 
	maxf = max(f.subs(x, interval[0]),f.subs(x, interval[1]))

	roots = [(xx, f.subs(x, xx)) for xx in sp.solve(fd1, x)]
	#print(roots)
	for i in roots:
	    if interval[0] <= i[0] <= interval[1]:
	        minf = min(minf, f.subs(x, i[0]))
	        maxf = max(maxf, f.subs(x, i[0]))
	
	#print(minf)
	#print(maxf)

	midpt = (t+u)/2.0

	if(maxf < 0.001):
		#print("\nfor {t} to {u}".format(t=t,u=u))
		print("elif (x <= {u}): \n    return {m} * x + {b}".format(u=u,m = m, b=b) )
		#print("y = {m}x + {b}\n".format(m = m, b=b))
		return 0
	
	else:
		fun(t,midpt)
		fun(midpt,u)

	#add linear regression line to scatterplot 
	#plt.plot(x, m*x+b)
	#plt.show()

fun(-6.906754778648554, 0.0)
fun(0.0, 6.906754778648554)

# Model of Disease state from Sagar
'''
x(t): Hidden disease state
    x' = a*A(t)(1-x) - (b + c*T(t))*x

g(t): Observable GI symptom burden
    g' = k*x*(1-g) - \pi*g 

d(t): Long term chronic degredation
    d' = \eps*[m*(g - g_c)_+(1-d) - \eta*T(t)*d]

A(t): Adverse load from observed data (fit from data)
    A(t) = w_1 * BP_d(t) + w_2 * HR_d(t) + w_3 * BM_d(t) + 
         + w_4 * W_d(t) + w_5 * C
'''

import numpy as np
import scipy as sp


def model(t, y, a, b, c, k, pi, eps, m, gc, n, A, T): 
    x = y[0]
    g = y[1]
    d = y[2]

    dx = (a*A*(1-x) - (b + c*T)*x)
    dg = (k*x*(1-g) - pi*g)
    dd = (eps*(m*max(0, g-gc)*(1-d) - n*T*d))

    z = np.zeros(3,)
    z[0] = dx
    z[1] = dg
    z[2] = dd

    return z

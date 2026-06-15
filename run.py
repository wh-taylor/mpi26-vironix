import numpy as np
from  scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import DSM as *


if __name__ == "__main__":

    a = 1
    b = 1
    c = 1
    pi = 1

    k = 1
    r = 1
    C = 1

    eps = 1
    m = 1 
    n = 1
    gc = 1

    BP = 2
    R = 2
    BM = 2
    W = 2

    t0 = 0
    tf = 10
    t_span = (t0, tf)

    y0 = np.array([0.5, 0.5, 0.5])
    A = 5
    T = 2
    sol = solve_ivp(ode, t_span, y0, args=(a, b, c, k, pi, eps, m, n, gc, BP, R, BM, W))



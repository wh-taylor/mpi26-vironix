# Model of Disease state from Sagar
'''
x(t): Hidden disease state
    x' = a*A(t)(1-x) - (b + c*T(t))*x

g(t): Observable GI symptom burden
    g' = k*x*(1-g) - \pi*g 

d(t): Long term chronic degredation
    d' = \eps*[m*(g - g_c)_+(1-d) - \eta*T(t)*d]

A(t): Adverse load from observed data
    A(t) = w_1 * BP_d(t) + w_2 * HR_d(t) + w_3 * BM_d(t) + 
         + w_4 * W_d(t) + w_5 * C

    Where: 



    


'''
function dydt = odefun(t, y)
    % STATES
    x = y(1);
    g = y(2);
    d = y(3);

    % CONSTANTS
    w1 = 1; w2 = 1; w3 = 1; w4 = 1; w5 = 1;
    a = 1; b = 1; c = 1;
    k = 0.5;
    r = 1;
    C=1;
    epsilon = 1;
    m = 1;
    n = 1;
    gc = 0.5;

    % INPUTS
    BP = 2;
    R  = 2;
    BM = 2;
    W  = 2;
    T = 1;

    % AUXILIARY
    A = w1*BP + w2*R + w3*BM + w4*W +w5*C;
    
    % ODES
    dxdt = a*A*(1 - x) - (b + c*T)*x;
    dgdt = k*x*(1 - g) - r*g;
    dddt = epsilon * ( m*max(g - gc,0)*(1 - d) - n*T*d );
    dydt = [dxdt; dgdt; dddt];
end
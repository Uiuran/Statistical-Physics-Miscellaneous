import numpy as np

def lingevin(x, t, func, seed= 1002384610730):

    n = len(x);
    
    np.random.seed(seed);
    noise = np.random.normal(size=n);
    y = np.zeros(n);
    
    y = func(x, t) + noise;

    return y;

def stochastic_int(func, x, t, param = ()):

    N = len(t);
    n = len(x);
    pn = len(param);     
    tau = t[1] - t[0];
    y = np.zeros((N,n));
    noise = np.random.normal(size=(N,n));
    couple = [];
    y[0,:] = x + tau*func(x, 0, param[0], param[1]) + np.sqrt(tau)*np.random.normal(size=n);

    for i in range(N-1):        
        B = np.random.normal(loc = 1.0, size = (6,6))*param[1];
        B = (B + B.T)/2.0;
        couple.append(B);
        y[i+1,:] = y[i,:] + tau*func(y[i,:], t[i], param[0], couple[i]) + np.sqrt(tau)*noise[i,:];
    
    return (y, couple);
    





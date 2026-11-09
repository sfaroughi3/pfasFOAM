import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root

x0 = 63396
x1 = -118250
x2 = 54854

ρb = 1500
αs = 2.49e-1
αaw = 2.49e-1
Fs = 0.5
Kf = 5.5e-5
Faw = 0.5
Kaw = 1
θ = 0.0199977
θs = 0.3
N = 0.85
R = 8.314
T = 293.15
a = 0.004
b = 0.107
σ0 = 0.071

ux = 0

Dx = 1.5e-2

Sw = θ/θs
Aaw = Sw * Sw * x2 + Sw * x1 + x0

def c_residual_fvm (Cvalnp1, Cvals, Cs2, Caw2, dC, dt, dx):
    res = Cvals.copy()
    Kawnew = (σ0 * b) / (R * T * (a + Cvalnp1))
    Aaw = Sw * Sw * x2 + Sw * x1 + x0
    Kawold =(σ0 * b) / (R * T * (a + Cvals))

    Cs2np1 = dt / (1 + αs * dt) * αs* (1 - Fs) * Kf * Cvalnp1**N + 1/(1+αs*dt) * Cs2
    Caw2np1 = dt / (1 + αaw * dt) * αaw * (1 - Faw) * Aaw * Kawnew * Cvalnp1  + 1/(1+αaw*dt) *  Caw2

    ddtterm = θ * (Cvalnp1 - Cvals) / dt
    diffval = Dx * θ * (Cvalnp1[2:] - 2 * Cvalnp1[1:-1] + Cvalnp1[:-2])/ dx /dx
    advval = ux * (Cvalnp1[2:] - Cvalnp1[:-2]) / dx

    cs1term = ρb * Fs * Kf * (Cvalnp1 ** N - Cvals ** N) / dt
    cs2term = ρb * αs * ( (1 - Fs) * Kf * Cvalnp1 ** N - Cs2np1)

    caw1term = Faw * Aaw * (Kawnew *Cvalnp1 - Kawold *Cvals) / dt
    caw2term = αaw * ( (1 - Faw) * Aaw * Kawnew * Cvalnp1 - Caw2np1)

    res[1:-1] = (ddtterm[1:-1] + cs1term[1:-1] + cs2term[1:-1] + caw1term[1:-1] + caw2term[1:-1]
                - diffval + advval)
    res[0] = (ddtterm[0] + cs1term[0] + cs2term[0] + caw1term[0] + caw2term[0]
                + ( -Dx*θ* ( Cvalnp1[1] -  Cvalnp1[0]) / dx /dx + θ*Dx*dC /dx ))
    res[-1] = (ddtterm[-1] + cs1term[-1] + cs2term[-1] + caw1term[-1] + caw2term[-1]
                - ( Dx*θ* (Cvalnp1[-2] - 3 * Cvalnp1[-1]) / dx /dx ))

    return res

def getcsol (dt, dC, nx):
    tarr = np.arange(0,300+dt,dt)
    x = np.linspace(0, 1, nx)
    dx = x[1] -x[0]
    cc = x[:-1]+dx/2

    Csolmat = np.zeros((len(tarr)+1, len(cc)))
    Csolmat[0,0] = 0
    Cs2solmat = np.zeros_like(Csolmat)
    Caw2solmat = np.zeros_like(Csolmat)

    for i in range(len(tarr)):
        cold = Csolmat[i,:]
        csold = Cs2solmat[i,:]
        cawold = Caw2solmat[i,:]
        sol = root(c_residual_fvm, cold, args=(cold, csold, cawold, dC, dt, dx))
        Csolmat[i+1,:] = sol.x
        Kaw = (σ0 * b) / (R * T * (a + Csolmat[i+1,:]))
        Cs2solmat[i+1,:] = dt/(1 + αs * dt) * αs* (1 - Fs) * Kf * Csolmat[i+1,:]**N + 1/(1+αs*dt) * csold
        Caw2solmat[i+1,:] = dt/(1 + αaw * dt) * αaw*(1 - Faw)*Aaw*Kaw * Csolmat[i+1,:] + 1/(1+αaw*dt)*cawold

    return Csolmat, Cs2solmat, Caw2solmat

Csol01, Cs2sol01, Caw2sol01 = getcsol(0.5, -0.5, 51)
Csol05, Cs2sol05, Caw2sol05 = getcsol(0.5, -5.0, 51)
Csol20, Cs2sol20, Caw2sol20 = getcsol(0.5, -10.0, 51)
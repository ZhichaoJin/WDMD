from os import X_OK
import time

from numpy.linalg import eigvals
import functions_dmd as dfc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from _wind import windowDMD
from scipy.integrate import odeint
import decimal

dt = 60/(4000*256) # time step of signal, 4000 rev/min with 256 steps
snapshots = pd.read_csv(r"\\nas.tu-clausthal.de\win-home$\\zj19\Desktop\\MA\\snapshotbase_shroud.csv", delimiter = " ", header = None).values
#snapshots = pd.read_csv(r"\\nas.tu-clausthal.de\win-home$\zj19\Desktop\optdmd\atilde.csv", delimiter = " ", header = None).values
X, Y = dfc.timeShift(snapshots) # create X(k+1) = Atilde * X(k), where X(k+1) is X2 and X(k) is X
s = np.shape(snapshots)[0]
U,S,V=dfc.SVD(X)
r=95
U = U[:,:r]
S_diag = S[:r,:r]
V = V[:r,:]

print(S_diag.shape,U.shape,V.shape)
Atilde=dfc.calcAtilde(U, S_diag, V, Y)
np.savetxt("atilde.csv", Atilde)

tspan = np.linspace(0, 60/4000, 256)
t = tspan[1:]
#print(t)
x_low,y_low=dfc.timeShift(Atilde)
n, m = len(x_low[:, 0]), len(x_low[0, :])
print (n,m)
w = 1
evalswindowDMD = np.empty((n, m), dtype=complex)
wdmd = windowDMD(n, w, 1)
Xw = x_low[:, :w]
Yw = y_low[:, :w]
print(np.linalg.matrix_rank(Xw),Yw.shape)
wdmd.initialize(x_low[:, :w], y_low[:, :w])
start = time.time()
for k in range(w, m):
    wdmd.update(x_low[:, k], y_low[:, k])
    evalswindowDMD[:, k] = np.log(np.linalg.eigvals(wdmd.A))/dt
end = time.time()
print("Window DMD, w=10, weighting = 1, time = " + str(end-start) + " secs")
np.savetxt("Wdmd_evalswindowDMD.csv", evalswindowDMD)

evals, W = wdmd.computemodes()
print("Y",Y.shape,"V",V.shape,"S",S_diag.shape,"W",W.shape)
Phi = dfc.dmdModes(V,S_diag,W,Y)
Phi, evals = dfc.cutModes(Phi, evals)
re =  np.real(evals)
im =  np.imag(evals)
np.savetxt("Wdmd_eigenvalues_real.csv", re)
np.savetxt("Wdmd_eigenvalues_imaginary.csv", im)
np.savetxt("Wdmd_mode_real.csv", np.real(Phi))
np.savetxt("Wdmd_mode_imaginary.csv", np.imag(Phi))

Phi_abs = dfc.scaleModes(Phi,snapshots)
np.savetxt("Wdmd_mode_abs.csv", np.real(Phi_abs))
print(Phi.shape)
DMDfreqs, DMDpower = dfc.dmdSpectrum(evals, dt, Phi, snapshots, s)
rowsPhi_abs = len(Phi_abs)
colsPhi_abs = len(Phi_abs[0])
Phi_mean = np.zeros(colsPhi_abs)
for i in range (colsPhi_abs):
    Phi_mean[i] = sum(Phi_abs[:,i])/rowsPhi_abs  
np.savetxt("wDMDfreqs_real.csv", np.real(DMDfreqs))
np.savetxt("wDMDfreqs_imaginary.csv", np.imag(DMDfreqs))
np.savetxt("wDMDpower.csv", DMDpower)
np.savetxt("wDMDpower_Phi.csv", Phi_mean)

import numpy as np
import matplotlib.pyplot as plt

from scipy.fft import fft

N = 1000
dt = 6/N
t = np.linspace(0,dt*N,N)

A = np.zeros((N,10))
for j in range(1,11):
    A[:,j-1] = np.sin(2*np.pi*0.02*j**2*t)

# B = np.matmul(A.transpose(),A)
# print(B)
# plt.imshow(B)
# plt.show()

func = lambda t: np.cos(2*np.pi*(3*(t-2)))*np.exp(-0*np.pi*(t-0)**2)
ftsine = np.fft.rfft(func(t),n=N)
f = np.linspace(0.0, 1.0/(2*dt), N//2)
f = np.fft.rfftfreq(N,dt)
mean = dt*np.sum(np.exp(-2*np.pi*1j*(3*t))*func(t))
print(mean)
fig = plt.figure()
plt.plot(t,(func(t)))
fig = plt.figure()
plt.plot(t,(np.exp(-2*np.pi*1j*(3*t))*func(t)).real)
plt.plot(t,(np.exp(-2*np.pi*1j*(3*t))*func(t)).imag)

fig = plt.figure()
plt.plot(f,2/N* ftsine[:N].real)
plt.plot(f,2/N*ftsine[:N].imag)
plt.plot(f,2/N* np.abs(ftsine[:N]))
plt.show()

# Number of sample points

N = 600

# sample spacing

T = 1.0 / 800.0

x = np.linspace(0.0, N*T, N)

y = np.sin(50.0 * 2.0*np.pi*x) + 0.5*np.sin(80.0 * 2.0*np.pi*x)

yf = fft(y)

xf = np.linspace(0.0, 1.0/(2.0*T), N//2)

import matplotlib.pyplot as plt

plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]))

plt.grid()

plt.show()

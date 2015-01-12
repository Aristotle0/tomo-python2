import numpy as np
import matplotlib.pyplot as plt
s = np.load('x.npy')
nt, ns = s.shape
print(nt, ns)
nt = nt // 2
sx = s[:nt, :]
xmax = np.amax(sx)
plt.figure()
plt.subplot(211)
for i in range(ns):
    plt.plot(sx[:,i]/xmax*5 + i, 'r')

t = np.load('T_raw.npy')
plt.subplot(212)
tmax = np.amax(t)
for i in range(ns):
    plt.plot(t[:, i]/tmax*5 + i, 'b')


v = np.load('v.npy')
plt.figure()
plt.plot(np.abs(v), 'b')
plt.show()

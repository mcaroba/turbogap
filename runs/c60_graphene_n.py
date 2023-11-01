import numpy as np
import matplotlib.pyplot as plt

mlg = np.loadtxt("graphene_n_mbd.dat")
c60_mlg = np.loadtxt("c60_graphene_n.dat")

c60 = -535.84852370

print(mlg)
print(c60_mlg)

n = c60_mlg[:,0]

c60_mlg = c60_mlg[:,1]

c60_mlg = c60_mlg - mlg - c60

plt.plot(n,c60_mlg)
plt.show()

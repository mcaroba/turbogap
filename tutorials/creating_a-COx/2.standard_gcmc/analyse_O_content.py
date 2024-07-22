import numpy as np
import matplotlib.pyplot as plt

# Read the mc.log file to see the acceptance and get averages

n_O = np.loadtxt("mc.log", usecols=(6,))

acceptance = np.loadtxt("mc.log", usecols=(2,), dtype=str)

N = len(acceptance)

accepted = np.zeros( ( N, ), dtype=int )

current=0
for i, (n,a) in enumerate( zip(n_O, acceptance) ):
    if a == "T":
        current = i

    accepted[i] = current

for i,(n,a) in enumerate(zip( n_O[accepted], accepted )):
    print(i, n, a)

equib_O_content = (1/N) * np.sum( n_O[accepted] )

print(" Equib. O content = ", equib_O_content)

fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_xlabel("MC Step")
ax.set_ylabel("# Oxygens")
ax.plot( np.arange(N), n_O[accepted] )
plt.show()

# with open("mc.log", "r") as f:

#     # We want obtain averages from MC, so we can read whether the move was accepted, and then

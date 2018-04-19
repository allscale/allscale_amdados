import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt("scalability_size.txt")
plt.plot(data[:,0],data[:,1])
plt.grid(True)
plt.xlabel("problem size")
plt.ylabel("execution time")
plt.title("Scalability in absolute values")
plt.savefig("abs-size-scalability.png", bbox_inches="tight")
plt.show()

plt.plot(data[:,0]/data[-1,0],
         data[:,1]/data[-1,1])
plt.grid(True)
plt.xlabel("(problem size) / (maximum problem size)")
plt.ylabel("(execution time) / (maximum execution time)")
plt.title("Scalability in relative values")
plt.axis("equal")
plt.axis([0, 1, 0, 1])
plt.savefig("rel-size-scalability.png", bbox_inches="tight")
plt.show()


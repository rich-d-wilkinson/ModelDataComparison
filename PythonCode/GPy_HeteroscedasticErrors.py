import GPy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

xlocs = np.arange(-10,10,1)
X = np.repeat(xlocs,20)
Xvars = np.repeat(np.arange(-10,10,1)**2/10.,20)+1
y = X + np.random.normal(0, Xvars)
plt.scatter(X, y)
plt.show()


k = GPy.kern.WhiteHeteroscedastic(1,X.size, Xvars)
print(k)

m = GPy.models.GPRegression(X.reshape(-1,1), y.reshape(-1,1), k)
m.optimize()

mu = m.predict(np.arange(-12,12,0.1).reshape(-1,1))

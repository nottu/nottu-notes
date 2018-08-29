import numpy as np
import matplotlib.pyplot as plt

x = np.arange(-2, 2., 0.1)

plt.plot(x, np.log(1+np.exp(-x)), label="Regresión")
plt.plot(x, np.exp(-x), label="Boosting")
plt.plot(x, [max(0, 1-v) for v in x], label="SVM")
plt.legend(loc='upper left')
plt.ylabel('C(x, y)')
plt.xlabel('yf(x)')
plt.ylim(-0.1,2)
plt.show()
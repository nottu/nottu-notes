# ./ej2 convex1 hestness convex1.txt 1E-10 1E-18 1E-10

fletcher = np.loadtxt('out.txt')
polak = np.loadtxt('out.txt')
hestness = np.loadtxt('out.txt')

fig, ax = plt.subplots()
plt.title(r'$f_k$')
ax.plot(fletcher[:,:1] , label='Fletcher-Reeves')
ax.plot(polak[:,:1] , label='Polak-Ribiere')
ax.plot(hestness[:,:1] , label='Hestenes-Stiefel')
legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')
plt.show()

fig, ax = plt.subplots()
plt.title(r'$||d_k||$')
ax.plot(fletcher[:,1:2] , label='Fletcher-Reeves')
ax.plot(polak[:,1:2] , label='Polak-Ribiere')
ax.plot(hestness[:,1:2] , label='Hestenes-Stiefel')
legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')
plt.show()

fig, ax = plt.subplots()
plt.title(r'$||g_k||$')
ax.plot(fletcher[:,2:3] , label='Fletcher-Reeves')
ax.plot(polak[:,2:3] , label='Polak-Ribiere')
ax.plot(hestness[:,2:3] , label='Hestenes-Stiefel')
legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')
plt.show()


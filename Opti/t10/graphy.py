import numpy as np
import matplotlib.pyplot as plt
import struct

def readMatrix(input_file):
    in_file = open(input_file,"rb")
    nrow = struct.unpack('i', in_file.read(4))[0]
    ncol = struct.unpack('i', in_file.read(4))[0]
    vt   = np.fromfile(in_file, dtype='float64')
    in_file.close()
    return(vt.reshape((nrow,ncol)))

datos = readMatrix('datos10.bin')
plt.plot(data[:,0], data[:,1], 'bx', label='datos')

initial = '0 0 15 -2 1'
p = [float(x) for x in initial.split(' ' )]
plt.plot(x, p[0]*x + p[1] + p[2]*np.exp(p[3] * (p[4] - x)**2), 'g--', label='X0')


model = '2.52803 1.69816 12.0305 -0.75611 1.49662'
p = [float(x) for x in model.split(' ' )]
plt.plot(x, p[0]*x + p[1] + p[2]*np.exp(p[3] * (p[4] - x)**2), 'r-', label='X*')
plt.legend(loc='upper left')

plt.savefig('grafica', dpi=300, figsize=(1, 2))
plt.clf()
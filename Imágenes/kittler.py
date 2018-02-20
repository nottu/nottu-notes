
import skimage
from skimage import io
import os
import matplotlib.pyplot as plt
import numpy as np
from skimage.color import rgb2gray
from skimage import img_as_ubyte

from math import sqrt
from math import log

def new_hist_T(p, a, v):
    _a = (1/v[0] - 1/v[1])
    _b = -2 * (a[0]/v[0] - a[1]/v[1])
    _c = 0
    _c += pow(a[0], 2)/v[0] - pow(a[1], 2)/v[1]
    _c += 2*(log(sqrt(v[0])) - log(sqrt(v[1])))
    _c -= 2*(log(p[0]) - log(p[1]))
    
    pos = ( -_b + sqrt( pow(_b, 2) - 4*_a*_c))/(2*_a)
    return pos

def get_params(img, T, mask=None):
    mask = img > T
    mask_i = img < T
    p = (float(img[mask_i].shape[0])/(img.shape[0] * img.shape[1]), float(img[mask].shape[0])/(img.shape[0] * img.shape[1]))
    print('img[mask_i]', img[mask_i].size)
    a = (img[mask_i].mean(), img[mask].mean())
    v = (img[mask_i].var(), img[mask].var())
    if v[0] == 0 : v[0] = 0.01
    if v[1] == 0 : v[1] = 0.01
    
    return(p,a,v)

def buscaUmbral(img, T):
    new_t = 0
    while( True ):
        new_t = new_hist_T(*get_params(img, T))
        if int(T) == int(new_t): break
        else: T = int(new_t)
    return int(new_t)

def ubralizaImg(img, T):
    mask = img > T
    img_u = np.zeros_like(img)
    img_u[mask] = 255
    return img_u

# if __name__ == '__main__':
#     filename = './info/celulas/1.jpg'
#     cel = img_as_ubyte(rgb2gray(io.imread(filename)))

#     plt.imshow(cel, cmap='gray')
#     plt.show()

#     T = buscaUmbral(cel, 180)
#     plt.imshow(ubralizaImg(cel, T), cmap='gray')
#     plt.show()
#     print('Umbral: ', T)
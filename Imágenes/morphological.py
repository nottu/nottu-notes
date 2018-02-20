import numpy as np
from numba import jit

@jit(nopython = True)
def getDiamondSE(size):
    size = (size * 2) + 1
    elem = np.zeros((size, size), dtype=np.uint8)
    mid = int(size / 2)
    for i in range(size):
        elem[i, mid] = 255
        if i < size-mid:
            for j in range(i):
                elem[i, mid - (j + 1)] = 255
                elem[i, mid + (j +1 )] = 255
        else:
            for j in range(size - i-1):
                elem[i, mid - (j + 1)] = 255
                elem[i, mid + (j + 1)] = 255
    return elem

@jit(nopython = True, nogil=True, parallel=True)
def __evalFunc(func, img, selem, i, j):
    _val = img[i, j]
    _pos = int(selem.shape[0]/2)
    for x in range(selem.shape[0]):
        val_y = x - _pos
        if val_y + i < 0 or  val_y + i >= img.shape[0] : continue
        for y in range(selem.shape[1]):
            val_x = y - _pos
            if val_x + j < 0 or val_x + j >= img.shape[1]: continue
            if selem[x, y] :
                if func == 0 and img[val_y + i, val_x + j] > _val:
                    _val = img[val_y + i, val_x + j]
                if func == 1 and img[val_y + i, val_x + j] < _val:
                    _val = img[val_y + i, val_x + j]
    return _val

@jit(nopython = True, nogil=True, parallel=True)
def __doOperation(func, img, selem):
    _img = np.zeros_like(img)
    
    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            _img[i, j] = __evalFunc(func, img, selem, i,j)
    return _img

def dilation(img, selem):
    _img = __doOperation(0, img, selem)
    return _img

def erosion(img, selem):
    _img = __doOperation(1, img, selem)
    return _img

def opening(img, selem):
    img2 = erosion(img, selem)
    return dilation(img2, selem)

def closing(img, selem):
    img2 = dilation(img, selem)
    return erosion(img2, selem)

def white_topHat(img, selem):
    img2 = opening(img, selem)
    for i in range(img.shape[0]) :
        for j in range(img.shape[1]) :
            v = int(img[i,j]) - int(img2[i,j])
            if v < 0 : v = 0
            img2[i,j] = v
    return img2

def black_topHat(img, selem):
    img2 = closing(img, selem)
    for i in range(img.shape[0]) :
        for j in range(img.shape[1]) :
            v =  int(img2[i,j]) - int(img[i,j])
            if v < 0 : v = 0
            img2[i,j] = v
    return img2



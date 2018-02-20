import os
from os import listdir
from morphological import *
from skimage import io
from skimage import img_as_ubyte
from skimage.color import rgb2gray

def main():
  mypath = './DRIVE'
  from os.path import isfile, join
  onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
  outdir = './THoutput'
  if not os.path.exists(outdir):
    os.makedirs(outdir)
  for filename in onlyfiles:
    fname = os.path.join(mypath, filename)
    img = img_as_ubyte(rgb2gray(io.imread(fname)))
    img2 = black_topHat(img, getDiamondSE(8))
    fname = os.path.join(outdir, filename)
    io.imsave(fname, img2)
if __name__ == '__main__':
  main()
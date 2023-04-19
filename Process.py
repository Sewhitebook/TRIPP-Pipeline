import time
import sep
import glob
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import astroalign as aa
from astropy.io import fits
from astropy.coordinates import SkyCoord
from photutils import centroids, CircularAperture
from photutils.aperture import aperture_photometry
from astropy.wcs.utils import skycoord_to_pixel, pixel_to_skycoord
from astropy.coordinates import ICRS, FK5
from regions import CirclePixelRegion, PixCoord

def mask(img, srcs, T):  # masking function. Takes science image and list of sources.
    b = sep.Background(img)
    dat = img-b.back() #background subtraction for app_phot
    for i, s in enumerate(srcs):
        x_cen = s['x']
        y_cen = s['y']
        radius = ((s['xmax'] - s['xmin'])/2)
        center_2 = PixCoord(x_cen, y_cen)
        circle = CirclePixelRegion(center_2, radius)
        mask = circle.to_mask()  # this and the previous line make a "circle" array with ones
        app2 = CircularAperture([x_cen, y_cen], radius)
        app_phot2 = aperture_photometry(dat, app2)  # aperture sum
        blend_val = float(app_phot2['aperture_sum'][0] / (np.pi * radius ** 2))# averaging the sum value over the area of the source
        if T == True:
            temp_blends.append((x_cen, y_cen, blend_val))
        d = mask.bbox.extent
        #print(d)
        if d[0] < 0 or d[1] < 0 or d[2] < 0 or d[3] < 0:
            continue
        l1, l2, l3, l4 = int(d[0]-.5), int(d[1]-.5), int(d[2]-.5), int(d[3]-.5)

        try:
            fill = dat[l3:l4, l1:l2] * (np.ones(mask.data.shape) - np.array(mask.data)) # used to fill noise data back in (esentially making the mask a circle again)
            dat[l3:l4, l1:l2] = mask.data*blend_val + fill #this just outright sets the values for each source to the blended value. note that the index slices are in y, x format. blame python
        except ValueError:
            continue
        #print(i)
    return dat

#files = sorted(glob.glob("C:/Users/lucaa/Documents/PipelineLite/M31S25/*.fz")) #file path for input data --> make sys.argv
files = sorted(glob.glob("C:/Users/lucaa/Documents/PipelineLite/Small/*.fz")) #file path for input data --> make sys.argv
outstr = "C:/Users/lucaa/Documents/PipelineLite/Residuals/{}.fz"
print(len(files))

hdus = [fits.open(f) for f in files]
data = [h[1].data for h in hdus]
aligned = [aa.register(i, data[0])[0] for i in data[0:]]
template = np.median(aligned, axis = 0)

bkg0 = sep.Background(template)
extracted0 = sep.extract(template-bkg0.back(), bkg0.globalrms*3, minarea =10, segmentation_map=False)

masked = [mask(img, extracted0, False) for img in aligned]
blend_template = mask(template, extracted0, False) #masks template --> look into true false thing, there was a reason for it at some point.

residuals = [np.subtract(i, blend_template) for i in masked]


counter = 0
for pics in residuals:
    bkg = sep.Background(pics)
    extracted = sep.extract(pics - bkg.back(), bkg.globalrms * 3, minarea=10, segmentation_map=False)
    header = fits.Header([fits.Card("History", "Extracted by sep")])
    hdul = fits.HDUList([fits.PrimaryHDU(pics), fits.BinTableHDU(data=extracted, header=header, name="SEP", ver=None)])
    hdul.writeto(outstr.format(counter), overwrite=True)
    plt.scatter(extracted['x'], extracted['y'], facecolors = 'none', edgecolors = 'r')
    counter += 1

plt.imshow(template, norm = LogNorm(vmin = 1, vmax = 70), origin = 'lower')
plt.show()


#plt.imshow(template, cmap="viridis", norm = LogNorm(vmin = 67, vmax = 212))
#plt.show()


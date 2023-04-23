import glob
import sep
import astropy.units as u
import astropy.coordinates as coord
import matplotlib.pyplot as plt
import numpy as np
import astroalign as aa
from astropy.io import fits
from astroquery.ipac.irsa import Irsa
from astroquery.sdss import SDSS
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord
from astropy.coordinates import SkyCoord
from matplotlib.colors import LogNorm
from photutils.aperture import aperture_photometry
from photutils import centroids, CircularAperture


files = sorted(glob.glob("/Users/lucaangeleri/Documents/LCO/sec32/*.fz"))
hdus = [fits.open(f) for f in files] #opens fits files so we can access header data
data = [h[1].data for h in hdus] #array for image data
aligned = [aa.register(i, data[0])[0] for i in data[0:]]
template = np.mean(aligned, axis = 0)

bkg_phot = sep.Background(template)
extracted_phot = sep.extract(template - bkg_phot.back(),  bkg_phot.globalrms*3, minarea =20, segmentation_map=False) #find sources in image
w = WCS(hdus[0][1].header) #WCS matrix object

reference = []
mags = []
for src in extracted_phot: #indexes extracted sources by try number to find reference stars
    x = src['x']
    y = src['y']
    coord = pixel_to_skycoord(x, y, w).transform_to('icrs')
    search = SDSS.query_crossid(coord, fields = ['ra', 'dec', 'psfMag_g', 'psfMagErr_g'], radius = 15*u.arcsec, region = False)
    if search: #if SDSS query returned results, continue
        if search['psfMag_g'] < 16 and search['type'] == 'STAR':
            reference.append([search['ra'], search['dec'], x, y, src['xmin'], src['xmax'], search['psfMag_g'], search['psfMagErr_g']]) #do we need to keep ra and dec?
            print(search['ra', 'dec', 'psfMag_g'])

for pnt in reference:
        app = CircularAperture([pnt[2], pnt[3]], (pnt[5] - pnt[4]) / 2) # takes instrumental mag at point
        app_phot2 = aperture_photometry(template - bkg_phot.back(), app)
        mag_instrumental = -2.5 * np.log10(app_phot2['aperture_sum'][0])
        mags.append([mag_instrumental, pnt[6][0], pnt[7][0]])
        print(mag_instrumental, pnt[6][0], pnt[7][0], pnt[0][0], pnt[1][0])
        plt.scatter(pnt[2], pnt[3], facecolors='none', edgecolors='r')



Minst = [m[0] for m in mags]
Mpsf =  [m[1] for m in mags]

p = np.polyfit(Minst, Mpsf, deg= 1)
x = np.arange(-15, 15)
y = [p[0]*i + p[1] for i in x]
print(p[0], p[1] )


plt.imshow(data[0], cmap = 'viridis', norm = LogNorm(vmin = 137, vmax = 320), origin='lower')
plt.show()

plt.plot(x, y)
for pnt in mags:

    plt.errorbar(pnt[0], pnt[1], yerr= pnt[2], marker  = 'o', linestyle ='')
    #plt.scatter(pnt[0], pnt[1])
plt.show()

"""

for i, img in enumerate(data):

    reference = []
    for src in extracted_phot: #indexes extracted sources by try number to find reference stars
        x = src['x']
        y = src['y']
        coord = pixel_to_skycoord(x, y, w).transform_to('icrs')
        search = SDSS.query_crossid(coord, fields = ['ra', 'dec', 'psfMag_g', 'psfMagErr_g'], radius = 15*u.arcsec, region = False)
        if search: #if SDSS query returned results, continue
            if search['psfMag_g'] < 15 and search['type'] == 'STAR':
                reference.append([search['ra'], search['dec'], x, y, src['xmin'], src['xmax'], search['psfMag_g'], search['psfMagErr_g']]) #do we need to keep ra and dec?
                print(search['ra', 'dec', 'psfMag_g'], i)
    for pnt in reference:
        app = CircularAperture([pnt[2], pnt[3]], (pnt[5] - pnt[4]) / 2) # takes instrumental mag at point
        app_phot2 = aperture_photometry(img - bkg_phot.back(), app)
        mag_instrumental = -2.5 * np.log10(app_phot2['aperture_sum'][0])
        mags.append([mag_instrumental, pnt[6][0], pnt[7][0]])
        print(mag_instrumental, pnt[6][0], pnt[7][0], pnt[0][0], pnt[1][0])
        plt.scatter(pnt[2], pnt[3], facecolors='none', edgecolors='r')
"""




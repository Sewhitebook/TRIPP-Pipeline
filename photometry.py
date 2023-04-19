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


files = sorted(glob.glob("C:/Users/lucaa/Documents/PipelineLite/Small/*.fz")) #file path for input data
hdus = [fits.open(f) for f in files] #opens fits files so we can access header data
data = [h[1].data for h in hdus] #array for image data


#prepare image for reference search:
bkg_phot = sep.Background(data[0])
extracted_phot = sep.extract(data[0] - bkg_phot.back(),  bkg_phot.globalrms*3, minarea =50, segmentation_map=False) #find sources in image
w = WCS(hdus[0][1].header) #WCS matrix object

#We want to scan an image for reference sources, as of writing script, search area is limited to 3 arcmin so we will have to do multiple searches (most likely)
reference_target = 10 #number of reference stars we want
reference_count = 0
attempt = 0

reference = []

while reference_count != reference_target: #indexes extracted sources by try number to find reference stars
    x = extracted_phot[attempt]['x']
    y = extracted_phot[attempt]['y']
    coord = pixel_to_skycoord(x, y, w).transform_to('icrs')
    search = SDSS.query_crossid(coord, fields = ['ra', 'dec', 'psfMag_g', 'psfMagErr_g'], radius = 15*u.arcsec, region = False)
    if search: #if SDSS query returned results, continue
        if search['psfMag_g'] < 14 and search['type'] == 'STAR':
            reference_count += 1
            reference.append([search['ra'], search['dec'], x, y, extracted_phot[attempt]['xmin'], extracted_phot[attempt]['xmax'], search['psfMag_g'], search['psfMagErr_g']]) #do we need to keep ra and dec?
            print(search['ra', 'dec', 'psfMag_g'])
    attempt +=1

mags= []
for pnt in reference:
    app = CircularAperture([pnt[2], pnt[3]], (pnt[5]-pnt[4])/2) #takes instrumental mag at point
    app_phot2 = aperture_photometry(data[0] - bkg_phot.back(), app)
    mag_instrumental = -2.5*np.log10(app_phot2['aperture_sum'][0])
    mags.append([mag_instrumental, pnt[6][0], pnt[7][0]])
    print(mag_instrumental, pnt[6][0], pnt[7][0], pnt[0][0], pnt[1][0])
    #plt.scatter(pnt[2], pnt[3], facecolors = 'none', edgecolors='r')

#plt.imshow(data[0], cmap = 'viridis', norm = LogNorm(vmin = 1, vmax = 100), origin='lower')
#plt.show()

for pnt in mags:
    plt.errorbar(pnt[0], pnt[1], yerr= pnt[2]*100, marker  = 'o', linestyle ='')
    #plt.scatter(pnt[0], pnt[1])
plt.show()






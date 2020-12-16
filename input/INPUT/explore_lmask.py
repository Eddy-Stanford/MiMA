import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf
#from cartopy import config
#import cartopy.crs as ccrs

with netcdf.netcdf_file('lmask.T42.nc', 'r') as cdf_data:
    #import pdb; pdb.set_trace()
    lmask = cdf_data.variables['lmask'][:][0]
    lmask = np.asarray(lmask)
    lat = cdf_data.variables['lat'][:]

    lon = cdf_data.variables['lon'][:]

    #ax = plt.axes(projection=ccrs.PlateCarree())
    plt.contourf(lon, lat, lmask)#, transform=ccrs.PlateCarree())

    #ax.coastlines()

    #plt.imshow(lmask, cmap='gray')
    plt.savefig('lmask.png')


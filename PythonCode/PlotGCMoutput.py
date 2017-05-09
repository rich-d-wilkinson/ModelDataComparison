import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib
import os

####################
# Load a GCM

#from R: HadCM3grid =  c(-180,178.75,1.25,-89.375,89.375,1.25)

lats = np.arange(-89.375,89.375,1.25)
longs = np.arange(-180,178.75+0.1, 1.25)
longgrid, latgrid = np.meshgrid(longs, lats)



GCM_dir = os.getenv("HOME")+'/Documents/Code/ModelDataComparison/DMC/Scripts/Model_data/CO2_anom/'
file_name = 'tczyi.txt'
gcm1 = np.genfromtxt(GCM_dir+file_name)
gcm_mask = np.genfromtxt(GCM_dir+'mask.txt', dtype='int')



############## Orthographic projection ####################################
# set up orthographic map projection with
# perspective of satellite looking down at 50N, 100W.
# use low resolution coastlines.
map = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
# draw coastlines, country boundaries, fill continents.
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.fillcontinents(color='coral',lake_color='aqua')
# draw the edge of the map projection region (the projection limb)
map.drawmapboundary(fill_color='aqua')
# draw lat/lon grid lines every 30 degrees.
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(-90,90,30))

x_ortho,y_ortho = map(longgrid, latgrid)
map.scatter(x_ortho.flatten()[gcm_mask-1], y_ortho.flatten()[gcm_mask-1], 15, marker='+', color='k')

yplot = np.zeros(longgrid.size)
yplot[gcm_mask-1] = gcm1
CS = map.contour(x_ortho,y_ortho,yplot.reshape(lats.size,longs.size),15,linewidths=1.5)
plt.show()


############### Equidistant Cylindrical Projection ####################################
# The simplest projection, just displays the world in latitude/longitude coordinates.
m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
m.drawcoastlines()
#m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
m.drawparallels(np.arange(-90.,91.,30.))
m.drawmeridians(np.arange(-180.,181.,60.))
m.drawmapboundary(fill_color='aqua')
plt.xlabel('lon')
plt.ylabel('lat')
levels = np.arange(-10,10,1)
from cmocean import cm as cm
m.contourf(longgrid,latgrid,yplot.reshape(lats.size,longs.size),15,levels=levels,
    cm = cm.thermal, linewidths=1.5)
m.colorbar()
plt.show()

import GPy
import sys
import os
sys.path.append(os.getenv("HOME") + "/Documents/Code/Emulation/GPyDifferentMetrics/")
from HaversineDist import Exponentialhaversine
from haversine import haversine
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib

obs_dir = os.getenv("HOME")+'/Documents/Code/ModelDataComparison/DMC/Scripts/Observation_data/P3+_SST_anom/'
file = 'lambda_10.txt'
observations = np.genfromtxt(obs_dir+file, skip_header=1)

X_obs = observations[:,0:2]
y_obs = observations[:,2].reshape(-1,1)

X_tmp = np.copy(X_obs)
X_tmp[:,0]=X_obs[:,1]
X_tmp[:,1] = X_obs[:,0]

### WARNING - TRANSPOSE lat and longs


###################
k1 = Exponentialhaversine(2, lengthscale=2000)
# how do we specify a different noice variance at each locations

m = GPy.models.GPRegression(X_obs, y_obs, k1)
m.optimize_restarts(10)
print(m)


#########################
#
#  Predict at GCM grid locations.
#


lats = np.arange(-89.375,89.375,1.25)
longs = np.arange(-180,178.75+0.1, 1.25)
longgrid, latgrid = np.meshgrid(longs, lats)

GCM_dir = os.getenv("HOME")+'/Documents/Code/ModelDataComparison/DMC/Scripts/Model_data/CO2_anom/'
file_name = 'tczyi.txt'
gcm1 = np.genfromtxt(GCM_dir+file_name)
gcm_mask = np.genfromtxt(GCM_dir+'mask.txt', dtype='int')

Xpred = np.column_stack((longgrid.flatten()[gcm_mask-1], latgrid.flatten()[gcm_mask-1]))
ypred, Vpred = m.predict_noiseless(Xpred)

yplot = np.zeros(longgrid.size)
yplot[gcm_mask-1] = ypred

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
m.scatter(X_obs[:,0],X_obs[:,1] )
plt.show()

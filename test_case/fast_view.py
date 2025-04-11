

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import numpy.ma as ma
import gsw 
import matplotlib.cm as cm
import xarray as xr
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec
from matplotlib import image
from matplotlib.dates import date2num
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cmocean as cmo
from matplotlib.colors import Normalize 
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from datetime import datetime, timedelta
import matplotlib.colors as colors
import matplotlib.dates as mdates
import matplotlib.colors as colors
import netCDF4 as nc
import glob
from numpy import save
from numpy import load
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset, inset_axes
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import cartopy.feature as cf
from cartopy.io import shapereader
from matplotlib.patheffects import Stroke
import shapely.geometry as sgeom


traj_all = f'OUTPUT_OPENDRIFT_NAME.nc'  # Use names for file retrieval
traj = xr.open_dataset(traj_all, chunks={'time': 'auto'})
lon_t = traj.lon.compute()
lat_t = traj.lat.compute()
stat_f = traj.status.compute()

#age_days = traj.age_seconds.compute() / 86400
# Filter particles by status
lon_active = lon_t.where(stat_f == 0, drop=True)
lat_active = lat_t.where(stat_f == 0, drop=True)
lon_stranded = lon_t.where(stat_f == 1, drop=True)
lat_stranded = lat_t.where(stat_f == 1, drop=True)



# Latitude and longitude boundaries
lat_min, lat_max = -30, -21
lon_min, lon_max = -49, -39
# Apply filters
mask_active = (
	(lat_active >= lat_min) & (lat_active <= lat_max) &
	(lon_active >= lon_min) & (lon_active <= lon_max) #&
	#(age_days <= 90)
)
mask_stranded = (
	(lat_stranded >= lat_min) & (lat_stranded <= lat_max) &
	(lon_stranded >= lon_min) & (lon_stranded <= lon_max)
)
# Apply masks
lon_active_filtered = lon_active.where(mask_active, drop=True)
lat_active_filtered = lat_active.where(mask_active, drop=True)
#age_active_filtered = age_days.where(mask_active, drop=True)
lon_stranded_filtered = lon_stranded.where(mask_stranded, drop=True)
lat_stranded_filtered = lat_stranded.where(mask_stranded, drop=True)


##############_-- plot the scatter over the ses grid --######################################## 
def plot_data_grid(ax, data_file,z):
	data = xr.open_dataset(data_file)
	lon = data['lon']
	lat = data['lat']
	lon = ma.masked_where(lon ==0, lon)
	lat = ma.masked_where(lat ==0, lat)
	lat_m = np.ma.masked_invalid(lat)
	lon_m = np.ma.masked_invalid(lon)
	ax.plot(lon_m,lat_m,'k',alpha=0.1,zorder = z)
	ax.plot(lon_m.T,lat_m.T,'k',alpha=0.1,zorder = z)


def prop_plot(ax,lon,lat,prop,cmap,bar_title,zorder_value,norm = None,specs_x = False,specs_y = False,colorbar = False):
	cf1 = ax.pcolormesh(lon,lat,prop,norm = norm,cmap = cmap,shading='nearest',zorder = zorder_value)
	if norm is None:
		norm = cm.colors.Normalize(vmin=vmin,vmax=vmax)
		cf1 = ax.contourf(lon,lat,prop,np.arange(vmin,vmax,step),extend = 'max',cmap = cmap,norm = norm,zorder=zorder_value)

	ax.patch.set_edgecolor('black')
	ax.patch.set_linewidth('0.5') 

	if colorbar is True:
		cb1 = plt.colorbar(cf1,ax=ax,fraction=0.04, pad=0.04,extend='both')
		cb1.set_label(bar_title, size='small')

	if specs_x is True:
		gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, 
						linewidth=0.3, color='k', alpha=0.7, linestyle='--')
		gl.top_labels = False
		gl.left_labels = False
		gl.right_labels = False
		gl.bottom_labels = True
		gl.xlines = False
		gl.ylines = False

		gl.xformatter = LONGITUDE_FORMATTER
		gl.yformatter = LATITUDE_FORMATTER
		gl.xlabel_style = {'size': 5, 'color': 'k'}
		gl.ylabel_style = {'size': 5, 'color': 'k'}

	if specs_y is True:
		gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
						linewidth=0.3, color='k', alpha=0.7, linestyle='--')
		gl.top_labels = False
		gl.left_labels = True
		gl.right_labels = False
		gl.bottom_labels = False
		gl.xlines = False
		gl.ylines = False

		gl.xformatter = LONGITUDE_FORMATTER
		gl.yformatter = LATITUDE_FORMATTER
		gl.xlabel_style = {'size': 5, 'color': 'k'}
		gl.ylabel_style = {'size': 5, 'color': 'k'}
	
	return ax

	

ecomfile = 'ses_jan_2020.nc'
data = xr.open_dataset(ecomfile)
lon = np.asarray(data.lon)
lat = np.asarray(data.lat)
depth = data.depth

df_ses = pd.read_csv('OPENDRFIT_DIRECTORY_GRID_LAND_POINTS/model_grid_SSVBES_withLandPoints_999_work.csv',delimiter='\t+')

lon_ses = df_ses['lon'].to_numpy()
lat_ses = df_ses['lat'].to_numpy()
depgrid = df_ses['depgrid'].to_numpy()

lon_land = lon_ses[np.where(depgrid < 0)]
lat_land = lat_ses[np.where(depgrid < 0)]

####ploting
plt.close('all')
fig=plt.figure(figsize=(7,5))
gs = gridspec.GridSpec(nrows=1,ncols=1,width_ratios=[1], height_ratios=[1])
gs.update(left=0.09, right=0.98, wspace=0.2, hspace=0.2, top =1, bottom = 0)
vmin = 0
vmax = 25
cmap = cmo.cm.deep

levels = 30
bar_title = "Depth [m]"
ax = plt.subplot(gs[0,0],projection=ccrs.PlateCarree())
ax.set_ylim(bottom = -24.15,top = -23.85)
ax.set_xlim(left = -46.5,right = -46.05)
tt_0 =' '
ax.set_title(tt_0,loc='center',fontsize='medium')
prop_plot(ax,lon,lat,(depth),cmap,tt_0,bar_title,30,vmin,vmax,specs_x = True,specs_y = True,colorbar = False)
#land points and releasing points
plot_data_grid(ax,ecomfile,10)
ax.scatter(-46.30653,-23.99264, marker = '*',c='orange',s= 30, zorder= 30,label = "INF 1")
ax.scatter(-46.37754,-23.9736, marker = 's',c='blue',s= 30, zorder= 30,label = "INF 2")
ax.scatter(-46.30752,-23.94844, marker = 'd',c='magenta',s= 30, zorder= 30,label = "MID")
ax.scatter(-46.35142,-24.02514, marker = 'o',c='k',s= 30, zorder= 30,label = "Cubatão")

ax.text(-46.297521,-23.966686 + 0.01, "Santos Harbor", fontsize = 8, zorder= 30)
ax.text(-46.222965,-23.908304 + 0.01, "Bertioga Channel", fontsize = 8, zorder= 30)
ax.text(-46.419616,-23.958953 + 0.01, "São Vicente Channel", fontsize = 8, zorder= 30)
ax.text(-46.375414,-23.881597 + 0.01, "Piaçaguera Channel", fontsize = 8, zorder= 30)

#chanels
lon_c = [-46.345687,-46.339669,-46.329811,-46.323199,-46.317632,-46.312012]
lat_c = [-23.971564,-23.971990,-23.974045,-23.976469,-23.980284,-23.985568]
t_c = ["C1","C2","C3","C4","C5", "C6"]
ax.scatter(lon_c,lat_c,marker = '.',c='red',s= 30, zorder= 30,label = "Channels")
for i, txt in enumerate(t_c):
	ax.annotate(txt, (lon_c[i], lat_c[i]),zorder = 50,fontsize = 'small')

ax.legend(title ='Seed points', loc ='lower center')

##########---> Here we scatter the particles modelled from opendrift <-----######
ax.scatter(
	lon_active_filtered, lat_active_filtered, s=0.25, c='purple', alpha=0.6, transform=ccrs.PlateCarree(), label='Active Particles'
)

if lon_stranded_filtered.count() > 0:
	ax.scatter(
		lon_stranded_filtered, lat_stranded_filtered, s=0.25, color='gold',
		alpha=0.6, transform=ccrs.PlateCarree(), label='Stranded Particles'
	)


#colorbar
cbar_1 = inset_axes(ax, width="70%", height="3%", loc=2) 
norm_1 = mpl.colors.Normalize(vmin=vmin,vmax=vmax)
cb1 = mpl.colorbar.ColorbarBase(cbar_1, cmap=cmap,norm=norm_1,extend = 'max',orientation='horizontal')
cb1.set_label(bar_title, size='medium')
cbar_1.xaxis.set_ticks_position('bottom')

plt.savefig("partilces_and_map.png", dpi = 300)











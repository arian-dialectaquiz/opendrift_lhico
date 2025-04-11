"""
Simulação anual para os 7 pontos na Baía (6 canais, 2 desembocaduras) para análise de exportação.
"""

#####################################
open_dir = '/data1/Deriva/opendrift_lhico/opendrift_ecom/' #diretório da pasta opendrift monterrey
#open_dir = '/home/arian/Deriva/GIT/opendrift' #diretório da pasta opendrift tripoli
import sys
sys.path.append(open_dir)

import numpy as np
import xarray as xr
from opendrift.readers import reader_ECOM_S2Z_many_grids
from opendrift.models.oceandrift import OceanDrift
from opendrift.models.plastic_model_ECOM import PlastDrift

from datetime import datetime, timedelta

from matplotlib.colors import ListedColormap
import cartopy.feature as cfeature
cmap = ListedColormap((cfeature.COLORS['water'],
										   cfeature.COLORS['land']))
#################### --- MAIN FUNCTION TO RUN AND SAVE --- #########

def run_single_place(reader_pcse,dt,lon_seed,lat_seed,zpos,radius_seed,number_seed,steps,time_steps,time_step_output,outfile_name,var_out):

		## -- general configurations -- ##
		    
		o = PlastDrift(loglevel=0)
		o.add_reader([reader_pcse])
		o.set_config('general:use_auto_landmask', False)
		var_out = var_out
		o.set_config('drift:current_uncertainty_uniform',0.05)
		o.set_config('general:coastline_action', 'stranding')
		o.set_config('environment:fallback:x_sea_water_velocity', None)
		o.set_config('environment:fallback:y_sea_water_velocity', None)
		## -- wavesconfigurations --#
		o.get_config('drift:stokes_drift') #check where in the plastic_ecom_model this is activate
		for i in range(0, len(dt),1):
				o.seed_elements(lon=lon_seed,lat=lat_seed,z= zpos,radius=radius_seed,number=number_seed,
												time=reader_pcse.start_time + timedelta(hours=int(dt[i])))

		o.run(steps=steps,time_step=time_steps,time_step_output=time_step_output,outfile = 'TRAJ_' + outfile_name,export_variables = var_out)
		o.write_netcdf_density_map(filename='DENS_'+ outfile_name)



##############------ 2020 ----###############################################
ecomfile_stokes_path = '/data1/Deriva/opendrift_lhico/opendrift_ecom/ecom_stokes_19_03.nc'
nc = xr.open_dataset(ecomfile_stokes_path)
#reader importation
reader_pcse = reader_ECOM_S2Z_many_grids.Reader(ecomfile_stokes_path)

## -- time variables -- ##
ecom_t = np.asarray(nc['time']) #time array of netcdf file with 8784 steps
ecom_dt = ecom_t[1] - ecom_t[0]
sec_dt = ecom_dt.astype('timedelta64[s]')

## -- general configurations -- ##
hours = 8688 #362 for the certantiy of running
time_steps       = 60 # timestep of 1 min in seconds
time_step_output = 86400 # timestep of 24 h in seconds

hours_to_minutes = (3600/time_steps)
steps =  hours * hours_to_minutes
#(steps/time_step_output)*time_steps = N days

## --- locations --###
lon_seed = [-46.30653,-46.37754,-46.345687,-46.339669,-46.329811,-46.323199,-46.317632,-46.312012]
lat_seed = [-23.99264,-23.9736,-23.971564,-23.971990,-23.974045,-23.976469,-23.980284,-23.985568]
t_seed = ["INF1","INF2","C1","C2","C3","C4","C5", "C6"]
radius_seed,number_seed = 10,42
var_out = ['latitude','longitude']
####################### ------------------- z0 - -----------------########
########### ------ INF1 --------############
zpos = np.linspace(-5,0,number_seed)

outfile_name_inf1  = t_seed[0] + "_20_STOKES.nc"

# 366912 to be seeded
dt_inf1 = np.arange(0,hours,1)
run_single_place(reader_pcse,dt_inf1,lon_seed[0],lat_seed[0],zpos,radius_seed,number_seed,steps,time_steps,time_step_output,outfile_name_inf1,var_out)


#########------ INF2 ----#############
outfile_name_inf2  = t_seed[1] + "_20_STOKES.nc"

dt_inf2= np.arange(0,hours,1)
run_single_place(reader_pcse,dt_inf2,lon_seed[1],lat_seed[1],zpos,radius_seed,number_seed,steps,time_steps,time_step_output,outfile_name_inf2,var_out)



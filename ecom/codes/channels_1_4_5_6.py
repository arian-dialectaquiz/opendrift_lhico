"""
Simulação anual para os 6 canais para análise de exportação.
Canais considerados abertos por 48h após a identificação do dia com pluviosidade > média + 2sigma
lançamentos nesses momentos, porém durando o ano todo como ativo.
"""

#####################################
open_dir = '/phd/Outros_trabalhos/Deriva' #diretório da pasta opendrift
import sys
sys.path.append(open_dir)

import numpy as np
import xarray as xr
from opendrift.readers import reader_ECOM_S2Z_many_grids
from opendrift.models.plastic_model_ECOM import PlastDrift

from datetime import datetime, timedelta

from matplotlib.colors import ListedColormap
import cartopy.feature as cfeature
cmap = ListedColormap((cfeature.COLORS['water'],
					   cfeature.COLORS['land']))
dan = '/home/danilo/phd/Outros_trabalhos/Deriva/Plastic_case/data'
#################### --- MAIN FUNCTION TO RUN AND SAVE --- #########
## --- locations --###
lon_seed = [-46.30653,-46.37754,-46.345687,-46.339669,-46.329811,-46.323199,-46.317632,-46.312012]
lat_seed = [-23.99264,-23.9736,-23.971564,-23.971990,-23.974045,-23.976469,-23.980284,-23.985568]
t_seed = ["INF1","INF2","C1","C2","C3","C4","C5", "C6"]
radius_seed,number_seed = 10,42
var_out = ['latitude','longitude']

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
ecomfile_stokes_path = '/home/danilo/phd/Outros_trabalhos/Deriva/Plastic_case/T2/ecom_with_ww3/2020/*.nc'
ncin = xr.open_mfdataset(ecomfile_stokes_path,chunks={'time': 1}, concat_dim='time') #from 2020-01-01T00:52:30.000000000 to 2021-01-01T23:52:30.000000000


#reader importation
reader_pcse = reader_ECOM_S2Z_many_grids.Reader(ecomfile_stokes_path)

## -- time variables -- ##
ecom_t = np.asarray(ncin['time']) #time array of netcdf file with 8784 steps
ecom_dt = ecom_t[1] - ecom_t[0]
sec_dt = ecom_dt.astype('timedelta64[s]')

## -- general configurations -- ##
time_steps       = 60 # timestep of 1 min in seconds
time_step_output = 86400 # timestep of 24 h in seconds
hours_to_minutes = (3600/time_steps)

######## --- CHANELS ---####
#here depends the "dt" variables (seeding moments) on the amount of hours that the channels were activated

open_20 = np.load('/home/danilo/phd/Outros_trabalhos/Deriva/Plastic_case/T2/t20_open_doors_index.npy')
ecom_t_seed_channels = ecom_t[open_20]

zpos = np.linspace(-5,0,number_seed)

##C1 - ok
outfile_name_c1  = t_seed[2] + "_20_STOKES_new.nc"
dt_c1 = open_20
##hours in between open doors
hours_open = np.diff([ecom_t_seed_channels[0],ecom_t_seed_channels[-1]]).astype('timedelta64[h]')[0].astype('int')+ 1
steps_open =  hours_open * hours_to_minutes
run_single_place(reader_pcse,dt_c1,lon_seed[2],lat_seed[2],zpos,radius_seed,number_seed,steps_open,time_steps,time_step_output,outfile_name_c1,var_out)


##C2 - ok

#outfile_name_c2  = t_seed[3] + "_20_STOKES.nc"
#dt_c2 = open_20
#
##hours in between open doors
#hours_open = np.diff([ecom_t_seed_channels[0],ecom_t_seed_channels[-1]]).astype('timedelta64[h]')[0].astype('int') + 1
#steps_open =  hours_open * hours_to_minutes
#run_single_place(reader_pcse,dt_c2,lon_seed[3],lat_seed[3],zpos,radius_seed,number_seed,steps_open,time_steps,time_step_output,outfile_name_c2,var_out)

#020-11-19 23:51:30
#c3 - ok

#outfile_name_c3  = t_seed[4] + "_20_STOKES.nc"
#dt_c3 = open_20
#
##hours in between open doors
#hours_open = np.diff([ecom_t_seed_channels[0],ecom_t_seed_channels[-1]]).astype('timedelta64[h]')[0].astype('int')+ 1
#steps_open =  hours_open * hours_to_minutes
#run_single_place(reader_pcse,dt_c3,lon_seed[4],lat_seed[4],zpos,radius_seed,number_seed,steps_open,time_steps,time_step_output,outfile_name_c3,var_out)

#c4 - ok

outfile_name_c4  = t_seed[5] + "_20_STOKES.nc"
dt_c4 = open_20

##hours in between open doors
hours_open = np.diff([ecom_t_seed_channels[0],ecom_t_seed_channels[-1]]).astype('timedelta64[h]')[0].astype('int')+ 1
steps_open =  hours_open * hours_to_minutes
run_single_place(reader_pcse,dt_c4,lon_seed[5],lat_seed[5],zpos,radius_seed,number_seed,steps_open,time_steps,time_step_output,outfile_name_c4,var_out)

#c5 - ok

outfile_name_c5  = t_seed[6] + "_20_STOKES.nc"
dt_c5 = open_20

#hours in between open doors
hours_open = np.diff([ecom_t_seed_channels[0],ecom_t_seed_channels[-1]]).astype('timedelta64[h]')[0].astype('int')+ 1
steps_open =  hours_open * hours_to_minutes
run_single_place(reader_pcse,dt_c5,lon_seed[6],lat_seed[6],zpos,radius_seed,number_seed,steps_open,time_steps,time_step_output,outfile_name_c5,var_out)
#
##c6 ok 
#
outfile_name_c6  = t_seed[7] + "_20_STOKES.nc"
dt_c6 = open_20

#hours in between open doors
ours_open = np.diff([ecom_t_seed_channels[0],ecom_t_seed_channels[-1]]).astype('timedelta64[h]')[0].astype('int')+ 1
#steps_open =  hours_open * hours_to_minutes

run_single_place(reader_pcse,dt_c6,lon_seed[7],lat_seed[7],zpos,radius_seed,number_seed,steps_open,time_steps,time_step_output,outfile_name_c6,var_out)



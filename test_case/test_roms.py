	

#####################################
open_dir = '/data1/Deriva/opendrift_lhico/roms/' #diretório da pasta opendrift
import sys
sys.path.append(open_dir)
import glob as glob
import numpy as np
from opendrift.readers import reader_ROMS_native ,reader_ROMS_lhico
from opendrift.models.oceandrift import OceanDrift
from datetime import datetime, timedelta
import xarray as xr 
import os

directory = "/data1/Deriva/rafaela/"
nc_files = glob.glob(os.path.join(directory, "05_roms_his_redo.nc")) #your roms file

for file in nc_files:
	filename = os.path.basename(file)
	nc = xr.open_dataset(file)
	num_days = int((np.max(nc.ocean_time) - np.min(nc.ocean_time)) / np.timedelta64(1, 'D'))
	
	o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information
	reader = reader_ROMS_lhico.Reader(file)
	o.add_reader(reader)
	o.set_config('general:coastline_action', 'stranding')
	o.set_config('general:use_auto_landmask', False)
	#o.set_config('seed:wind_drift_factor', 0.02)
	o.set_config('drift:current_uncertainty_uniform',0.01)
	
	var_out = var_out = ['latitude','longitude','age_seconds','x_sea_water_velocity','y_sea_water_velocity']
	
	
	lons = [-50.05581349, 
		-50.03813582, 
	 	-50.02045815, 
	 	-50.00278048, 
	 	-49.98510281, 
		-49.96742514,
		-49.94974747, 
		-49.9320698 , 
		-49.91439213, 
		-49.89671446, 
		-49.87903679, 
		-49.86135912, 
		-49.84368145, 
		-49.82600378, 
		-49.80832611, 
		-49.79064844, 
		-49.77297077, 
		-49.7552931 , 
		-49.73761543,
		-49.71993776, 
		-49.70226009, 
		-49.68458243, 
		-49.54316107, 
		-49.5254834 , 
		-49.50780573, 
		-49.49012806, 
		-49.47245039, 
		-49.45477272, 
		-49.43709505, 
		-49.41941738, 
		-49.40173971, 
		-49.38406204, 
		-49.36638437, 
		-49.3487067 , 
		-49.33102903,
		-49.31335137, 
		-49.2956737 , 
		-49.27799603,
		-49.26031836,
		-49.24264069,
		-49.22496302,
		-49.20728535]
	
	lats = [0.76344065,
    	0.74576298,
     	0.72808531,
     	0.71040764,
     	0.69272997,
    	0.6750523,
    	0.65737463,
    	0.63969696,
    	0.62201929,
    	0.60434162,
    	0.58666395,
    	0.56898628,
    	0.55130861,
    	0.53363094,
    	0.51595328,
    	0.49827561,
    	0.48059794,
    	0.46292027,
    	0.4452426,
    	0.42756493,
    	0.40988726,
    	0.39220959,
    	0.25078823,
    	0.23311056,
    	0.21543289,
    	0.19775522,
    	0.18007755,
    	0.16239988,
    	0.14472222,
    	0.12704455,
    	0.10936688,
    	0.09168921,
    	0.07401154,
    	0.05633387,
    	0.0386562,
    	0.02097853,
    	0.00330086,
    	-0.01437681,
    	-0.03205448,
    	-0.04973215,
    	-0.06740982,
    	-0.08508749]
	
	#scheduling run
	time_step = 60 #1 minutes
	hours_to_run = num_days*24
	hours_to_minutes = (3600/time_step)
	steps =  hours_to_run * hours_to_minutes
	#defining saving
	time_step_output = 60 * time_step #saving hourly
	start_time = reader.start_time + timedelta(hours=1)
	end_time = start_time + timedelta(hours = hours_to_run)
	
	#seeding
	for s in range(len(lons)):
	
		o.seed_elements(lon=lons[s], lat=lats[s], radius=100, radius_type='uniform',number=500,z=0, time= start_time)
	
	
	o.run(steps=steps,time_step=time_step,time_step_output=time_step_output, outfile = f"OPENDRIFT_{filename}",export_variables=var_out)
	print(o)
	o.plot(filename = f"FIGURE_{filename}.png")


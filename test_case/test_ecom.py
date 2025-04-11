

#####################################
open_dir = '/data1/Deriva/opendrift_lhico/ecom/' #add your opendrift-eulerian-model directory

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

#5)  Importando resultado hidrodinâmico
#ecomfile_stokes_path = '/data1/Deriva/opendrift_lhico/opendrift_ecom/ecom_stokes_19_03.nc' #add your ecom file for SES, SBB or Cananeia
#nc = xr.open_dataset(ecomfile_stokes_path)
##reader importation
#reader_pcse = reader_ECOM_S2Z_many_grids.Reader(ecomfile_stokes_path)

ecom_file = 'ses_jan_2020_crop.nc'
nc = xr.open_dataset(ecom_file)
#reader importation
reader_pcse = reader_ECOM_S2Z_many_grids.Reader(ecom_file)


#6) Aplicar o reader ao arquivo
o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information
o.add_reader(reader_pcse)


#7) Configurações livres do OceanDrift
o.set_config('general:coastline_action', 'stranding')
o.set_config('drift:current_uncertainty_uniform', 0.05)	
o.set_config('general:use_auto_landmask', False) #selecionar a máscara do modelo

#8) Criando a rodada
#variáveis para salvamento
var_out = ['latitude', 'longitude', 'age_seconds','z']
output_name = 'test_1.nc'
#8.1) 
time_step = 300  # 5min  in seconds
duration_ = timedelta(hours = 24) #quanto tempo rodará
# Define output saving interval (every 24 hours)
time_step_output = 86400/4 # in seconds
#9) quando começar
start_time = reader_pcse.start_time + timedelta(hours=5*24)

#10) Onde lançar partículas:
## --- locations --###
lon_seed = [-46.30653,-46.37754,-46.345687,-46.339669,-46.329811,-46.323199,-46.317632,-46.312012]
lat_seed = [-23.99264,-23.9736,-23.971564,-23.971990,-23.974045,-23.976469,-23.980284,-23.985568]
#t_seed = ["INF1","INF2","C1","C2","C3","C4","C5", "C6"]

# Seeding particles only where depth >= 180 meters
#for lon, lat, d in zip(lon_seed, lon_seed, z0):
	#o.seed_elements(lon=lon, lat=lat, radius=20, radius_type='uniform', number=21, z=z0, time=start_time)


n_ = 1000 #número de partículas
o.seed_elements(lon=lon_seed[5], lat=lat_seed[5], radius=20, radius_type='uniform', number=n_, z=-5, time=start_time)


o.run(duration=duration_, time_step=time_step, time_step_output=time_step_output)#,outfile=output_name, export_variables=var_out)

o.plot(filename='fig1.png')

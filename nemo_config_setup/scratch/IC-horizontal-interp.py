import numpy as np
import xarray as xr
import gsw
import sys
sys.path.append('/home/users/birgal/')
from nemo_python.interpolation import interp_latlon_cf, neighbours, neighbours_z, extend_into_mask
from nemo_python.utils import fix_lon_range
from tqdm import tqdm

# File locations on Jasmin:
folder_SOSE = '/gws/nopw/j04/terrafirma/birgal/NEMO_AIS/B-SOSE/climatology/'
folder_NEMO = '/gws/nopw/j04/terrafirma/birgal/NEMO_AIS/bathymetry/'

# Load files:
nemo = xr.open_dataset(f'{folder_NEMO}coordinates_AIS.nc').squeeze()
sose = xr.open_dataset(f'{folder_SOSE}SALT_climatology_m01.nc')

# Helper function to convert an xarray dataset with 3D T and S to TEOS10
# - inputs: practical salinity and potential temperature
# - returns: absolute salinity and conservative temperature
def convert_to_teos10(dataset):
    # Convert to TEOS10

    # Need 3D lat, lon, pressure at every point, so if 1D or 2D, broadcast to 3D
    if dataset.lon.values.ndim <= 2:
        lon   = xr.broadcast(dataset['lon'], dataset['potT'])[0]
    if dataset.lat.values.ndim <= 2:
        lat   = xr.broadcast(dataset['lat'], dataset['potT'])[0]
    if dataset.depth.values.ndim <= 2:
        # Need pressure in dbar at every 3D point: approx depth in m
        press = np.abs(xr.broadcast(dataset['depth'], dataset['potT'])[0])
    else:
        press = np.abs(dataset['depth'])
    
    # Get absolute salinity
    absS  = gsw.SA_from_SP(dataset['pracS'], press, lon, lat)
    
    # Get conservative temperature
    consT  = gsw.CT_from_t(absS, dataset['potT'], press)

    return absS, consT

# input:
# - source_file = '/home/birgal/Documents/antarctic/data/B-SOSE/SALT_climatology_m01.nc'
# - ndim=2 or 3 --> tells whether need to loop over depths
# - convert_TEOS10=False
salt_file = f'{folder_SOSE}SALT_climatology_m01.nc'
temp_file = f'{folder_SOSE}THETA_climatology_m01.nc'
dataset   = 'SOSE'

### Horizontal interpolation
datasets = []
# 
print('Horizontally interpolating each depth level in source dataset to NEMO grid')
# Loop over all SOSE depth levels:
for dl in tqdm(range(sose.Z.size)):
    if dataset == 'SOSE':
        name_remapping = {'XC':'lon', 'YC':'lat', 'Z':'depth'}
            
        # Read temperature and salinity from January as a single dataset, and slice to latitude range that I want to interpolate to reduce size
        SOSE_salt = xr.open_dataset(f'{salt_file}').rename(name_remapping).sel(lat=slice(-90, -48))
        SOSE_temp = xr.open_dataset(f'{temp_file}').rename(name_remapping).sel(lat=slice(-90, -48))
        
        # For SOSE, convert longitudes from 0-360 to -180 to 180 
        SOSE_temp['lon'] = fix_lon_range(SOSE_temp.lon)
        SOSE_salt['lon'] = fix_lon_range(SOSE_salt.lon)
        SOSE_temp = SOSE_temp.sortby('lon')
        SOSE_salt = SOSE_salt.sortby('lon')

        # convert to TEOS10:
        SOSE = xr.Dataset({'lon':SOSE_temp['lon'], 'lat':SOSE_temp['lat'], 'depth': SOSE_temp['depth'],
                           'potT':SOSE_temp['THETA'], 'pracS':SOSE_salt['SALT']})
        SOSE_AS, SOSE_CT = convert_to_teos10(SOSE)

        #theta_source = xr.where(SOSE_temp.maskC.isel(depth=dl)==1, SOSE_CT.isel(depth=dl), np.nan)
        #theta_source = xr.where(theta_source==0, np.nan, theta_source)
        salt_source  = xr.where(SOSE_salt.maskC.isel(depth=dl)==1, SOSE_AS.isel(depth=dl), np.nan)
        salt_source  = xr.where(salt_source==0, np.nan, salt_source)
        # Now wrap up into a new Dataset
        source = xr.Dataset({'lon':SOSE['lon'], 'lat':SOSE['lat'], 'SALT': salt_source}) 
        # 'salt':SOSE_SA.isel(depth=0)})
        
        # Interpolate slices of depth levels along lat-lon (horizontally)
        interp_src = interp_latlon_cf(source, nemo, pster_src=False, periodic_src=True, periodic_nemo=True, method='conservative')
        
        datasets.append(interp_src)

SOSE_interpolated = xr.concat(datasets, dim='z').assign_coords(z=np.abs(SOSE.depth.values[0:dl+1]))

SOSE_interpolated.to_netcdf('/gws/nopw/j04/terrafirma/birgal/NEMO_AIS/initial-conditions/salt-horizontal-interp.nc')


import xarray as xr
import numpy as np
import cmocean
import matplotlib.pyplot as plt
import glob
import sys
sys.path.append('/home/users/birgal/')
from nemo_python_git.ocean_calcs import convective_resistance

# Input parameters:
try:
   run_dir = sys.argv[1]
except:
   raise Exception('Please pass a run folder name')

base_dir = '/gws/nopw/j04/anthrofail/birgal/NEMO_AIS/'
nemo_mesh= f'{base_dir}bathymetry/mesh_mask-20240305.nc'
weddell_overall = {'x': slice(850, 1200, None), 'y': slice(50, -1, None)}

gridT_files  = np.sort(glob.glob(f'{run_dir}files/*grid_T*'))
icemod_files = np.sort(glob.glob(f'{run_dir}files/*icemod*'))
dsT     = xr.open_mfdataset(gridT_files).rename({'y_grid_T':'y', 'x_grid_T':'x'}).isel(**weddell_overall)
ice     = xr.open_mfdataset(icemod_files).isel(**weddell_overall)
ds_mesh = xr.open_dataset(nemo_mesh).squeeze().isel(**weddell_overall).rename({'nav_lev':'deptht'})

CR = convective_resistance(dsT, ds_mesh, H=600)
CR = CR.where(CR!=0)

x_inds = [110, 150, 182, 250]
y_inds = [220, 275, 190, 280]
colors = ['#ff6f69','#96ceb4','#ffcc5c','#5c8fff']

sea_ice_area    = ice.siconc * dsT.area_grid_T
annual_seaice_area  = sea_ice_area.sum(['y','x']).groupby('time_counter.year')

def correlate_ice_CR(CR, annual_seaice_area, run_dir):   

    fig, ax = plt.subplots(1,1, figsize=(8, 5), dpi=100)
    for i, (xind, yind) in enumerate(zip(x_inds[0:3], y_inds[0:3])):
        print(i)
        yearly_min = CR.isel(x=xind, y=yind).groupby('time_counter.year').min()
        ax.scatter(annual_seaice_area.min()*1e-12, yearly_min, c=colors[i])
    ax.set_xlabel('Annual minimum sea ice extent (millions of km2)')
    ax.set_ylabel('Annual minimum convective resistance')    
 
    fig.savefig(f'{run_dir}figures/convective_resistance_corr_sea_ice.jpg') 

    return

def spatial_viz(dsT, CR, run_dir):

    fig, ax = plt.subplots(1,2, figsize=(16,4), dpi=100, gridspec_kw={'width_ratios':[1, 2.5]})

    cm = ax[0].pcolormesh(dsT.x, dsT.y, CR.where(CR!=0).isel(time_counter=4), cmap=cmocean.cm.turbid, vmin=0, vmax=0.5)
    ax[0].scatter(x_inds, y_inds, c=colors, marker='*', s=70)
    ax[0].set_ylim(50,350)
    ax[0].set_title(f'To depth: 600 m')
    fig.colorbar(cm, label='Convective resistance, m2/s2', extend='both', ax=ax[0])

    for i, (xind, yind) in enumerate(zip(x_inds, y_inds)):
        yearly_min = CR.isel(x=xind, y=yind).groupby('time_counter.year').min()
        yearly_max = CR.isel(x=xind, y=yind).groupby('time_counter.year').max()
        ax[1].plot(CR.isel(x=xind, y=yind).time_counter, CR.isel(x=xind, y=yind), c=colors[i], linewidth=0.4)
        ax[1].plot((yearly_min.time_counter - np.timedelta64(334,'D')), yearly_min, c=colors[i], linewidth=2, linestyle='-')
        ax[1].plot((yearly_max.time_counter - np.timedelta64(126,'D')), yearly_max, c=colors[i], linewidth=2, linestyle='--')

    ax[1].set_ylim(0, 0.8)
    ax[1].set_ylabel('Convective resistance, m2/s2')

    fig.savefig(f'{run_dir}figures/convective_resistance_spatial.jpg') 

    return 

#spatial_viz(dsT, CR, run_dir)
correlate_ice_CR(CR, annual_seaice_area, run_dir)

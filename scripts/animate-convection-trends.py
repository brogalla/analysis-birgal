import xarray as xr
import numpy as np
import cmocean
import matplotlib.pyplot as plt
import glob
import tqdm
import sys
sys.path.append('/home/users/birgal/')
from nemo_python_git.plots import create_animation
from datetime import datetime

# Input parameters:
try:
   run_dir = sys.argv[1]
except:
   raise Exception('Please pass a run folder name')

base_dir   = '/gws/nopw/j04/anthrofail/birgal/NEMO_AIS/'
nemo_bathy = xr.open_dataset(f'{base_dir}bathymetry/mesh_mask-20240305.nc').squeeze().bathy_metry

for varname1, varname2 in [('thetao', 'ConsTemp')]: # [('so', 'AbsSal')]:#, ('thetao', 'ConsTemp')]:
    print(varname1, varname2)
    # varname1 is the variable name in my simulation file, varname2 is the variable name in the initial conditions file

    # load initial conditions file for comparison
    ini = xr.open_dataset(f'{base_dir}initial-conditions/WOA23-{varname2}-initial-conditions-20241023.nc').squeeze()
    ini[varname2] = xr.where(abs(ini[varname2] > 1e3), np.nan, ini[varname2])

    # check which files have already been analyzed in the past and only loop over the files that still need to be done
    gridT_files = glob.glob(f'{run_dir}files/*grid_T*')
    done_files  = np.sort(glob.glob(f'{run_dir}animations/frames/intermediate_trend_{varname2}_*.png'))
    if len(done_files)!=0:
        last_date   = datetime.strptime(done_files[-1].split(f'{varname2}_')[1].split('.png')[0], '%Y-%m-%d').year
        for f, fileT in enumerate(gridT_files):
            if datetime.strptime(fileT.split('1d_')[1].split('0101_')[0], '%Y').year == last_date:
                print(f'starting from year {last_date}')
                gridT_files = gridT_files[f:]
                break
    for file in tqdm.tqdm(gridT_files):
        sim = xr.open_dataset(file).rename({'y_grid_T':'y', 'x_grid_T':'x'})
    
        # calculate the difference between the surface value and the average over 500-2000 m depth
        surface_sim      = sim[varname1].isel(deptht=0)
        intermediate_sim = sim[varname1].sel(deptht=slice(500,2000)).mean(dim='deptht')
        surface_ini      = ini[varname2].isel(deptht=0)
        intermediate_ini = ini[varname2].sel(deptht=(sim.deptht < 2000)*(sim.deptht>500)).mean(dim='deptht') 
    
        # mask regions where the bathymetry depth is less than 2000 m
        intermediate_sim = xr.where(nemo_bathy < 2000, np.nan, intermediate_sim)
        intermediate_ini = xr.where(nemo_bathy < 2000, np.nan, intermediate_ini)
        # and calculate the anomaly relative to the initial conditions
        surface_anomaly      = surface_sim - surface_ini
        intermediate_anomaly = intermediate_sim - intermediate_ini
    
        # create figure for each time in results file
        for time in range(0,surface_sim.time_counter.size):
            fig1, ax1 = plt.subplots(1,2, figsize=(18,5), dpi=300)
            fig2, ax2 = plt.subplots(1,2, figsize=(18,5), dpi=300)
            
            if varname2=='AbsSal':
                kwags1a = {'vmin':32, 'vmax':35, 'ylim':(150,410), 'xlim':(820,None)}
                kwags1b = {'vmin':34.6, 'vmax':34.9, 'ylim':(150,410), 'xlim':(820,None)}
            else:
                kwags1a = {'vmin':-2, 'vmax':5 , 'ylim':(150,410), 'xlim':(820,None)}
                kwags1b = {'vmin':-1, 'vmax':3 , 'ylim':(150,410), 'xlim':(820,None)}
            kwags2a = {'vmin':-1.5,'vmax':1.5,'cmap':cmocean.cm.balance, 'ylim':(150,410), 'xlim':(820,None)}
            kwags2b = {'vmin':-0.05,'vmax':0.05,'cmap':cmocean.cm.balance, 'ylim':(150,410), 'xlim':(820,None)}

            surface_sim.isel(time_counter=time).plot(ax=ax1[0], **kwags1a)
            surface_anomaly.isel(time_counter=time).plot(ax=ax1[1], **kwags2a)
            intermediate_sim.isel(time_counter=time).plot(ax=ax2[0], **kwags1b)
            intermediate_anomaly.isel(time_counter=time).plot(ax=ax2[1], **kwags2b)
            
            date = np.datetime_as_string(surface_sim.isel(time_counter=time).time_counter, unit='D')
            for ax in [ax1, ax2]:
                ax[0].set_title(f'{varname2} in simulation, {date}')
                ax[1].set_title(f'{varname2} anomaly from IC, {date}')
                
            fig1.set_tight_layout(True); fig2.set_tight_layout(True);
            
            fig1.savefig(f'{run_dir}animations/frames/surface_trend_{varname2}_{date}.png', facecolor='white', transparent=False, format='png')
            fig2.savefig(f'{run_dir}animations/frames/intermediate_trend_{varname2}_{date}.png', facecolor='white', transparent=False, format='png')
            plt.close(fig1); plt.close(fig2);

create_animation(glob.glob(f'{run_dir}animations/frames/surface_trend_AbsS*'), out_file=f'{run_dir}animations/animation_surface_trend_AbsSal.mp4')
#create_animation(glob.glob(f'{run_dir}animations/frames/surface_trend_ConsT*'), out_file=f'{run_dir}animations/animation_surface_trend_ConsTemp.mp4')
create_animation(glob.glob(f'{run_dir}animations/frames/intermediate_trend_AbsS*'), out_file=f'{run_dir}animations/animation_intermediate_trend_AbsSal.mp4')
#create_animation(glob.glob(f'{run_dir}animations/frames/intermediate_trend_ConsT*'), out_file=f'{run_dir}animations/animation_intermediate_trend_ConsTemp.mp4')



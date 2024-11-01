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

for varname1, varname2 in [('so', 'AbsSal'), ('thetao', 'ConsTemp')]:
    print(varname1, varname2)
    # varname1 is the variable name in my simulation file, varname2 is the variable name in the initial conditions file

    # load initial conditions file for comparison
    ini = xr.open_dataset(f'{base_dir}initial-conditions/WOA23-{varname2}-initial-conditions-20241023.nc').squeeze()
    ini[varname2] = xr.where(abs(ini[varname2] > 1e3), np.nan, ini[varname2])
    
    gridT_files = glob.glob(f'{run_dir}files/*grid_T*')
    done_files  = np.sort(glob.glob(f'{run_dir}animations/frames/delta_{varname2}_*.png'))
    if len(done_files)!=0:
        last_date   = datetime.strptime(done_files[-1].split(f'{varname2}_')[1].split('.png')[0], '%Y-%m-%d').year
        for f, fileT in enumerate(gridT_files):
            if datetime.strptime(fileT.split('1m_')[1].split('0101_')[0], '%Y').year == last_date:
                print(f'starting from year {last_date}')
                gridT_files = gridT_files[f:]
                break

    for file in tqdm.tqdm(gridT_files):
        sim = xr.open_dataset(file).rename({'y_grid_T':'y', 'x_grid_T':'x'})
    
        # calculate the difference between the surface value and the average over 500-2000 m depth
        delta_sim = (sim[varname1].sel(deptht=slice(500,2000)).mean(dim='deptht') - sim[varname1].isel(deptht=0))
        delta_ini = (ini[varname2].sel(deptht=(sim.deptht < 2000)*(sim.deptht>500)).mean(dim='deptht') - ini[varname2].isel(deptht=0))
    
        # mask regions where the bathymetry depth is less than 2000 m
        delta_sim = xr.where(nemo_bathy < 2000, np.nan, delta_sim)
        delta_ini = xr.where(nemo_bathy < 2000, np.nan, delta_ini)
        # and calculate the anomaly relative to the initial conditions
        delta_anomaly = delta_sim - delta_ini
    
        # create figure for each time in results file
        for time in range(0,delta_sim.time_counter.size):
            fig, ax = plt.subplots(1,2, figsize=(17,4), dpi=300)
            
            kwags = {'vmin':-3,'vmax':3,'cmap':cmocean.cm.balance, 'ylim':(220,410), 'xlim':(820,None) } # 200, None
            delta_sim.isel(time_counter=time).plot(ax=ax[0], **kwags)
            delta_anomaly.isel(time_counter=time).plot(ax=ax[1], **kwags)
            
            date = np.datetime_as_string(delta_sim.isel(time_counter=time).time_counter, unit='D')
            ax[0].set_title(f'$\delta${varname2} (500-2000 m minus surface) in simulation, {date}')
            ax[1].set_title(f'$\delta${varname2} (500-2000 m minus surface) anomaly from IC, {date}')
            for axis in ax.ravel():
                axis.set_xlabel(''); axis.set_ylabel('');
            fig.set_tight_layout(True)
            fig.savefig(f'{run_dir}animations/frames/delta_{varname2}_{date}.png', facecolor='white', transparent=False, format='png')
            plt.close(fig)

create_animation(glob.glob(f'{run_dir}animations/frames/delta_AbsS*'), out_file=f'{run_dir}animations/animation_delta_AbsSal.mp4')
create_animation(glob.glob(f'{run_dir}animations/frames/delta_ConsT*'), out_file=f'{run_dir}animations/animation_delta_ConsTemp.mp4')


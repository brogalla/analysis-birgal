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

varname = 'mldr10_1'
    
gridT_files = glob.glob(f'{run_dir}files/*grid_T*')
done_files  = np.sort(glob.glob(f'{run_dir}animations/frames/mld_*.png'))
if len(done_files)!=0:
   last_date   = datetime.strptime(done_files[-1].split(f'mld_')[1].split('.png')[0], '%Y-%m-%d').year
   print(last_date)
   for f, fileT in enumerate(gridT_files):
       if datetime.strptime(fileT.split('1d_')[1].split('0101_')[0], '%Y').year == last_date:
           print(f'starting from year {last_date}')
           gridT_files = gridT_files[f:]
           break

for file in tqdm.tqdm(gridT_files):
   sim = xr.open_dataset(file).rename({'y_grid_T':'y', 'x_grid_T':'x'})
    
   # create figure for each time in results file
   for time in range(0,sim.time_counter.size):
       fig, ax = plt.subplots(1,1, figsize=(8,4), dpi=300)
            
       kwags = {'vmin':0,'vmax':600,'cmap':cmocean.cm.deep, 'ylim':(150,410), 'xlim':(820,None) } # 200, None
       sim[varname] = xr.where(sim[varname]==0, np.nan, sim[varname])
       sim[varname].isel(time_counter=time).plot(ax=ax, **kwags)
            
       date = np.datetime_as_string(sim.isel(time_counter=time).time_counter, unit='D')
       ax.set_title(date)
       ax.set_xlabel(''); ax.set_ylabel('');
       fig.set_tight_layout(True)
       fig.savefig(f'{run_dir}animations/frames/mld_{date}.png', facecolor='white', transparent=False, format='png')
       plt.close(fig)

create_animation(glob.glob(f'{run_dir}animations/frames/mld_*'), out_file=f'{run_dir}animations/animation_MLD.mp4')


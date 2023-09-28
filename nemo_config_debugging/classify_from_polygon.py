import numpy as np
import xarray as xr
import pyproj
from multiprocessing import Process, Manager
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from shapely.ops import transform
from joblib import Parallel, delayed

def create_transformation(polygon):
    # transform  points from polygon into x,y coordinates, so that the .contains method works properly.
    wgs84 = pyproj.CRS('EPSG:4326')
    utm = pyproj.CRS('EPSG:32618')

    # from lat lon to UTM
    project       = pyproj.Transformer.from_crs(wgs84, utm, always_xy=True).transform
    utm_polygon   = transform(project, polygon)

    # # and back to lat lon
    # project_back  = pyproj.Transformer.from_crs(utm, wgs84, always_xy=True).transform
    # polar_isobath = transform(project_back, utm_isobath)

    return project, utm_polygon

def check_contain(x_points, y_points, project, utm_polygon):  # the managed list `L` passed explicitly.
    # Function takes x_points array of longitude values and y_points array of latitudes, and checks whether the coordinate 
    # pairs fall within the 2000 m Anarctic isobath 
    
    # Make coordinate pairs
    coords = np.array(list(zip(x_points, y_points)))
    
    # Create a mask where true is points within the 2000 m isobath, and false is outside of the isobath
    contains  = np.vectorize(lambda p: utm_polygon.contains(transform(project, Point(p))), signature='(n)->()')
    
    # Write booleans to a list that is shared between processes
    return list(contains(coords))

def check_contain_parallel(x_points, y_points, polygon):  # the managed list `L` passed explicitly.
    # Function takes x_points array of longitude values and y_points array of latitudes, and checks whether the coordinate 
    # pairs fall within the 2000 m Anarctic isobath 
    project, utm_polygon = create_transformation(polygon)
    
    # Make coordinate pairs
    coords = np.array(list(zip(x_points, y_points)))
    
    # Create a mask where true is points within the 2000 m isobath, and false is outside of the isobath
    contains  = np.vectorize(lambda p: utm_polygon.contains(transform(project, Point(p))), signature='(n)->()')
    
    # Write booleans to a list that is shared between processes
    return list(contains(coords))

def coord_in_polygon(x_coord, y_coord, polygon, n_chunks=16):
    # ---- Main function ----
    # Uses n_chunks chunks to check if x_coord array of longitudes and y_coord array of latitudes are
    # within the provided shapely polygon
    
    project, utm_polygon = create_transformation(polygon)
    
    i_start=0; i_end=0;
    boolean_array = []
    for i in range(n_chunks):
        if i==n_chunks:
            print('last chunk')
            i_end = -1
        else:
            i_end += x_coord.shape[0]/n_chunks
        boolean_list = check_contain(x_coord[int(i_start):int(i_end)], y_coord[int(i_start):int(i_end)], \
                                     project, utm_polygon)  # Pass a list of lat lon points to test
        i_start += x_coord.shape[0]/n_chunks
        
        boolean_array = boolean_array + boolean_list

    np_boolean = np.array(boolean_array)
                    
    return np_boolean

def coord_in_polygon_parallel(x_coord, y_coord, polygon, n_chunks=16):
    # ---- Main function ----
    # Uses n_chunks chunks to check if x_coord array of longitudes and y_coord array of latitudes are
    # within the provided shapely polygond
    
    i_start=0; i_end=0;
    for i in range(n_chunks):
        i_end += x_coord.shape[0]/n_chunks
        if i==0:
            chunks = slice(int(i_start),int(i_end))
        else:
            chunks = np.append(chunks, slice(int(i_start),int(i_end)))
        i_start += x_coord.shape[0]/n_chunks
         
    boolean = Parallel(n_jobs=n_chunks, backend= 'multiprocessing')(delayed(check_contain_parallel)(x_coord[chunk], y_coord[chunk], polygon) for chunk in chunks)
    
    return np.array(boolean)


# def check_contain(L, x_points, y_points, project, utm_polygon):  # the managed list `L` passed explicitly.
#     # Function takes x_points array of longitude values and y_points array of latitudes, and checks whether the coordinate 
#     # pairs fall within the 2000 m Anarctic isobath 
    
#     # Make coordinate pairs
#     coords = np.array(list(zip(x_points, y_points)))
    
#     # Create a mask where true is points within the 2000 m isobath, and false is outside of the isobath
#     contains  = np.vectorize(lambda p: utm_polygon.contains(transform(project, Point(p))), signature='(n)->()')
    
#     # Write booleans to a list that is shared between processes
#     L.append(list(contains(coords)))

# def coord_in_polygon(x_coord, y_coord, polygon, n_chunks=16):
#     # ---- Main function ----
#     # Uses n_chunks processors to check if x_coord array of longitudes and y_coord array of latitudes are
#     # within the provided shapely polygon
#     # Only works in chunks mutiples of 16
    
#     project, utm_polygon = create_transformations(polygon)
    
#     with Manager() as manager:
#         L = manager.list()  # <-- list can be shared between processes.
#         processes = []
            
#         i_start=0; i_end=0;
#         for i in range(n_chunks):
#             i_end += x_coord.shape[0]/n_chunks
#             p = Process(target=check_contain, args=(L, x_coord[int(i_start):int(i_end)],
#                                                        y_coord[int(i_start):int(i_end)],
#                                                        project, utm_polygon))  # Pass a list of lat lon points to test
#             i_start += x_coord.shape[0]/n_chunks
#             p.start()
#             processes.append(p)
            
#         for p in processes:
#             p.join()
    
#         boolean_array = np.array(L, dtype='object')

#         # Convert to a more manageable numpy array format. I'm sure there's a better way to write the shared processes to an array instead 
#         # of a list to avoid this step, but I don't immediately know it...
#         boolean_array_total = boolean_array[0]
#         for i in range(1,n_chunks):
#             boolean_array_total = boolean_array_total + boolean_array[i]
    
#             np_boolean = np.array(boolean_array_total)
                    
#         return np_boolean
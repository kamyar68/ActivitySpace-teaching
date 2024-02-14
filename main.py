# IREM tool in free Python
import pathlib
import geopandas as gpd
#from osgeo import gdal, osr, ogr
from rasterio.transform import from_origin
import rasterio
from rasterio.enums import Resampling
from rasterio import Affine
from shapely.geometry import Point, MultiPoint, MultiPolygon, linestring, Polygon
import networkx as nx
import osmnx as ox
from shapely.geometry import mapping
import warnings
from rasterio.mask import mask
import tempfile
import shapely
from shapely.geometry import Point, LineString
import matplotlib.pyplot as plt
from functions import generate_points_along_routes, generate_points_along_polygon
import math
import pandas as pd
from scipy.interpolate import griddata
import numpy as np
import os
import tempfile
from rasterio.enums import Resampling
from rasterio.windows import Window
from rasterio.features import geometry_mask
warnings.filterwarnings('ignore')

NOTEBOOK_PATH = pathlib.Path().resolve()
DATA_DIRECTORY = NOTEBOOK_PATH / "data"
output_DIRECTORY = NOTEBOOK_PATH / "output"
temp_DIRECTORY = NOTEBOOK_PATH / "output" / "temp"
raster_DIRECTORY = NOTEBOOK_PATH / "output/raster"

# Set your parameters
#home_file = DATA_DIRECTORY / "home20.shp"
#print (home_file)
#activity_file = DATA_DIRECTORY / "eep20.shp"
output_file = output_DIRECTORY /"Output.gpkg"
#route_file = DATA_DIRECTORY /"shortest_routes.shp"
uid_field = 'uid'  # Replace with your actual uid field name
D1 = 500  # Replace with your desired buffer distance for home location
D2 = 200   # Replace with your desired buffer distance for activity points

def nb_noroute(home_file, activity_file, output_file, uid_field, D1, D2):
    # Read home and activity shapefiles
    gdf_home = gpd.read_file(home_file)
    gdf_activity = gpd.read_file(activity_file)
    gdf_activity=gdf_activity.to_crs("EPSG:3067")
    gdf_home=gdf_home.to_crs("EPSG:3067")

    # Merge based on the common uid field
    gdf = gdf_home.merge(gdf_activity, on=uid_field, suffixes=('_home', '_activity'))

    # Iterate over each unique uid
    unique_uids = gdf[uid_field].unique()
    convex_hulls = []

    for uid in unique_uids:
        # Select data for the current uid
        current_data = gdf[gdf[uid_field] == uid]

        # Buffer activity points
        current_data['geometry'] = current_data['geometry_activity'].apply(lambda x: x.buffer(D2))
        current_data1 = gpd.GeoDataFrame(current_data)
        # Buffer home location
        current_data['geometry'] = current_data['geometry_home'].apply(lambda x: x.buffer(D1))

        current_data2= gpd.GeoDataFrame(current_data)
        # Merge the two buffered GeoDataFrames
        print(current_data)
        current_data_union = gpd.overlay(current_data1, current_data2, how='union')

        # Create convex hull
        convex_hull = current_data_union.unary_union.convex_hull

        # Append to the list
        convex_hulls.append({'uid': uid, 'geometry': convex_hull})

    # Create GeoDataFrame from the list of convex hulls
    gdf_result = gpd.GeoDataFrame(convex_hulls, geometry='geometry')
    gdf_result.crs=gdf_home.crs

    # Save the result to a new shapefile
    gdf_result.to_file(output_file)


def nb_withRoute(home_file, activity_file, routes_file, output_file, uid_field, D1, D2):
    # Read home and activity shapefiles
    gdf_home = gpd.read_file(home_file)
    gdf_activity = gpd.read_file(activity_file)
    gdf_routes = gpd.read_file(routes_file)
    gdf_activity=gdf_activity.to_crs("EPSG:3067")
    gdf_home=gdf_home.to_crs("EPSG:3067")
    gdf_routes=gdf_routes.to_crs("EPSG:3067")

    # Merge based on the common uid field
    gdf = gdf_home.merge(gdf_activity, on=uid_field, suffixes=('_home', '_activity'))

    # Iterate over each unique uid
    unique_uids = gdf[uid_field].unique()
    convex_hulls = []

    for uid in unique_uids:
        # Select data for the current uid
        current_data = gdf[gdf[uid_field] == uid]

        # Buffer activity points
        current_data['geometry'] = current_data['geometry_activity'].apply(lambda x: x.buffer(D2))
        current_data1 = gpd.GeoDataFrame(current_data)

        # Buffer home location
        current_data['geometry'] = current_data['geometry_home'].apply(lambda x: x.buffer(D1))
        current_data2 = gpd.GeoDataFrame(current_data)

        # Merge the two buffered GeoDataFrames
        current_data_union = gpd.overlay(current_data1, current_data2, how='union')

        # Buffer routes and include them in the union
        routes_for_uid = gdf_routes[gdf_routes[uid_field] == uid]
        routes_buffered = routes_for_uid.copy()
        routes_buffered['geometry'] = routes_for_uid['geometry'].buffer(10)  # Buffer routes by 10 meters

        current_data_union_with_routes = gpd.overlay(current_data_union, routes_buffered, how='union')

        # Create convex hull
        convex_hull = current_data_union_with_routes.unary_union.convex_hull

        # Append to the list
        convex_hulls.append({'uid': uid, 'geometry': convex_hull})

    # Create GeoDataFrame from the list of convex hulls
    gdf_result = gpd.GeoDataFrame(convex_hulls, geometry='geometry')
    gdf_result.crs = gdf_home.crs

    # Save the result to a new shapefile
    gdf_result.to_file(output_file)

#def IREM(home_file, activity_file, route_file, output_file,)

#nb_with_routes_osm(home_file, activity_file, output_file, uid_field, D1, D2)
#nb_withRoute(home_file, activity_file, route_file, output_file, uid_field, D1, D2)
#nb_noroute(home_file, activity_file, output_file, uid_field, D1, D2)
def IREM(home_file,activity_file,output_file,route_file, pixel_size):
    maxw = 30  # Replace with your actual maxw value
    poi = gpd.read_file(activity_file)
    nb= gpd.read_file(output_file)
    routes = gpd.read_file(route_file)
    home = gpd.read_file(home_file)
    poi = poi.to_crs('EPSG:3067')
    routes = routes.to_crs('EPSG:3067')
    home = home.to_crs('EPSG:3067')
    
    all_files = os.listdir(raster_DIRECTORY)
    
    # Filter files with a .tif extension
    existingrasters = [file for file in all_files if file.lower().endswith('.tif')]
    existinguidlist=[]
    for inRaster in existingrasters:
        uid1=inRaster .rsplit('_', 1)[1]
        uid=int(uid1.rsplit('.', 1)[0])
        existinguidlist.append(uid)
    
    def weight(x):
        x = float(x)
        maxw=30
        o = float(-1 * x / maxw)
        value = 1 / (1 + math.exp(o))
        value = round(value, 3)
        return value
    def calculate_w_road(row):
        mod = row['tmod']
        w1 = row['w']
    
        if mod == 'walking':
            rw = math.sqrt(w1 * w2)
        elif mod == 'bike':
            rw = (math.sqrt(w1 * w2)) / 3.4
        else:
            rw = (math.sqrt(w1 * w2)) / 10.0
    
    
    
        return rw
    w2 = float(1/(1+ math.exp(-1)))
    home['w'] = w2
    poi['w'] = poi['Freq'].apply(lambda x: weight(x))
    joined_paths = gpd.GeoDataFrame(pd.merge(routes, poi[['DESTid', 'Freq', 'tmod', 'w']], on="DESTid"))
    print("Converting routes to points...")
    pointed_paths=generate_points_along_routes(joined_paths, temp_DIRECTORY/"pointed_path.shp", 400)
    print("Converting home ranges to points...")
    pointed_nb=generate_points_along_polygon(gpd.read_file(output_file), temp_DIRECTORY/"pointed_nb.shp", 100)
    print("Points created")
    print("Now calculating weights...")
    pointed_paths['w_road'] = pointed_paths.apply(calculate_w_road, axis=1)
    pointed_paths['w'] = pointed_paths['w_road']
    pointed_nb['w'] = 0.05
    # List of fields to keep
    fields_to_keep = ['uid', 'w', 'geometry']
    
    # Concatenate GeoDataFrames as rows and keep only specified fields
    concatenated_df = pd.concat([pointed_paths[fields_to_keep], home[fields_to_keep], poi[fields_to_keep], pointed_nb[fields_to_keep]], ignore_index=True)
    concatenated_gdf = gpd.GeoDataFrame(concatenated_df, geometry='geometry')
    #print (concatenated_gdf.head())
    #concatenated_gdf.to_file("D:/DeveloperZone/MyWorkScripts/ActivitySpaces/PurePython/output/temp/concat.shp")
    
    for uid in concatenated_gdf['uid'].unique():
    
        print (f"Making IREM raster for uid:{uid}")
        if uid in existinguidlist:
            print ("raster already exists")
        else:
            # Filter data for the current uid
            subset = concatenated_gdf[concatenated_gdf['uid'] == uid]
    
            # Get bounds for the raster
            xmin, ymin, xmax, ymax = subset.total_bounds
    
            # Set pixel size (adjust as needed)
            #pixel_size = 50
    
            # Create an empty raster
            cols = int((xmax - xmin) / pixel_size)
            rows = int((ymax - ymin) / pixel_size)
    
            transform = from_origin(xmin, ymax, pixel_size, pixel_size)
            raster_data = np.zeros((rows, cols))
    
            # Perform IDW interpolation
            for i in range(rows):
                for j in range(cols):
                    x = xmin + j * pixel_size
                    y = ymax - i * pixel_size
                    weights = 1 / np.sqrt(((subset.geometry.x - x) ** 2 + (subset.geometry.y - y) ** 2))
                    raster_data[i, j] = np.sum(weights * subset['w']) / np.sum(weights)
            #raster_filepath = os.path.join(raster_DIRECTORY, f'nb_{uid}.tif')
            raster_filepath = raster_DIRECTORY / f'nb_{uid}.tif'
            print (raster_filepath)
            # Save the raster to a GeoTIFF file
            raster_profile = {
                'driver': 'GTiff',
                'count': 1,
                'dtype': raster_data.dtype,
                'width': cols,
                'height': rows,
                'crs': subset.crs,
                'transform': transform,
            }
    
            with rasterio.open(raster_filepath, 'w', **raster_profile) as dst:
                dst.write(raster_data, 1)
    
    print("IDW rasters created successfully.")
    #print (nb.head())
    # Iterate through each raster file in the folder
    for raster_filename in os.listdir(raster_DIRECTORY):
        if raster_filename.endswith('.tif'):
            # Extract the UID from the raster filename
            uid = int(raster_filename.split('_')[1].split('.')[0])
            print(f"clipping IREM for uid: {uid}")
    
            # Filter the polygon GeoDataFrame for the current UID
            polygon_subset = nb[nb['uid'] == uid]
            #print (polygon_subset)
    
            # Load the raster
            raster_path =  raster_DIRECTORY / f'{raster_filename}'
            with rasterio.open(raster_path) as src:
                # Clip the raster using the polygon
                out_image, out_transform = mask(src, polygon_subset.geometry.apply(mapping), crop=True)
                out_meta = src.meta.copy()
    
            # Update metadata with new transform and dimensions
            out_meta.update({"driver": "GTiff",
                             "height": out_image.shape[1],
                             "width": out_image.shape[2],
                             "transform": out_transform,
                             "nodata": 0})
    
            # Overwrite the original raster with the clipped version
            with rasterio.open(raster_path, "w", **out_meta) as dest:
                dest.write(out_image)
    
    print("IREM is completed")
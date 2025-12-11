import xarray as xr
import pandas as pd
import numpy as np
import glob
import os
import geopandas as gpd
import utm
import mgrs

# ----------------------------- Functions ------------------------------------

def read_pixel_gamma(path):
    # Get date from filename
    filename = os.path.basename(path)
    date_str = filename.split('_')[1]
    dt = pd.to_datetime(date_str, format="%Y%m%dT%H%M%S")
    
    try:
        # Open specified group
        with xr.open_zarr(path, group=group_path, consolidated=False) as ds:
            
            # Check if gamma0 exists in this file
            if 'gamma0' not in ds:
                return None
            
            # Select the pixel and .load() it - drop acquisition date because it causes many NaNs
            pixel_val = ds['gamma0'].drop_vars('acquisition_date').sel(y=northing, x=easting, method="nearest").load()
            
            # Only use file if there is data at the pixel (non-zero)
            vv_value = pixel_val.sel(polarization='V:V').item()
            if vv_value == 0:
                print(f"Skipping {filename}: V:V value at pixel is zero.")
                return None
            
            # Convert to dB
            pixel_val = 10 * np.log10(pixel_val)
            
            # Remove attributes to prevent error
            if '_eopf_attrs' in pixel_val.attrs:
                del pixel_val.attrs['_eopf_attrs']
                
            # Attach the time coordinate
            pixel_val = pixel_val.expand_dims(time=[dt])
            
            print(f"Successfully read {filename}")
            return pixel_val
    except: 
        return None

# ------------------------- End of Functions ---------------------------------

## Configuration        
# Specify path containing zarr directories
zarr_dir = r"\\projectdata.eurac.edu\SNOWCOP\S1_ESA" 
zarr_path = os.path.join(zarr_dir, "*.zarr")
file_paths = sorted(glob.glob(zarr_path))

# Read shapefile & get lon, lat (skip this if you already have lon and lat)
shape_dir = r"\\projectdata.eurac.edu\SNOWCOP\Dati\data-merged\SW_radiation.shp" 
gdf = gpd.read_file(shape_dir)
laguna = gdf[gdf["sta_name"] == "Laguna Negra"].iloc[0]
easting, northing = laguna.geometry.x, laguna.geometry.y
lat, lon = utm.to_latlon(easting, northing, zone_number=19, northern=False) # Specify UTM zone & hemisphere (https://hls.gsfc.nasa.gov/products-description/tiling-system/)
print(f"Point: East {easting} North {northing} | Lon {lon} Lat {lat}")

# Read MGRS tile name based on lon, lat
# lat, lon = 
easting, northing, _, _ = utm.from_latlon(lat, lon)
m = mgrs.MGRS()
tile_id = m.toMGRS(lat, lon, MGRSPrecision=0)
print(f"Tile: {tile_id} | Easting: {easting}, Northing: {northing}")
group_path = f"{tile_id}/measurements"


## Pixel extraction
# Process all the files
pixel_list = [read_pixel_gamma(f) for f in file_paths]

# Filter out Nones
pixel_list = [p for p in pixel_list if p is not None]

if pixel_list:
    # Combine into time series
    ds_timeseries = xr.concat(pixel_list, dim='time').sortby('time')
    ds_timeseries = ds_timeseries.squeeze()
    
    # Save as netcdf
    output_filename = f"output_{tile_id}.nc"
    ds_timeseries.to_netcdf(output_filename)
    
    # Open netcdf
    dataset = xr.open_dataset(output_filename)
    print(dataset)
    print(f"Exact pixel coordinates: Easting: {dataset.coords['x'].values} Northing: {dataset.coords['y'].values}")
else:
    print("No data found, double check your tile_id and coordinates")
        
    

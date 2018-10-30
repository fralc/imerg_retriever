from gpm_utils import *

# EXAMPLE 1:
# retrieve a file from ftp repository (the start date is derived from a filename)
filename = '3B-HHR.MS.MRG.3IMERG.20180104-S000000-E002959.0000.V05B.HDF5'
start_date = get_start_time(filename)
filename = get_imerg(start_date)

# EXAMPLE 2:
# Compare tif and hdf5 version of same IMERG data in order to verify the georeferencing matrix
filename = '3B-HHR.MS.MRG.3IMERG.20180201-S090000-E092959.0540.V05B.HDF5'

# Read the hdf5
datasets = read_hdf5(filename)
# plt.imshow(datasets['IRprecipitation'].T, vmin=0)

# Read the reference tif
src_ds = gdal.Open('3B-HHR-E.MS.MRG.3IMERG.20180201-S090000-E092959.0540.V05B.30min.tif')
myarray = np.array(src_ds.GetRasterBand(1).ReadAsArray(), dtype=np.float64)
geot = src_ds.GetGeoTransform()

# Plot both hdf5 and tif precipitation maps transforming hdf5 map to be congruent with the tiff map.
precip = datasets['precipitationCal']
myarray[myarray == np.max(myarray)] = np.nan
precip[precip == np.min(precip)] = np.nan
plt.subplot(2, 1, 1)
plt.imshow(myarray)
plt.subplot(2, 1, 2)
plt.imshow(np.flipud(precip.T))
plt.show()

# EXAMPLE 3:
# verify extraction function
# Read the reference tif
src_ds = gdal.Open('data/3B-HHR-E.MS.MRG.3IMERG.20180201-S090000-E092959.0540.V05B.30min.tif')
myarray = np.array(src_ds.GetRasterBand(1).ReadAsArray(), dtype=np.float64)
geot = src_ds.GetGeoTransform()

raster = Raster(myarray, geot)
italy_lat_range = (35.49370, 47.09178)
italy_lon_range = (6.62662, 18.52038)
raster_italy = raster.extract(lat_range=italy_lat_range,
                              lon_range=italy_lon_range)
raster_italy.save_as_tiff('italy.tif')

# EXAMPLE 4:
# Read hdf5 as raster and extract italy area only
filename = '3B-HHR.MS.MRG.3IMERG.20180201-S090000-E092959.0540.V05B.HDF5'

# Read the hdf5
GEO_TRANSFORMATION_PARAMS = imerg_info['GEO_TRANSFORMATION_PARAMS']
datasets = read_hdf5(filename)
geot = build_geot(**GEO_TRANSFORMATION_PARAMS)
rasters = datasets_to_rasters(datasets, geot)

italy_lat_range = (35.49370, 47.09178)
italy_lon_range = (6.62662, 18.52038)
rasters_ita = {}
for k, r in rasters.items():
    rasters_ita[k] = r.extract(lat_range=italy_lat_range,
                               lon_range=italy_lon_range)






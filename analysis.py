from gpm_utils import *

TEMPLATE_FILENAME = '3B-HHR-L.MS.MRG.3IMERG.20181001-S000000-E002959.0000.V05B.RT-H5'

# we assume that the geotransform matrix is equal to that available in imerg_info
geot = build_geot(**imerg_info['GEO_TRANSFORMATION_PARAMS'])

# we take the field keys for building filename from a template filename
finfo = filename_info(TEMPLATE_FILENAME)

# retrieving the file relative to a given time-stamp
start_time = datetime.strptime("201810211730", "%Y%m%d%H%M")
# fn = get_imerg_web(finfo, start_time)
fn = '3B-HHR-L.MS.MRG.3IMERG.20181021-S173000-E175959.1050.V05B.RT-H5'

# reading the file
ras = hdf5_to_rasters(os.path.join('data', fn), geot)
# ras['precipitationCal'].save_as_tiff('data/precipCal.tif')
# ras['probabilityLiquidPrecipitation'].save_as_tiff('data/probLiq.tif')

# extracting AOI data for a raster
italy_lat_range = (35.49370, 47.09178)
italy_lon_range = (6.62662, 18.52038)
raster = ras['precipitationCal']
raster_italy = raster.extract(lat_range=italy_lat_range,
                              lon_range=italy_lon_range)
raster_italy.save_as_tiff('data/italy.tif')
print('ok')
#########################################################

# fn = r'C:\Users\LEI00025\gits\arpa\gribs\anomalieTemperatureMassime_0102201800002402.grib'
fn_arpa = r'C:\Users\LEI00025\Downloads\precipitazioni_2110201800000148.grib'
lyrs_arpa = GribData(fn_arpa)

fn_imerg = 'data/3B-HHR-L.MS.MRG.3IMERG.20181021-S173000-E175959.1050.V05B.RT-H5'
lyrs_imerg = H5Data(fn_imerg)



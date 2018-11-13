from gpm_utils import *

TEMPLATE_FILENAME = '3B-HHR-L.MS.MRG.3IMERG.20181001-S000000-E002959.0000.V05B.RT-H5'

# we assume that the geotransform matrix is equal to that available in imerg_info
geot = build_geot(**imerg_info['GEO_TRANSFORMATION_PARAMS'])

# we take the field keys for building filename from a template filename
finfo = filename_info(TEMPLATE_FILENAME)

# retrieving the file relative to a given time-stamp
start_time = datetime.strptime("201811041730", "%Y%m%d%H%M")
fn = get_imerg_web(finfo, start_time)
start_time = datetime.strptime("201806041730", "%Y%m%d%H%M")
fn = get_imerg(start_time)

ras = hdf5_to_rasters(os.path.join('data', fn), geot)
ras['probabilityLiquidPrecipitation']
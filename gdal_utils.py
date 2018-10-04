from osgeo import gdal, ogr, osr
import numpy as np

# web references:
# https://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html
# https://gis.stackexchange.com/questions/57834/how-to-get-raster-corner-coordinates-using-python-gdal-bindings
# https://pcjericks.github.io/py-gdalogr-cookbook/


def raster2array(filename):
    src_ds = gdal.Open(filename)
    raster = np.array(src_ds.GetRasterBand(1).ReadAsArray(), dtype=np.float64)
    return raster


def GetExtent(gt, cols, rows):
    ''' Return list of corner coordinates from a geotransform

        @type gt:   C{tuple/list}
        @param gt: geotransform
        @type cols:   C{int}
        @param cols: number of columns in the dataset
        @type rows:   C{int}
        @param rows: number of rows in the dataset
        @rtype:    C{[float,...,float]}
        @return:   coordinates of each corner
    '''
    ext = []
    xarr = [0, cols]
    yarr = [0, rows]

    for px in xarr:
        for py in yarr:
            x = gt[0] + (px * gt[1]) + (py * gt[2])
            y = gt[3] + (px * gt[4]) + (py * gt[5])
            ext.append([x, y])
            print(x, y)
        yarr.reverse()
    return ext


def ReprojectCoords(coords, src_srs, tgt_srs):
    ''' Reproject a list of x,y coordinates.

        @type geom:     C{tuple/list}
        @param geom:    List of [[x,y],...[x,y]] coordinates
        @type src_srs:  C{osr.SpatialReference}
        @param src_srs: OSR SpatialReference object
        @type tgt_srs:  C{osr.SpatialReference}
        @param tgt_srs: OSR SpatialReference object
        @rtype:         C{tuple/list}
        @return:        List of transformed [[x,y],...[x,y]] coordinates
    '''
    trans_coords = []
    transform = osr.CoordinateTransformation(src_srs, tgt_srs)
    for x, y in coords:
        x, y, z = transform.TransformPoint(x, y)
        trans_coords.append([x, y])
    return trans_coords


def array2raster(newRasterfn, geot, array):

    cols = array.shape[1]
    rows = array.shape[0]

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Byte)
    outRaster.SetGeoTransform(geot)

    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)

    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(4326)

    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()


def pixelCoords_to_mapCoords(x, y, geot):
    # xgeo = geot[0] + x * geot[1] + y * geot[2]
    # ygeo = geot[3] + x * geot[4] + y * geot[5]
    xgeo, ygeo = np.dot(np.reshape(geot, [2, 3]), np.array([1, x, y]).T)
    return xgeo, ygeo


def mapCoords_to_pixelCoords(geot, xgeo=None, ygeo=None):
    if geot[2] != 0 or geot[4] != 0:
        raise NotImplementedError('Not implemented for rotated images.')
    if xgeo is not None and ygeo is not None:
        x = int((xgeo - geot[0]) / geot[1])
        y = int((ygeo - geot[3]) / geot[5])
        return x, y
    elif xgeo is not None:
        x = int((xgeo - geot[0]) / geot[1])
        return x
    elif ygeo is not None:
        y = int((ygeo - geot[3]) / geot[5])
        return y


def build_geot(originX, originY, pixelWidth, pixelHeight):
    geot = (float(originX), float(pixelWidth), 0.,
            float(originY), 0., float(pixelHeight))

    return geot


class Raster(object):
    """
    Instances of this class represent a georeferenced raster along with some methods to
        manage it.
    Args:
        array (numpy.ndarray): 2D array with raster values.
        geot (tuple): GDAL geotransform tuple.
    Attributes:
        array (numpy.ndarray): 2D array with raster values.
        geot (tuple): GDAL geotransform tuple relative to the raster array.
    """
    def __init__(self, array, geot):
        self._array = array
        self._geot = geot

    @property
    def array(self):
        return self._array

    @array.setter
    def array(self, array_value):
        raise ValueError('array property cannot be changed because it would loose its'
                         'dependency to geotransform information.')

    @property
    def geot(self):
        return self._geot

    @geot.setter
    def geot(self, geot_value):
        raise ValueError('geot property cannot be changed because it would loose its'
                         'dependency to array property.')

    def extract(self, lat_range, lon_range):
        x_range = [mapCoords_to_pixelCoords(geot=self._geot, xgeo=x) for x in lon_range]
        y_range = [mapCoords_to_pixelCoords(geot=self._geot, ygeo=y) for y in lat_range]

        array_extracted = self._array[y_range[1]:y_range[0], x_range[0]:x_range[1]]
        geot_extracted = build_geot(originX=lon_range[0] - self._geot[1] / 2.,
                                    originY=lat_range[1] - self._geot[5] / 2.,
                                    pixelWidth=self._geot[1],
                                    pixelHeight=self._geot[5])

        return Raster(array_extracted, geot_extracted)

    def save_as_tiff(self, filename):
        array2raster(filename, self._geot, self._array)


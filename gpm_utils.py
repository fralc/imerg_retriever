# -*- coding: utf-8 -*-
"""Functions to retrieve GPM files from PPS repositories
"""

from ftplib import FTP
import h5py
import gdal
import numpy as np
from datetime import datetime
import logging
import matplotlib.pyplot as plt
from gdal_utils import Raster, build_geot

# TODO read credentials from yaml file. Evalute even other parameters
_USERNAME = 'your_email'
_PASSWORD = 'your_email'

_FIELD_KEYS = ['dataType',
               'satellite',
               'instrument',
               'algorithmName',
               'startDate-SstartTime-EendTime',
               'sequenceIndicator',
               'VdataVersion',
               'extension']

_DATA_LEVELS = {'1A': 'Instrument count, geolocated, at instantaneous field of view (IFOV)',
                '1B': 'Geolocated, calibrated T b or radar power at IFOV',
                '1C': 'Intercalibrated brightness temperatures Tc at IFOV',
                '2A': 'Geolocated geophysical parameters at IFOV from a single instrument',
                '2B': 'Geolocated geophysical parameters at IFOV from multiple instruments',
                '3A': 'Space/time averaged geophysical parameters from a single instrument',
                '3B': 'Space/time averaged geophysical parameters from multiple instruments',
                '4': 'Combined satellite, ground and/or model data'}

_ACCUMULATIONS = {'HR': 'The product accumulates data for 1 hour',
                  'HHR': 'The product accumulates data every half hour',
                  'DAY': 'The product accumulates data for a single day',
                  'PENT': 'The product accumulates data for a 5-day period',
                  '7DAY': 'The product accumulates data for a 7-day period',
                  'MO': 'The product accumulates data for a designated month',
                  'ORBIT': 'The product accumulates data for each orbit'}

_REPOSITORIES = {'final': 'arthurhou.pps.eosdis.nasa.gov',
                 'early': 'jsimpson.pps.eosdis.nasa.gov',
                 'late': 'jsimpson.pps.eosdis.nasa.gov'}

# 3IMERGHH data fields, variable names and data units (IMERG_DOC, pag 19)
_IMERGHH_FIELDS = [
    ('name', 'Units', 'Data field Variable'),
    ('precipitationCal', 'mm/hr', 'Multi-satellite precipitation estimate with gauge calibration (recommended for '
                                  'general use)'),
    ('precipitationUncal', 'mm/hr', 'Multi-satellite precipitation estimate'),
    ('randomError', 'mm/hr', 'Random error for gauge-calibrated multi-satellite precipitation'),
    ('HQprecipitation', 'mm/hr', 'Merged microwave-only precipitation estimate'),
    ('HQprecipSource', 'index values', 'Microwave satellite source identifier'),
    ('HQobservationTime', 'min. into half hour', 'Microwave satellite observation time'),
    ('IRprecipitation', 'mm/hr', 'IR-only precipitation estimate'),
    ('IRkalmanFilterWeight', 'percent', 'Weighting of IR-only precipitation relative to the morphed merged '
                                        'microwave-only precipitation'),
    ('probabilityLiquidPrecipitation', 'percent', 'Probability of liquid precipitation phase'),
    ('precipitationQualityIndex', 'n/a', 'Quality Index for precipitationCal field')
]

_HQprecipSource_values = {0: 'no observation',1: 'TMI',
                          2: '(unused)',
                          3: 'AMSR',
                          4: 'SSMI',
                          5: 'SSMIS',
                          6: 'AMSU',
                          7: 'MHS',
                          8: 'SAPHIR',
                          9: 'GMI',
                          10: '(unused)',
                          11: 'ATMS',
                          12: 'AIRS',
                          13: 'TOVS',
                          14: 'CRIS',
                          15: 'Spare CSI 1',
                          16: 'Spare CSI 2',
                          17: 'Spare CSI 3',
                          18: 'Spare CSI 4',
                          19: 'Spare CSI 5',
                          20: 'Spare CTSS 1',
                          21: 'Spare CTSS 2',
                          22: 'Spare CTSS 3',
                          23: 'Spare CTSS 4',
                          24: 'Spare CTSS 5'}

_GEO_TRANSFORMATION_PARAMS = {'originX': -180,
                              'originY': 90,
                              'pixelWidth': 0.10000000149011612,
                              'pixelHeight': -0.10000000149011612}


def get_imerg(start_time, latency='final', ftp=None):
    """
    Download an IMERG file.
    Args:
        start_time (datetime.datetime): the starting date of the observation
        latency (str): can be 'early', 'late', 'final'
        ftp (ftplib.FTP): the ftp connection

    Returns:
        str: the downloaded filename
    """

    latency = latency.lower()

    # TODO latency='early' and latency='late' are not implemented
    if latency in ['early', 'late']:
        raise NotImplementedError('EARLY and LATE latencies still are not implemented.')
    elif latency != 'final':
        raise ValueError('latency must be "EARLY", "LATE", or "FINAL".')

    folder = '/gpmdata/{}/{}/{}/imerg/'.format(start_time.strftime('%Y'),
                                               start_time.strftime('%m'),
                                               start_time.strftime('%d'))

    # Get connection
    if ftp is None:
        ftp = _ftp_connect(_USERNAME, _PASSWORD, _REPOSITORIES[latency], folder)
        close_ftp = True
    else:
        close_ftp = False

    # Check folder
    if ftp.pwd() != folder:
        ftp.cwd(folder)

    # Get the filename with the smaller timedelta w.r.t. start_time
    content = ftp.nlst()
    filename = _get_closer_filename(content, start_time)

    # Retrieve the file
    _ftp_download(filename, ftp)
    logging.info('file {} has been downloaded.'.format(filename))

    if close_ftp:
        ftp.quit()

    return filename


def filename_info(filename, deep=False):
    """
    Provides information derived from a GPM file name.
    Args:
        filename (str): GPM data filename.
        deep (bool): if True more detailed information are derived (splitting sub-fields), if False only
            main information fields will be returned.

    Returns:
        dict: field values indexed by field names.

    """

    field_names = filename.split('.')

    if len(field_names) != 8:
        raise ValueError('The provided filename has not 8 fields')
    field_dict = {k: v for k, v in zip(_FIELD_KEYS, field_names)}

    if not deep:
        return field_dict

    infos = {}

    dataTypes = field_dict['dataType'].split('-')
    infos['level'] = dataTypes[0]
    if len(dataTypes) > 1:
        infos['accumulation'] = dataTypes[1]
        if len(dataTypes) > 2:
            infos['latency'] = dataTypes[2]

    infos['satellite'] = field_dict['satellite']
    infos['instrument'] = field_dict['instrument']
    infos['algorithmName'] = field_dict['algorithmName']

    times = field_dict['startDate-SstartTime-EendTime'].split('-')
    infos['timeUtcStart'] = datetime.strptime(times[0] + times[1], '%Y%m%dS%H%M%S')
    infos['timeUtcEnd'] = datetime.strptime(times[0] + times[2], '%Y%m%dE%H%M%S')

    infos['sequenceIndicator'] = field_dict['sequenceIndicator']
    infos['VdataVersion'] = field_dict['VdataVersion']
    infos['extension'] = field_dict['extension']

    return infos


def read_hdf5(hdf5_filename):
    """
    Open a IMERG HDF5 file and returns gridded datasets.
    Args:
        hdf5_filename (str): IMERG HDF5 filename.

    Returns:
        dict: dataset variables indexed with the name available in the HDF5 file.

    """
    # TODO verify if only IMERG have data into GRID structure
    f = h5py.File(hdf5_filename, 'r')
    datasets = {}
    for variable in f['Grid'].keys():
        datasets[variable] = f['Grid/' + variable].value
    f.close()
    return datasets


def datasets_to_rasters(datasets, geot):
    """
    Transform a selection of datasets (precipitationCal, precipitationUncal, probabilityLiquidPrecipitation,
        IRprecipitation) into rasters.
    Args:
        datasets (dict): a dictionary containing array datasets, e.g., as that obtained by read_hdf5.
        geot (tuple): GDAL geotransform tuple.

    Returns:
        dict: Raster instances indexed

    """
    # TODO verify that these are the only grid of interest
    precip_vars = [
     'precipitationCal',
     'precipitationUncal',
     'probabilityLiquidPrecipitation',
     'IRprecipitation'
    ]

    rasters = {}
    for k, v in datasets.items():
        if k in precip_vars:
            rasters[k] = Raster(np.flipud(v.T), geot)

    return rasters


def hdf5_to_rasters(hdf5_filename, geot):
    datasets = read_hdf5(hdf5_filename)
    rasters = datasets_to_rasters(datasets, geot)
    return rasters


def extract(datasets, grid_name, lat_range, lon_range):
    """
    Extract a portion of grid coming from a GPM dataset using 'lat' and 'lon' dataset variables
    Args:
        datasets (dict): gpm data as obtained from read_hdf5
        grid_name (str): name of the grid of interest
        lat_range (tuple): min and max latitude of the area to extract
        lon_range (tuple): min and max longitude of the area to extract

    Returns:
        numpy.ndarray: extracted grid

    """
    lat_mask = [all([i, s]) for i, s in zip(datasets['lat'] > lat_range[0], datasets['lat'] < lat_range[1])]
    lon_mask = [all([i, s]) for i, s in zip(datasets['lon'] > lon_range[0], datasets['lon'] < lon_range[1])]

    lat_pos = [e for e, x in enumerate(lat_mask) if x]
    lon_pos = [e for e, x in enumerate(lon_mask) if x]

    grid_cut = datasets[grid_name][min(lon_pos):max(lon_pos), min(lat_pos):max(lat_pos)]

    return grid_cut


def get_start_time(filename):
    times = filename.split('.')[4].split('-')
    return datetime.strptime(times[0] + times[1], '%Y%m%dS%H%M%S')


def _ftp_connect(username, password, repository, folder=None):
    ftp = FTP(repository, username, password)
    if folder is not None:
        ftp.cwd(folder)
    return ftp


def _ftp_download(filename, ftp):
    with open(filename, "wb") as f:
        ftp.retrbinary('RETR {}'.format(filename), f.write)
    logging.info('{} downloaded as {}'.format(filename, filename))


def _get_closer_filename(content, start_time):
    deltas = []
    for file in content:
        if file.split('.')[-1] == 'HDF5':
            deltas.append(abs(get_start_time(file) - start_time))
        else:
            deltas.append(np.nan)
    filename = content[np.argmin(deltas)]
    return filename







# -*- coding: utf-8 -*-
"""Functions to retrieve GPM files from PPS repositories
"""

from yaml import load
import os
from ftplib import FTP
import urllib.request
import h5py
import gdal
import numpy as np
from datetime import datetime, timedelta
import logging
import matplotlib.pyplot as plt
from gdal_utils import Raster, GribData, build_geot

with open('credentials.yaml', "r") as f:
    _credentials = load(f)
_USERNAME = _credentials['username']
_PASSWORD = _credentials['password']

with open('imerg_info.yaml', "r") as f:
    imerg_info = load(f)
_REPOSITORIES = imerg_info['REPOSITORIES']
_FIELD_KEYS = imerg_info['FIELD_KEYS']


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

    if latency not in ['early', 'late', 'final']:
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


def get_imerg_web(field_key_dict, start_time):
    # https://storm.pps.eosdis.nasa.gov/storm/NRT.jsp
    filename = build_filename(field_key_dict, start_time)
    url = "{}/NRT?email={}&filename=data/imerg/late/{}/{}".format(
        _REPOSITORIES['web_nrt'],
        _USERNAME,
        datetime.strftime(start_time, '%Y%m'),
        filename
    )
    print(url)
    urllib.request.urlretrieve(url, os.path.join('data', filename))
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


def build_filename(field_key_dict, start_time):
    field_key_dict['startDate-SstartTime-EendTime'] = "{}-S{}-E{}".format(
        datetime.strftime(start_time, '%Y%m%d'),
        datetime.strftime(start_time, '%H%M%S'),
        datetime.strftime(start_time + timedelta(minutes=29, seconds=59), '%H%M%S'),
    )
    field_key_dict['sequenceIndicator'] = str(int(
        (int(datetime.strftime(start_time, '%H')) + int(datetime.strftime(start_time, '%M')) / 60.) / .5 * 30
    ))
    filename = ".".join([field_key_dict[k] for k in _FIELD_KEYS])
    return filename


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


class H5Data(GribData):
    def __init__(self, filename):
        self._geot = build_geot(**imerg_info['GEO_TRANSFORMATION_PARAMS'])
        self._raw_data, self._raw_data_dict = H5Data.read_h5(filename, self._geot)
        for k, v in self._raw_data_dict.items():
            setattr(self, k, v)

    @staticmethod
    def read_h5(fn, geot):
        raster_dict = hdf5_to_rasters(fn, geot)
        infos = filename_info(fn, deep=True)
        layers = []
        layers_dict = {}
        for k, v in raster_dict.items():
            infos['Variable'] = k
            layers.append(
                (v,
                 infos['timeUtcStart'],
                 infos,
                 os.path.split(fn)[-1])
            )
            layers_dict[k] = layers[-1]
        return layers, layers_dict






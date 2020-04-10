#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Project specific function selecting observed data with svd decomposition.

:copyright:
    Amandine Sergeant (sergeant@lma.cnrs-mrs.fr), 2019
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""
import numpy as np
import obspy
import instaseis
from obspy.core.util import AttribDict
from obspy.io.xseed import Parser
from scipy import signal
import warnings

from lasif import LASIFError


def instaseis_synthetics_function(processing_info, iteration, db):  # NOQA
    """
    Function to compute the synthetics seismogram for one individual seismogram.
    It uses the same preprocessing parameters as for data seismograms
    defined in the iteration file.
    This is part of the project so it can change depending on the project.

    Please keep in mind that you will have to manually update this file to a
    new version if LASIF is ever updated.

    You can do whatever you want in this function as long as the function
    signature is honored. The file is read from ``"input_filename"`` and
    written to ``"output_filename"``.

    One goal of this function is to make sure that the data is available at the
    same time steps as the synthetics. The first time sample of the synthetics
    will always be the origin time of the event.

    Furthermore the data has to be converted to m/s.

    :param processing_info: A dictionary containing information about the
        file to be processed. It will have the following structure.
    :type processing_info: dict

    .. code-block:: python

        {'event_information': {
            'depth_in_km': 22.0,
            'event_name': 'GCMT_event_VANCOUVER_ISLAND...',
            'filename': '/.../GCMT_event_VANCOUVER_ISLAND....xml',
            'latitude': 49.53,
            'longitude': -126.89,
            'm_pp': 2.22e+18,
            'm_rp': -2.78e+18,
            'm_rr': -6.15e+17,
            'm_rt': 1.98e+17,
            'm_tp': 5.14e+18,
            'm_tt': -1.61e+18,
            'magnitude': 6.5,
            'magnitude_type': 'Mwc',
            'origin_time': UTCDateTime(2011, 9, 9, 19, 41, 34, 200000),
            'region': u'VANCOUVER ISLAND, CANADA REGION'},
         'input_filename': u'/.../raw/7D.FN01A..HHZ.mseed',
         'output_filename': u'/.../processed_.../7D.FN01A..HHZ.mseed',
         'process_params': {
            'dt': 0.75,
            'highpass': 0.007142857142857143,
            'lowpass': 0.0125,
            'npts': 2000},
         'station_coordinates': {
            'elevation_in_m': -54.0,
            'latitude': 46.882,
            'local_depth_in_m': None,
            'longitude': -124.3337},
         'station_filename': u'/.../STATIONS/RESP/RESP.7D.FN01A..HH*'}

    Please note that you also got the iteration object here, so if you
    want some parameters to change depending on the iteration, just use
    if/else on the iteration objects.

    >>> iteration.name  # doctest: +SKIP
    '11'
    >>> iteration.get_process_params()  # doctest: +SKIP
    {'dt': 0.75,
     'highpass': 0.01,
     'lowpass': 0.02,
     'npts': 500}

    Use ``$ lasif shell`` to play around and figure out what the iteration
    objects can do.

    """

    # =========================================================================
    # Step 1: Compute synthetic seismogram
    # =========================================================================
    st = db.get_seismograms(
        source=processing_info["event_information"]["source"],
        receiver=processing_info["receiver"],
        kind='velocity',
        components=processing_info["component"],
        remove_source_shift=True,
        dt=processing_info["process_params"]["dt"])

    # =========================================================================
    # Step 2: Preprocess the synthetic seismogram
    # =========================================================================
    # Use same preprocessing as for data seismograms
    # see the project preprocessing_function.py
    freqmin = processing_info["process_params"]["highpass"]
    freqmax = processing_info["process_params"]["lowpass"]

    st.detrend("linear")
    st.detrend("demean")
    st.taper(0.05, type="cosine")
    st.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=3,
              zerophase=False)
    st.detrend("linear")
    st.detrend("demean")
    st.taper(0.05, type="cosine")
    st.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=3,
              zerophase=False)

    starttime = processing_info["event_information"]["origin_time"]
    endtime = starttime + processing_info["process_params"]["dt"] * \
        (processing_info["process_params"]["npts"] - 1)
    st = st.trim(starttime, endtime)

    tr = st[0]
    tr.stats.coordinates = AttribDict({
        'latitude': processing_info["station_coordinates"]["latitude"],
        'elevation': processing_info["station_coordinates"]["elevation_in_m"],
        'longitude': processing_info["station_coordinates"]["longitude"]})
    tr.stats.station = processing_info["station"]

    # =========================================================================
    # Save processed synthetic and clean up.
    # =========================================================================
    # Convert to single precision to save some space.
    tr.data = np.require(tr.data, dtype="float32", requirements="C")
    tr.stats._format = processing_info["output_filename"].split('.')[-1]
    if hasattr(tr.stats, "mseed"):  # to be fixed
        tr.stats.mseed.encoding = "FLOAT32"

    tr.write(processing_info["output_filename"], format=tr.stats._format)

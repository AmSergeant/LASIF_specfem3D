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
import os
import numpy as np
import obspy
import warnings

from lasif import LASIFError


def stf_deconvolution(to_be_processed, output_folder, components=['E', 'N', 'Z'],):  # NOQA
    """
    Function to perform the actual preprocessing for one individual seismogram.
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

    def source_deconvolution_freq(
            stream_data,
            stream_green,
            lambd=0.001,
            recompute_syn=False):

        # Calculate STF:
        # deconvoluate the Green's functions from observed seismograms
        # following Pratt 1999, equation 17
        nfft = stream_data[0].stats.npts
        num = np.zeros(nfft, dtype=complex)
        den = np.zeros(nfft, dtype=complex)
        chi_obs = []
        for tr, sy in zip(stream_data, stream_green):
            tr_fft = np.fft.fft(tr.data, nfft)
            sy_fft = np.fft.fft(sy.data, nfft)

            num += np.conjugate(sy_fft) * tr_fft
            den += np.conjugate(sy_fft) * sy_fft
            chi_obs.append(np.sum(tr.data**2))
        chi_obs = 0.5 * np.sum(chi_obs)

        water_level = lambd * np.max(np.abs(den))
        s = num / (den + water_level)
        src = np.real(np.fft.ifft(s))
        stream_src = obspy.Stream()
        stream_src += tr.copy()
        stream_src[0].stats.station = ''
        stream_src[0].data = src
        residual = []
        stream_syn = obspy.Stream()

        # recompute synthetics with the estimated STF
        if recompute_syn:
            src_fft = np.fft.fft(src, nfft)
            chi_syn = []
            for tr, sy in zip(stream_data, stream_green):
                sy_fft = np.fft.fft(sy.data, nfft)
                cal = sy.copy()
                cal.data = np.real(np.fft.ifft(src_fft * sy_fft))
                stream_syn += cal
                res = tr.data - cal.data
                chi_syn.append(np.sum(res**2))
            chi_syn = 0.5 * np.sum(chi_syn)
            residual = chi_syn / chi_obs

        return stream_src, stream_syn, residual

    # =========================================================================
    # Entering the function
    # =========================================================================
    from matplotlib.dates import date2num, num2date
    SECONDS_PER_DAY = 3600 * 24

    process_params = to_be_processed[0]["processing_info"]["process_params"]
    seconds_prior_arrival = process_params["seconds_prior_arrival"]
    window_length_in_sec = process_params["window_length_in_sec"]

    for comp in components:
        # =========================================================================
        # Component selection
        # =========================================================================
        # !!!!!! replace first_P_arrival by phase of interest to be calculated in preprocess_data and given in process_info

        wav_file_list = [wav["processing_info"]["input_filename"]
                         for wav in to_be_processed
                         if comp in wav["processing_info"]["channel"]]
        syn_file_list = [wav["processing_info"]["output_filename"]
                         for wav in to_be_processed
                         if comp in wav["processing_info"]["channel"]]
        first_arrival = [wav["processing_info"]["first_P_arrival"]
                         for wav in to_be_processed
                         if comp in wav["processing_info"]["channel"]]

        idx_sigwin_start = int(
            np.ceil(
                (np.min(first_arrival) -
                 process_params["seconds_prior_arrival"]) /
                process_params["dt"]))
        idx_sigwin_end = int(
            np.ceil(
                (np.max(first_arrival) +
                 process_params["window_length_in_sec"]) /
                process_params["dt"]))

        Time = np.arange(
            0,
            process_params["npts"] *
            process_params["dt"],
            process_params["dt"])
        starttime = to_be_processed[0]["processing_info"]["event_information"]["origin_time"]
        t_start = num2date(((Time[idx_sigwin_start] / SECONDS_PER_DAY)
                            + date2num(starttime.datetime)))
        t_end = num2date(((Time[idx_sigwin_end] / SECONDS_PER_DAY)
                          + date2num(starttime.datetime)))
        startdate = obspy.UTCDateTime(
            t_start.year,
            t_start.month,
            t_start.day,
            t_start.hour,
            t_start.minute,
            t_start.second,
            t_start.microsecond)
        enddate = obspy.UTCDateTime(t_end.year, t_end.month, t_end.day,
                                    t_end.hour, t_end.minute, t_end.second,
                                    t_end.microsecond)

        # =========================================================================
        # read traces, window around phase of interest
        # =========================================================================
        st_wav = obspy.Stream()
        st_syn = obspy.Stream()
        for wav_file, syn_file in zip(wav_file_list, syn_file_list):
            wav = obspy.read(wav_file)
            syn = obspy.read(syn_file)
            wav.trim(startdate, enddate)
            syn.trim(startdate, enddate)
            #wav[0].data = wav[0].data[idx_sigwin_start:idx_sigwin_end]
            #syn[0].data = syn[0].data[idx_sigwin_start:idx_sigwin_end]
            wav[0].data /= np.max(wav[0].data)
            syn[0].data /= np.max(syn[0].data)
            st_wav += wav
            st_syn += syn

        # if no waveform selected at the previous step (snr criteria), quit the
        # process
        if not st_wav or not st_syn:
            raise LASIFError(
                "No data for this event, will skip the stf estimation")
        else:

            st_wav.taper(0.01)
            st_syn.taper(0.01)

            # =========================================================================
            # stf deconvolution
            # =========================================================================
            # stf, new_syn, residual = source_deconvolution_freq(st_wav, st_syn,
            # lambd=0.001, recompute_syn=True)
            stf, p, pp = source_deconvolution_freq(
                st_wav, st_syn, lambd=0.001, recompute_syn=False)

            '''
            src = obspy.read(wav_file_list[0])
            src[0].stats.station=''
            src_trace = np.zeros(process_params["npts"], dtype=float)
            src_trace[idx_sigwin_start:idx_sigwin_end] = stf[0].data
            src[0].data = src_trace.copy()
            stf = src.copy()
            '''

            # =========================================================================
            # write stf file
            # =========================================================================
            # Convert to single precision to save some space.
            tr = stf[0].copy()
            tr.data = np.require(tr.data, dtype="float32", requirements="C")
            tr.stats._format = wav[0].stats._format
            if hasattr(tr.stats, "mseed"):  # to be fixed
                tr.stats.mseed.encoding = "FLOAT32"

            # channel_id = [item["processing_info"]["channel"]
            #         for item in to_be_processed
            #         if comp in item["processing_info"]["channel"]][0]
            stf_filename = os.path.join(output_folder, "stf_%s__%s__%s"
                                        % (comp,
                                           to_be_processed[0]["processing_info"]["output_filename"].split('/')[-1].split('__')[-2],
                                           to_be_processed[0]["processing_info"]["output_filename"].split('/')[-1].split('__')[-1]))
            tr.write(stf_filename, format=tr.stats._format)

            """
            # plotting
            import matplotlib.pyplot as plt
            from lasif.visualization import plot_waveform_section
            st_green = st_syn.copy()
            '''
            st_wav = obspy.Stream()
            st_green = obspy.Stream()
            for wav_file, syn_file in zip(wav_file_list, syn_file_list):
                wav = obspy.read(wav_file)
                syn = obspy.read(syn_file)
                wav[0].data /= np.max(wav[0].data)
                syn[0].data /= np.max(syn[0].data)
                st_wav += wav
                st_green += syn
            '''

            ax1 = plt.subplot(121)
            plot_waveform_section(ax1,st_wav, [], scale=.5, colors='k')
            plot_waveform_section(ax1,st_green, [], scale=.5, colors='b')

            ax2 = plt.subplot(122)
            nfft = st_wav[0].stats.npts
            st_syn = obspy.Stream()
            src_fft = np.fft.fft(stf[0].data,nfft)
            for sy in st_green:
                sy_fft = np.fft.fft(sy.data,nfft)
                cal = sy.copy()
                cal.data = np.real(np.fft.ifft(src_fft*sy_fft))
                st_syn += cal
            plot_waveform_section(ax2,st_wav, [], scale=.5, colors='k')
            plot_waveform_section(ax2,st_syn, [], scale=.5, colors='r')
            plt.show()
            """

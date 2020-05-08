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
import warnings

from lasif import LASIFError


def stf_deconvolution(to_be_processed, output_folder, components = ['E','N','Z'],):  # NOQA
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

    def source_deconvolution_freq(stream_data, stream_green, lambd=0.001, recompute_syn=False):

        # Calculate STF:
        # deconvoluate the Green's functions from observed seismograms 
        # following Pratt 1999, equation 17
        nfft = stream_data[0].stats.npts
        num = np.zeros(nfft,dtype = complex)
        den = np.zeros(nfft,dtype = complex)
        chi_obs = []
        for tr, sy in zip(stream_data, stream_green):
            tr_fft = np.fft.fft(tr.data,nfft)
            sy_fft = np.fft.fft(sy.data,nfft)

            num += np.conjugate(sy_fft)*tr_fft
            den += np.conjugate(sy_fft)*sy_fft
            chi_obs.append(np.sum(tr.data**2))
        chi_obs = 0.5*np.sum(chi_obs)


        water_level = lambd*np.max(np.abs(den))
        s = num/(den + water_level)
        src = np.real(np.fft.ifft(s))
        stream_src = obspy.Stream()
        stream_src += tr.copy()
        stream_src[0].stats.station=''
        stream_src[0].data = src
        residual = []
        stream_syn = obspy.Stream()

        # recompute synthetics with the estimated STF
        if recompute_syn:
            src_fft = np.fft.fft(src,nfft)
            chi_syn = []
            for tr, sy in zip(stream_data, stream_green): 
                sy_fft = np.fft.fft(sy.data,nfft)
                cal = sy.copy()
                cal.data = np.real(np.fft.ifft(src_fft*sy_fft))
                stream_syn += cal
                res = tr.data - cal.data
                chi_syn.append(np.sum(res**2))
            chi_syn = 0.5*np.sum(chi_syn)
            residual = chi_syn/chi_obs

        return stream_src, stream_syn, residual



    # =========================================================================
    # Entering the function
    # =========================================================================
    from matplotlib.dates import date2num, num2date
    import obspy
    SECONDS_PER_DAY = 3600*24

    process_params = to_be_processed[0]["processing_info"]["process_params"]

    #process_params["window_length_in_sec"] = 30.
    seconds_prior_arrival = process_params["seconds_prior_arrival"]
    window_length_in_sec = process_params["window_length_in_sec"]


    evla = to_be_processed[0]["processing_info"]["event_information"]["latitude"]
    evlo = to_be_processed[0]["processing_info"]["event_information"]["longitude"]
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

        idx_sigwin_start = int(np.ceil((np.min(first_arrival) -
                                        process_params["seconds_prior_arrival"]) / process_params["dt"]))
        idx_sigwin_end = int(np.ceil((np.max(first_arrival) +
                                      (process_params["window_length_in_sec"] - 
                                       process_params["seconds_prior_arrival"]))/ process_params["dt"]))


        Time = np.arange(0,process_params["npts"]*process_params["dt"], process_params["dt"])
        starttime = to_be_processed[0]["processing_info"]["event_information"]["origin_time"]
        t_start = num2date( ((Time[idx_sigwin_start] / SECONDS_PER_DAY) 
                             + date2num(starttime.datetime)) )
        t_end = num2date( ((Time[idx_sigwin_end] / SECONDS_PER_DAY) 
                           + date2num(starttime.datetime)) )
        startdate = obspy.UTCDateTime(t_start.year, t_start.month, t_start.day, 
                                      t_start.hour, t_start.minute, t_start.second,
                                      t_start.microsecond)
        enddate = obspy.UTCDateTime(t_end.year, t_end.month, t_end.day, 
                                    t_end.hour, t_end.minute, t_end.second,
                                    t_end.microsecond)


        # =========================================================================
        # read traces, window around phase of interest 
        # =========================================================================
        st_wav = obspy.Stream()
        st_syn = obspy.Stream()
        filenames = []
        back_azimuths = []
        for wav_file, syn_file in zip(wav_file_list, syn_file_list):
            if os.path.isfile(wav_file) and os.path.isfile(syn_file):
                wav = obspy.read(wav_file)
                syn = obspy.read(syn_file)
            else:
                continue


            if comp == 'R':
                # need to load T data and synthetics for RT-> EN rotation
                channel_id = wav[0].stats.network+'.'+wav[0].stats.station+'.'+wav[0].stats.location+'.'+wav[0].stats.channel
                T_channel_id = channel_id.replace(wav[0].stats.channel, wav[0].stats.channel[:-1]+'T')
                T_wav_file = wav_file.replace(channel_id, T_channel_id)
                T_syn_file = syn_file.replace(channel_id, T_channel_id)
                #print("for R: %s, reading T component %s"%(wav_file, T_wav_file))
                if os.path.isfile(T_wav_file) and os.path.isfile(T_syn_file):
                    wav += obspy.read(T_wav_file)          
                    syn += obspy.read(T_syn_file)

                    stla = [proc["processing_info"]["station_coordinates"]["latitude"]
                            for proc in to_be_processed 
                            if proc["processing_info"]["input_filename"] == wav_file][0]
                    stlo = [proc["processing_info"]["station_coordinates"]["longitude"]
                            for proc in to_be_processed 
                            if proc["processing_info"]["input_filename"] == wav_file][0]
                    from obspy.geodetics.base import gps2dist_azimuth
                    epicentral_distance, azimuth, baz = gps2dist_azimuth(evla, evlo, stla, stlo)
                    wav.rotate('RT->NE', back_azimuth = baz)
                    syn.rotate('RT->NE', back_azimuth = baz)
                    back_azimuths.append(baz)

                    wav[0].data /= np.max(wav[0].data)
                    syn[0].data /= np.max(syn[0].data)
                    #wav[0].data = wav[0].data[idx_sigwin_start:idx_sigwin_end]
                    #syn[0].data = syn[0].data[idx_sigwin_start:idx_sigwin_end]
                    wav.trim(startdate,enddate)
                    syn.trim(startdate,enddate)

                    st_wav += wav
                    st_syn += syn
                    filenames.append(os.path.basename(wav_file))
            elif comp == 'T':
                # need to load R data and synthetics for RT-> EN rotation
                channel_id = wav[0].stats.network+'.'+wav[0].stats.station+'.'+wav[0].stats.location+'.'+wav[0].stats.channel
                R_channel_id = channel_id.replace(wav[0].stats.channel, wav[0].stats.channel[:-1]+'R')
                R_wav_file = wav_file.replace(channel_id, R_channel_id)
                R_syn_file = syn_file.replace(channel_id, R_channel_id)
                #print("for T: %s, reading R component %s"%(wav_file, R_wav_file))
                if os.path.isfile(R_wav_file) and os.path.isfile(T_syn_file):
                    wav += obspy.read(R_wav_file)          
                    syn += obspy.read(R_syn_file)

                    stla = [proc["processing_info"]["station_coordinates"]["latitude"]
                            for proc in to_be_processed 
                            if proc["processing_info"]["input_filename"] == wav_file][0]
                    stlo = [proc["processing_info"]["station_coordinates"]["longitude"]
                            for proc in to_be_processed 
                            if proc["processing_info"]["input_filename"] == wav_file][0]
                    from obspy.geodetics.base import gps2dist_azimuth
                    epicentral_distance, azimuth, baz = gps2dist_azimuth(evla, evlo, stla, stlo)
                    wav.rotate('RT->NE', back_azimuth = baz)
                    syn.rotate('RT->NE', back_azimuth = baz)
                    back_azimuths.append(baz)

                    wav[0].data /= np.max(wav[0].data)
                    syn[0].data /= np.max(syn[0].data)
                    #wav[0].data = wav[0].data[idx_sigwin_start:idx_sigwin_end]
                    #syn[0].data = syn[0].data[idx_sigwin_start:idx_sigwin_end]
                    wav.trim(startdate,enddate)
                    syn.trim(startdate,enddate)

                    st_wav += wav
                    st_syn += syn
                    filenames.append(os.path.basename(wav_file))

            else:
                wav[0].data /= np.max(wav[0].data)
                syn[0].data /= np.max(syn[0].data)
                #wav[0].data = wav[0].data[idx_sigwin_start:idx_sigwin_end]
                #syn[0].data = syn[0].data[idx_sigwin_start:idx_sigwin_end]
                wav.trim(startdate,enddate)
                syn.trim(startdate,enddate)

                st_wav += wav
                st_syn += syn
                filenames.append(os.path.basename(wav_file))

        # if no waveform selected at the previous step (snr criteria), quit the process
        if not st_wav or not st_syn:
            raise LASIFError("No data for this event, will skip the stf estimation")
        else:

            st_wav.taper(0.01)
            st_syn.taper(0.01)

            '''
            from lasif.visualization import plot_waveform_section
            import matplotlib.pyplot as plt
            fig0 = plt.figure()  
            ax1 = fig0.add_subplot(121)
            plot_waveform_section(ax1,st_wav.select(component='E'), [], scale=1, colors='k')
            plot_waveform_section(ax1,st_syn.select(component='E'), [], scale=1, colors='r')
            ax2 = fig0.add_subplot(122)
            plot_waveform_section(ax2,st_wav.select(component='N'), [], scale=1, colors='k')
            plot_waveform_section(ax2,st_syn.select(component='N'), [], scale=1, colors='r')
            plt.show()  
            '''


            # =========================================================================
            # stf deconvolution
            # =========================================================================
            #stf, new_syn, residual = source_deconvolution_freq(st_wav, st_syn,
            #                                               lambd=0.001, recompute_syn=True)
            if comp == 'Z' or comp == 'E' or comp == 'N':
                stf, st_syn_from_stf, residual = source_deconvolution_freq(st_wav.select(component = comp), 
                                                                           st_syn.select(component = comp),
                                                                           lambd=0.001, recompute_syn=True)


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
                if hasattr(tr.stats, "mseed"): # to be fixed
                    tr.stats.mseed.encoding = "FLOAT32"

                #channel_id = [item["processing_info"]["channel"] 
                #         for item in to_be_processed 
                #         if comp in item["processing_info"]["channel"]][0]


                stf_filename = os.path.join(output_folder, "stf_%s__%s__%s"\
                                            %(comp,
                                              to_be_processed[0]["processing_info"]["output_filename"].split('/')[-1].split('__')[-2],
                                              to_be_processed[0]["processing_info"]["output_filename"].split('/')[-1].split('__')[-1]))
                tr.write(stf_filename, format=tr.stats._format)
                print(stf_filename)

                # =========================================================================
                # write synthetic from stf files
                # =========================================================================
                for tr, filename in zip(st_syn_from_stf, filenames):
                    if hasattr(tr.stats, "mseed"): # to be fixed
                        tr.stats.mseed.encoding = "FLOAT32"

                    tr.write(os.path.join(output_folder, filename), format=tr.stats._format)


            # for R and T components, comute both East and North stf        
            elif comp == 'R' or comp == 'T':

                stf_E, st_syn_from_stf_E, residual = source_deconvolution_freq(st_wav.select(component = 'E'), 
                                                                               st_syn.select(component = 'E'),
                                                                               lambd=0.001, recompute_syn=True)
                stf_N, st_syn_from_stf_N, residual = source_deconvolution_freq(st_wav.select(component = 'N'), 
                                                                               st_syn.select(component = 'N'),
                                                                               lambd=0.001, recompute_syn=True)
                tr_E = stf_E[0].copy()
                tr_E.data = np.require(tr_E.data, dtype="float32", requirements="C")
                tr_E.stats._format = wav[0].stats._format
                tr_N = stf_N[0].copy()
                tr_N.data = np.require(tr_N.data, dtype="float32", requirements="C")
                tr_N.stats._format = wav[0].stats._format
                if hasattr(tr_N.stats, "mseed"): # to be fixed
                    tr_N.stats.mseed.encoding = "FLOAT32"
                    tr_E.stats.mseed.encoding = "FLOAT32"

                #channel_id = [item["processing_info"]["channel"] 
                #         for item in to_be_processed 
                #         if comp in item["processing_info"]["channel"]][0]


                E_stf_filename = os.path.join(output_folder, "stf_E__%s__%s"\
                                              %(to_be_processed[0]["processing_info"]["output_filename"].split('/')[-1].split('__')[-2],
                                                to_be_processed[0]["processing_info"]["output_filename"].split('/')[-1].split('__')[-1]))
                tr_E.write(E_stf_filename, format=tr_E.stats._format)
                print(E_stf_filename)
                N_stf_filename = os.path.join(output_folder, "stf_N__%s__%s"\
                                              %(to_be_processed[0]["processing_info"]["output_filename"].split('/')[-1].split('__')[-2],
                                                to_be_processed[0]["processing_info"]["output_filename"].split('/')[-1].split('__')[-1]))
                tr_N.write(N_stf_filename, format=tr_N.stats._format)
                print(N_stf_filename)

                #ST= obspy.Stream()      
                for tr_syn_from_stf_E, tr_syn_from_stf_N, baz, filename in zip(st_syn_from_stf_E, st_syn_from_stf_N, back_azimuths, filenames):
                    stream = obspy.Stream()
                    stream += tr_syn_from_stf_E 
                    stream += tr_syn_from_stf_N
                    #ST += tr_syn_from_stf_E.copy() #stream.select(component = 'E')
                    #ST += tr_syn_from_stf_N.copy() #stream.select(component = 'N')
                    stream = stream.rotate('NE->RT', back_azimuth = baz)


                    if comp == 'R':
                        tr = stream.select(component = 'R').copy()
                    elif comp == 'T':
                        tr = stream.select(component = 'T').copy()
                    tr = tr[0]
                    if hasattr(tr.stats, "mseed"): # to be fixed
                        tr.stats.mseed.encoding = "FLOAT32"
                    tr.write(os.path.join(output_folder, filename), format=tr.stats._format)
                    #ST += stream
                '''
                    from lasif.visualization import plot_waveform_section
                    import matplotlib.pyplot as plt
                    fig0 = plt.figure()  
                    ax1 = fig0.add_subplot(121)
                    plot_waveform_section(ax1,st_wav.select(component='E'), [], scale=1, colors='k')
                    plot_waveform_section(ax1,st_syn.select(component='E'), [], scale=1, colors='r')
                    plot_waveform_section(ax1,ST.select(component='R'), [], scale=1, colors='b')
                    ax2 = fig0.add_subplot(122)
                    plot_waveform_section(ax2,st_wav.select(component='N'), [], scale=1, colors='k')
                    plot_waveform_section(ax2,st_syn.select(component='N'), [], scale=1, colors='r')
                    plot_waveform_section(ax2,ST.select(component='T'), [], scale=1, colors='g')
                    plt.show()
                    '''




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
            fig0 = plt.figure()  
            ax1 = fig0.add_subplot(121)
            plot_waveform_section(ax1,st_wav, [], scale=.5, colors='k')
            plot_waveform_section(ax1,st_green, [], scale=.5, colors='b')
            
            ax2 = fig0.add_subplot(122)
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
            #plt.show()
            """


            '''
            # window each trace in the phase window
            st_green = st_syn.copy()
            #st_syn = st_syn_from_stf.copy()
            origin_time = to_be_processed[0]["processing_info"]["event_information"]["origin_time"]
            # event origin time is the starttime of data seismograms
            
            #Time_trace_in_sec = np.arange(0,process_params["npts"]*process_params["dt"], process_params["dt"])
            #Time_trace_datetime = num2date( ((Time_trace_in_sec / SECONDS_PER_DAY) 
            #                 + date2num(origin_time.datetime)) )
            
            # here st_wav and st_syn_from_stf are already windowed
            first_tt_arrival = np.min(first_arrival)
            #first_tt_arrival = stf[0].stats.starttime - origin_time
            st_wav_win = obspy.Stream()
            st_syn_win = obspy.Stream()
            st_green_win = obspy.Stream()
            for tr_wav, tr_syn, tr_green, tt_arrival in zip(st_wav,st_syn_from_stf, st_green, first_arrival):
                #print(tt_arrival, first_tt_arrival)
                idx_sigwin_start = int(np.ceil((tt_arrival - first_tt_arrival) / process_params["dt"]))
                idx_sigwin_end = idx_sigwin_start + int(process_params["window_length_in_sec"]/ process_params["dt"]) - 1
                #print(tr_wav.stats.npts)
                # pb t_start can become negative !!!!!! weird
                
                t_start = num2date( ((tr_wav.times()[idx_sigwin_start] / SECONDS_PER_DAY) 
                             + date2num(tr_wav.stats.starttime.datetime)) )
                t_end = num2date( ((tr_wav.times()[idx_sigwin_end] / SECONDS_PER_DAY) 
                                   + date2num(tr_wav.stats.starttime.datetime)) )
                #print(idx_sigwin_start, idx_sigwin_end, t_start, t_end)
                startdate = obspy.UTCDateTime(t_start.year, t_start.month, t_start.day, 
                                              t_start.hour, t_start.minute, t_start.second,
                                              t_start.microsecond)
                enddate = obspy.UTCDateTime(t_end.year, t_end.month, t_end.day, 
                                              t_end.hour, t_end.minute, t_end.second,
                                              t_end.microsecond)
                st_wav_win += tr_wav.trim(startdate,enddate)
                st_syn_win += tr_syn.trim(startdate,enddate)
                st_green_win += tr_green.trim(startdate,enddate)
                
            
           
            
            # compute time-shifts anomaly with cross-correlation
            def align_trace(timeserie,tau):
                """
                Parameters
                ----------
                timeserie : np.array, dtype = float
                    seismic trace
                tau : scalar, dtype = int
                    number of samples to shift the trace
                    tau > 0 whill shift the trace forward in time
                    tau <0 will shift the trace backward in time
        
                Returns
                -------
                align : np.array, dtype = float
                    shifted trace. 
                    output trace has same number of samples as input signal which has been filled with zero padding
        
                """
                tau=int(tau)
                align=np.zeros(timeserie.shape[0])
                if tau>=0:
                    #roll the signal forward
                    align[tau:-1]=timeserie[0:-1-tau]
                elif tau<0:
                    #roll the signal backward
                    align[0:-1-np.abs(tau)]=timeserie[np.abs(tau):-1]
                return align
        
            time_shift = []
            st_syn_recal = obspy.Stream()
            # if time_shift > 0, data P wave is faster than synthetics
            # if time_shift < 0, data P wave is slower than synthetics
        
            from collections import OrderedDict
            stations = OrderedDict()
            
            from obspy.geodetics.base import locations2degrees
            import obspy.signal.freqattributes
            npts = st_wav_win[0].stats.npts
            nfft = obspy.signal.util._npts2nfft(npts)
            goal_freq = 1./10
            evla = to_be_processed[0]["processing_info"]["event_information"]["latitude"]
            evlo = to_be_processed[0]["processing_info"]["event_information"]["longitude"]
            offsets= []
            for tr_wav, tr_syn, tr_green, wav_file, tt_arrival in zip(st_wav_win, st_syn_win, st_green_win, wav_file_list, first_arrival):
                cc = np.correlate(tr_syn.data,tr_wav.data, mode="full")
                ts = (cc.argmax() - tr_wav.stats.npts  + 1) * process_params["dt"]
                time_shift.append(ts)
                recal = tr_syn.copy()
                recal.data = align_trace(tr_syn.data, -int((cc.argmax() - tr_wav.stats.npts  + 1)))
                st_syn_recal += recal
                
                a = np.dot(tr_syn.data, tr_wav.data)
                b = np.linalg.norm(tr_syn.data,ord=None)*np.linalg.norm(tr_wav.data,ord=None)
                cc = a/b
                
                a = np.dot(recal.data, tr_wav.data)
                b = np.linalg.norm(tr_syn.data,ord=None)*np.linalg.norm(tr_wav.data,ord=None)
                cc_shift = a/b
                
                A = np.dot(tr_wav.data, np.transpose(tr_green.data))\
                        /np.dot(tr_green.data, np.transpose(tr_green.data))
                
            
                wav_fft = np.fft.fft(tr_wav.data, nfft)
                wav_fft = wav_fft[0:nfft/2]
                Freq = np.linspace(0, 0.5*tr_wav.stats.sampling_rate, len(wav_fft))
                phase_angle = np.angle(wav_fft[np.argmin(np.abs(Freq-goal_freq))], deg=True)
                
                channel_id = '%s.%s.%s.%s'%(tr_wav.stats.network,
                                            tr_wav.stats.station,
                                            tr_wav.stats.location,
                                            tr_wav.stats.channel)
                station_name = tr_wav.stats.station
                
                station_info = [info["processing_info"]["station_coordinates"] 
                                for info in to_be_processed
                                if channel_id in info["processing_info"]["channel_id"]][0]
                offset = locations2degrees(evla, evlo, 
                station_info["latitude"],station_info["longitude"])
                offsets.append(offset)
                stations[station_name] = {
                "channel_id": channel_id,
                "latitude": float(station_info["latitude"]),
                "longitude": float(station_info["longitude"]),
                "time shift": float(ts),
                "input_file": wav_file,
                "epicentral_distance": offset,
                "tt_arrival": tt_arrival,
                "phase_angle": phase_angle,
                "cc_shift": cc_shift,
                "cc": cc,
                "amp anomaly": A}
              
            
            """
            # plotting
            import matplotlib.pyplot as plt
            from lasif.visualization import plot_waveform_section
            
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            plot_waveform_section(ax1,st_wav_win, offsets, scale=0.2, colors='k')
            plot_waveform_section(ax1,st_syn_win, offsets, scale=0.2, colors='b')
            plot_waveform_section(ax1,st_syn_recal, offsets, scale=0.2, colors='r')
            
            #import lasif.components.visualizations
            #lasif.components.visualizations.plot_stations(0)
            #plt.show()
            
            print(stations)
            """
            
            return stations
            '''



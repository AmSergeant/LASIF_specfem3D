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


def data_svd_selection(to_be_processed, components = ['E','N','Z'], cc_threshold = 0.7):  # NOQA
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

    def construct_reference_wavefrom_from_dtt_inversion(stream):
        """
        construct the reference waveform for teleseismic seismograms based on a svd decomposition of the traces
        see Bekara, M., and M. Van der Baan (2007), Local singular value decomposition for signal enhancement of seismic data
        for details on the methodology
        Will align traces based on the dtt inverted for every combination of two traces (see invert_time_shift)

        Parameters
        ----------
        stream : stream of data. Prefer to use traces windowed around the phase of interest


        Returns
        -------
        reference_waveform : np.array of same length as the individual traces in stream
        stream_aligned : stream
                        the input stream but with individual traces aligned with each other

        """


        def invert_time_shift(stream):
            """
            invert for the best time-shifts for m indivual traces to align them with respect to each other
            will calculate the time-shifts for every combination of 2 traces (N) in the stream object contained in "d"
            and invert for the least-square solution "m" following: d = Gm
            The invert of the non-square matrix G: Ginv = inv(G'G)*G'
            with d = (t1_2, t1_3, t1_4, ..., t1_m, t2_3, t2_4, ..., tm-1_m) of size N x 1
                 m = (t1, t2, ..., tm) of size m x 1
                 G = [-1 1 0 0 0 0 ... 0 0
                      -1 0 1 0 0 0 ... 0 0
                      -1 0 0 1 0 0 ... 0 0
                      -1 0 0 0 0 0 ... 0 1
                      0 -1 1 0 0 0 ... 0 0
                      0 -1 0 1 0 0 ... 0 0
                            ....
                      0 0 0 0 0 0 ... 1 -1] of size N x m
            Following Brenguier et al. (2014, Science) Supp. Mat


            Parameters
            ----------
            stream : data stream of length m
                contains seismic traces to cac

            Returns
            -------
            model_time_shift : np.array of size m x 1

            """
            import itertools
            m = len(stream)
            all_combin = list(itertools.combinations(np.arange(m), 2))
            N = len(all_combin)
            time_shift_doublet=np.zeros((N,1),dtype = float)
            G=np.zeros((N,m),dtype = int) 
            for i, combin in enumerate(all_combin):
                i1 = combin[0]
                i2 = combin[1]
                cc = np.correlate(stream[i1].data,stream[i2].data, mode="full")
                time_shift = int(cc.argmax() - stream[0].stats.npts  + 1)
                time_shift_doublet[i] = time_shift
                G[i][i1] = -1
                G[i][i2] = 1
            Ginv = np.linalg.inv(np.transpose(G).dot(G)).dot(np.transpose(G))
            model_time_shift = Ginv.dot(time_shift_doublet)
            return model_time_shift


        # Invert the time-shift for individual traces
        time_shifts = invert_time_shift(stream)

        # align traces + construct array of aligned data
        data_matrix = np.zeros((len(stream),stream[0].stats.npts),dtype = float)
        i = 0
        stream_aligned = stream.copy()
        for i, tr in enumerate(stream_aligned):
            tr.data = align_trace(tr.data, time_shifts[i] )
            data_matrix[i] = tr.data  
            i += 1 

        #rank = np.linalg.matrix_rank(data_matrix)
        #print("!!!! Rank %d should be lower than number of traces %d !!!"%(rank, len(stream)))

        # svd decomposition of data_matrix
        U, s, V = np.linalg.svd(data_matrix)

        # the reference waveform is the first eigenimage
        reference_wav = - np.sign(s[0])*V[0] # check with the sign !

        return reference_wav, stream_aligned


    def construct_reference_wavefrom_from_dtt_stack(stream, stack):
        """
        construct the reference waveform for teleseismic seismograms based on a svd decomposition of the traces
        see Bekara, M., and M. Van der Baan (2007), Local singular value decomposition for signal enhancement of seismic data
        for details on the methodology
        Will align traces based on the dtt computed on the trace stack

        Parameters
        ----------
        stream : stream of data. Prefer to use traces windowed around the phase of interest


        Returns
        -------
        reference_waveform : np.array of same length as the individual traces in stream
        stream_aligned : stream
                        the input stream but with individual traces aligned with each other

        """
        # =========================================================================
        # Align individual seismograms with the stack, and construct data matrix
        # =========================================================================
        data_matrix = np.zeros((len(stream),stream[0].stats.npts),dtype = float)
        i = 0
        stream_aligned = stream.copy()
        for tr in stream_aligned:
            cc = np.correlate(stack,tr.data, mode="full")
            time_shift = int(cc.argmax() - stream[0].stats.npts + 1)

            tr.data = align_trace(tr.data, time_shift )
            data_matrix[i] = tr.data  
            i += 1
        # svd decomposition of data_matrix
        U, s, V = np.linalg.svd(data_matrix)
        # the reference waveform is the first eigenimage
        reference_wav = - np.sign(s[0])*V[0] # check with the sign !

        return reference_wav, stream_aligned




    # =========================================================================
    # Entering the function
    # =========================================================================
    process_params = to_be_processed[0]["process_params"]

    seconds_prior_arrival = process_params["seconds_prior_arrival"]
    window_length_in_sec = process_params["window_length_in_sec"]



    for comp in components:

        # =========================================================================
        # Selection on one component
        # =========================================================================
        # !!!!!! replace first_P_arrival by phase of interest to be calculated in preprocess_data and given in process_info

        file_list = [wav["output_filename"] 
                     for wav in to_be_processed 
                     if comp in wav["channel"]]
        arrival_times = [wav["first_P_arrival"] 
                         for wav in to_be_processed 
                         if comp in wav["channel"]]
        #print("%d files for %s"%(len(file_list),comp))

        # =========================================================================
        # read traces, window around phase of interest and construct a stack
        # =========================================================================
        st_win = obspy.Stream()
        stack = np.zeros(int(window_length_in_sec/process_params["dt"]),dtype=float)
        files_to_read = []
        for file_to_read, first_tt_arrival in zip(file_list, arrival_times):
            if os.path.isfile(file_to_read):
                tr = obspy.read(file_to_read)
                idx_sigwin_start = int(np.ceil((first_tt_arrival - seconds_prior_arrival) / process_params["dt"]))
                idx_sigwin_end = idx_sigwin_start + int(window_length_in_sec/ process_params["dt"])
                tr[0].data = tr[0].data[idx_sigwin_start:idx_sigwin_end]
                tr.detrend("linear")

                stack += tr[0].data
                st_win += tr
                files_to_read.append(file_to_read)
        file_list = files_to_read
        del files_to_read
        #print("%d files to process"%len(file_list))

        # if no waveform selected at the previous step (snr criteria), quit the process
        if not st_win:
            print(("No %s data selected for this event, will skip the svd selection"%comp))
            continue

        stack /= len(st_win)

        # =========================================================================
        # construct reference waveforms and aligned data
        # =========================================================================
        ## For testing
        ## choose the option for align the traces: inversion of time-shift across the traces, or use the stack
        ## it appears that both method yiel similar results
        ## So I suggest to use either one of the two options
        try:
            reference_wav, st_aligned = construct_reference_wavefrom_from_dtt_inversion(st_win)
        except:
            reference_wav, st_aligned = construct_reference_wavefrom_from_dtt_stack(st_win, stack)

        # =========================================================================         
        # waveform selection based on correlation coeff with reference waveform
        # =========================================================================
        st_win_select = obspy.Stream()
        for tr in st_win:
            a = np.dot(tr.data, np.transpose(reference_wav))
            b = np.linalg.norm(tr.data,ord=None)*np.linalg.norm(reference_wav,ord=None)
            cc = a/b
            if np.abs(cc)>=cc_threshold:
                st_win_select += tr

        '''
        # plotting
        from lasif.visualization import plot_waveform_section
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec 
        fig=plt.figure()
        spec = gridspec.GridSpec(ncols=1, nrows=6, figure=fig)
        ax=fig.add_subplot(spec[0:4, 0])
        plot_waveform_section(ax,st_win,[], scale=4, colors=(.14, .59, .78))
        if st_win_select:
            plot_waveform_section(ax,st_win_select,[], scale=4,colors='k')
        ax2=fig.add_subplot(spec[5, 0])
        #ax2.plot(st_win[0].times(),stack/max(stack))
        ax2.plot(st_win[0].times(),reference_wav/max(reference_wav),'r')
        ax2.set_xlim(0,window_length_in_sec)
        plt.show()
        '''



        # =========================================================================
        # Remove non-selected files
        # =========================================================================    
        if st_win_select:
            files_to_keep = []
            for tr in st_win_select:
                loc_id = "%s.%s"%(tr.stats.network,tr.stats.station)
                files_to_keep.append([wav["output_filename"] 
                                      for wav in to_be_processed 
                                      if (loc_id in wav["station_filename"]) 
                                      and (comp in wav["channel"])][0])
            files_to_rmv = list(set(file_list) - set(files_to_keep))
            #print("%d/%d %s files to rmv"%(len(files_to_rmv),len(file_list),comp))

        else:
            #print("no selected %s files, will remove all"%comp)
            files_to_rmv = file_list


        # remove non-selected files
        for the_file in files_to_rmv:
            os.system("rm %s"%the_file)


    '''    
    """ I do not use this anymore, I check for every files
    # =========================================================================
    # find other horizontal files to be processed based on selection on Z 
    # =========================================================================        
    if st_win_select:
        E_file_list = []
        N_file_list = []
        for tr in st_win_select:
            loc_id = "%s.%s"%(tr.stats.network,tr.stats.station)
            E_file_list.append([wav["processing_info"]["output_filename"] 
                 for wav in to_be_processed 
                 if (loc_id in wav["processing_info"]["station_filename"]) 
                 and ('E' in wav["processing_info"]["channel"])][0])
            N_file_list.append([wav["processing_info"]["output_filename"] 
                 for wav in to_be_processed 
                 if (loc_id in wav["processing_info"]["station_filename"]) 
                 and ('N' in wav["processing_info"]["channel"])][0])
        print("%d E files and %d N files to process"%(len(E_file_list), len(N_file_list)))
    """
    
    # =========================================================================
    # Same process on E component
    # =========================================================================
    E_file_list = [wav["processing_info"]["output_filename"] 
                 for wav in to_be_processed 
                 if 'E' in wav["processing_info"]["channel"]]
    E_arrival_times = [wav["processing_info"]["first_P_arrival"] 
                 for wav in to_be_processed 
                 if 'E' in wav["processing_info"]["channel"]]
    
    st_win = obspy.Stream()
    stack = np.zeros(int(window_length_in_sec/process_params["dt"]),dtype=float)
    E_files = []
    for file_to_read, first_tt_arrival in zip(E_file_list, E_arrival_times):
        if os.path.isfile(file_to_read):
            tr = obspy.read(file_to_read)
            idx_sigwin_start = int(np.ceil((first_tt_arrival - seconds_prior_arrival) / process_params["dt"]))
            idx_sigwin_end = idx_sigwin_start + int(window_length_in_sec/ process_params["dt"])
            tr[0].data = tr[0].data[idx_sigwin_start:idx_sigwin_end]
            tr.detrend("linear")
            
            stack += tr[0].data
            st_win += tr
            E_files.append(file_to_read)
    E_file_list = E_files
    print("%d E files to actually process"%len(E_file_list))
       
    
    if st_win:
        stack /= len(st_win)
        try:
            reference_wav, st_aligned = construct_reference_wavefrom_from_dtt_inversion(st_win)
            del stack
        except:
            reference_wav, st_aligned = construct_reference_wavefrom_from_dtt_stack(st_win, stack)
                  
        st_win_select = obspy.Stream()
        for tr in st_win:
            a = np.dot(tr.data, np.transpose(reference_wav))
            b = np.linalg.norm(tr.data,ord=None)*np.linalg.norm(reference_wav,ord=None)
            cc = a/b
            if np.abs(cc)>=cc_threshold:
                st_win_select += tr
                
        if st_win_select:
            E_files_to_keep = []
            for tr in st_win_select:
                loc_id = "%s.%s"%(tr.stats.network,tr.stats.station)
                E_files_to_keep.append([wav["processing_info"]["output_filename"] 
                     for wav in to_be_processed 
                     if (loc_id in wav["processing_info"]["station_filename"]) 
                     and ('E' in wav["processing_info"]["channel"])][0])
            E_files_to_rmv = list(set(E_file_list) - set(E_files_to_keep))
            print("%d/%d E files to rmb"%(len(E_files_to_rmv),len(E_file_list)))
            
        else:
            print("no selected E files, will remove all")
            E_files_to_rmv = E_file_list
         
        
        # remove non-selected Efiles
        for Efile in E_files_to_rmv:
            os.system("rm %s"%Efile)
           
            
            
    # =========================================================================
    # Same process on N component
    # =========================================================================
    N_file_list = [wav["processing_info"]["output_filename"] 
                 for wav in to_be_processed 
                 if 'N' in wav["processing_info"]["channel"]]
    N_arrival_times = [wav["processing_info"]["first_P_arrival"] 
                 for wav in to_be_processed 
                 if 'N' in wav["processing_info"]["channel"]]
    
    st_win = obspy.Stream()
    stack = np.zeros(int(window_length_in_sec/process_params["dt"]),dtype=float)
    N_files = []
    for file_to_read, first_tt_arrival in zip(N_file_list, N_arrival_times):
        if os.path.isfile(file_to_read):
            tr = obspy.read(file_to_read)
            idx_sigwin_start = int(np.ceil((first_tt_arrival - seconds_prior_arrival) / process_params["dt"]))
            idx_sigwin_end = idx_sigwin_start + int(window_length_in_sec/ process_params["dt"])
            tr[0].data = tr[0].data[idx_sigwin_start:idx_sigwin_end]
            tr.detrend("linear")
            
            stack += tr[0].data
            st_win += tr
            N_files.append(file_to_read)
    N_file_list = N_files
    print("%d N files to actually process"%len(N_file_list))
       
    
    if st_win:
        stack /= len(st_win)
        try:
            reference_wav, st_aligned = construct_reference_wavefrom_from_dtt_inversion(st_win)
            del stack
        except:
            reference_wav, st_aligned = construct_reference_wavefrom_from_dtt_stack(st_win, stack)
                  
        st_win_select = obspy.Stream()
        for tr in st_win:
            a = np.dot(tr.data, np.transpose(reference_wav))
            b = np.linalg.norm(tr.data,ord=None)*np.linalg.norm(reference_wav,ord=None)
            cc = a/b
            if np.abs(cc)>=cc_threshold:
                st_win_select += tr
                
        if st_win_select:
            N_files_to_keep = []
            for tr in st_win_select:
                loc_id = "%s.%s"%(tr.stats.network,tr.stats.station)
                N_files_to_keep.append([wav["processing_info"]["output_filename"] 
                     for wav in to_be_processed 
                     if (loc_id in wav["processing_info"]["station_filename"]) 
                     and ('N' in wav["processing_info"]["channel"])][0])
            N_files_to_rmv = list(set(N_file_list) - set(N_files_to_keep))
            print("%d/%d N files to rmb"%(len(N_files_to_rmv),len(N_file_list)))
            
        else:
            print("no selected N files, will remove all")
            N_files_to_rmv = N_file_list
        
        
        # remove non-selected Efiles
        for Nfile in N_files_to_rmv:
            os.system("rm %s"%Nfile)
                
    '''

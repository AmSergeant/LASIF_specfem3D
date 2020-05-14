#!/usr/bin/env python
# -*- coding: utf-8 -*-


import copy
import itertools
import math
import numpy as np
import os

from lasif import LASIFError, LASIFNotFoundError

from .component import Component
import matplotlib.pyplot as plt

class VisualizationsComponent(Component):
    """
    Component offering project visualization. Has to be initialized fairly
    late at is requires a lot of data to be present.

    :param communicator: The communicator instance.
    :param component_name: The name of this component for the communicator.
    """


    def plot_stations(self,plot_relief=True):
        """
        Plots the domain and locations for all stations on the map.
        """
        from lasif import visualization

        stations = self.comm.query.get_all_stations()
        if plot_relief:
            m = self.comm.project.domain.plot(skip_map_features=True)
            m.shadedrelief()
            min_lat = np.round(m.latmin, decimals=1)
            max_lat = np.round(m.latmax, decimals=1)
            min_lon = np.round(m.lonmin, decimals=0)
            max_lon = np.round(m.lonmax, decimals=0)
            width = np.abs(max_lon - min_lon)
            if width > 60.0:
                stepsize = 10.0
            elif 30.0 < width <= 60.0:
                stepsize = 5.0
            elif 10.0 < width <= 30.0:
                stepsize = 2.0
            else:
                stepsize = 1.0
            parallels = np.arange(min_lat,max_lat,stepsize)
            meridians = np.arange(min_lon,max_lon,stepsize)

            m.drawparallels(parallels, labels=[1,0,0,0], color='white', linewidth = 0.5)
            m.drawmeridians(meridians, labels=[0,0,0,1], color='white')
        else:
            m = self.comm.project.domain.plot()
        scatter_stations = visualization.plot_stations(map_object=m, station_dict=stations)
        fig = plt.gcf()
        ax = plt.gca()
        annot = ax.annotate("", xy=(0,0), xytext=(-10,10),textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w"),
                            arrowprops=dict(arrowstyle="->"))
        annot.set_visible(False)
        station_names = []
        for station_name in stations:
            station_names.append(station_name)

        def hover_stations(event):        
            vis = annot.get_visible()
            cont, ind = scatter_stations.contains(event)
            if cont:
                pos = scatter_stations.get_offsets()[ind["ind"][0]]
                annot.xy = pos
                text = [station_names[n] for n in ind["ind"]][0]
                annot.set_text(text)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

        fig.canvas.mpl_connect("button_press_event", hover_stations) 


    def plot_stations_for_event(self,event_name, iteration_name = "raw", plot_relief=True, color = 'red', ax=None):
        """
            Plots the domain and locations for all stations available for one event on the map.
            if no iteration_name given, will plot all stations available in raw data folder for this event
            else take stations avaiable in processed data folder
        """
        from lasif import visualization

        if plot_relief=="True":
            m = self.comm.project.domain.plot(skip_map_features=True, ax=ax)
            m.shadedrelief()
            min_lat = np.round(m.latmin, decimals=1)
            max_lat = np.round(m.latmax, decimals=1)
            min_lon = np.round(m.lonmin, decimals=0)
            max_lon = np.round(m.lonmax, decimals=0)
            width = np.abs(max_lon - min_lon)
            if width > 60.0:
                stepsize = 10.0
            elif 30.0 < width <= 60.0:
                stepsize = 5.0
            elif 10.0 < width <= 30.0:
                stepsize = 2.0
            else:
                stepsize = 1.0
                parallels = np.arange(min_lat,max_lat,stepsize)
                meridians = np.arange(min_lon,max_lon,stepsize)
                m.drawparallels(parallels, labels=[1,0,0,0], color='white', linewidth = 0.5)
                m.drawmeridians(meridians, labels=[0,0,0,1], color='white')
        else:
            m = self.comm.project.domain.plot(ax=ax)

        # Get the event and extract information from it.
        event_info = self.comm.events.get(event_name)
        if iteration_name == "raw":
            stations = self.comm.query.get_all_stations_for_event(event_name)
        else:
            stations = self.comm.query.get_all_stations_for_event_for_iteration(event_name, iteration_name)

        visualization.plot_stations_for_event(
            map_object=m, station_dict=stations,
            event_info=event_info, color = color)

        return m



    def plot_events(self, plot_type="map", azimuthal_projection=False, iteration_name = "raw", events = None, ax=None):
        """
        Plots the domain and beachballs for all events on the map.

        :param plot_type: Determines the type of plot created.
            * ``map`` (default) - a map view of the events
            * ``depth`` - a depth distribution histogram
            * ``time`` - a time distribution histogram
            * ``azimuth`` - an azimuth distribution histogram
        """
        from lasif import visualization
        import matplotlib.pyplot as plt
        if not ax:
            ax = plt.subplot(111)

        # get config as defined in the config file
        proj = self.comm.project
        ds = proj.config["download_settings"]
        config = ds["configuration"]


        # Get a dictionary containing all events
        if events is None:
            if iteration_name =='raw':
                events = list(self.comm.events.get_all_events().values())
            else:
                events = self.comm.query.get_all_events_for_processed_data(iteration_name)

        if plot_type == "map":
            beachballsize=0.02
            if config == "teleseismic":
                if azimuthal_projection is True:
                    m = self.comm.project.domain.plot(Teleseismic=True, azimuthal_projection=True, ax=ax)
                    beachballsize=0.04
                else:
                    m = self.comm.project.domain.plot(Teleseismic=True, azimuthal_projection=False, ax=ax)
            else:
                m = self.comm.project.domain.plot(ax=ax)

            visualization.plot_events(events, map_object=m, beachball_size=beachballsize)
            import matplotlib.pyplot as plt
            plt.title("%d events"%len(events))
            #stations = self.comm.query.get_all_stations()
            #visualization.plot_stations(map_object=m, station_dict=stations)
            return m

        elif plot_type == "depth":
            visualization.plot_event_histogram(events, "depth")
        elif plot_type == "time":
            visualization.plot_event_histogram(events, "time")
        elif plot_type == "azimuth":
            lat_0 = self.comm.project.domain.center.latitude
            lon_0 = self.comm.project.domain.center.longitude
            visualization.plot_event_histogram(events, "azimuth", lat_0 = lat_0, lon_0 = lon_0)
        else:
            msg = "Unknown plot_type"
            raise LASIFError(msg)

    def plot_beachball_arrivals(self,event_info, stations, Phases=['P'], beachcolor="red",ax=None):
        """
            plot seismic arrivals for phases on a beachball plot
        """
        
        import matplotlib.pyplot as plt
        if not ax:
            ax = plt.subplot(111)
        from lasif.tools.pyrocko import beachball, cake, orthodrome
        ax.set_axis_off()
        ax.set_xlim(-1.1, 1.1)
        ax.set_ylim(-1.1, 1.1)
        ax.set_aspect('equal')
        mt = [event_info["m_tt"], event_info["m_pp"], event_info["m_rr"], -event_info["m_tp"],
              event_info["m_rt"], -event_info["m_rp"]]
        beachball.plot_beachball_mpl(
            mt, ax,
            position=(0, 0),
            size=2.0,
            color_t=beachcolor,
            linewidth=1.,
            projection='lambert',
            size_units='data')

        # earth model and phase for takeoff angle computations
        mod = cake.load_model('ak135-f-average.vf')
        slat = event_info["latitude"]
        slon = event_info["longitude"]
        rdepth = 0.
        sdepth = event_info["depth_in_km"]*1e3
        lats = []
        lngs = []
        for station in stations:
            lats.append(stations[station]["latitude"])
            lngs.append(stations[station]["longitude"])
        for i, Phase in enumerate(Phases):
            xcoor = []
            ycoor = []
            if i>2:
                import random
                r = random.random()
                b = random.random()
                g = random.random()
                color = (r, g, b)
            elif i == 0:
                color = 'black'
            elif i == 1:
                color = 'blue'
            for rlat, rlon in zip(lats,lngs):
                distance = orthodrome.distance_accurate50m(slat, slon, rlat, rlon)
                rays = mod.arrivals(
                    phases=cake.PhaseDef(Phase),
                    zstart=sdepth, zstop=rdepth, distances=[distance*cake.m2d])
                if not rays:
                    continue      
                takeoff = rays[0].takeoff_angle()
                azi = orthodrome.azimuth(slat, slon, rlat, rlon)        
                # to spherical coordinates, r, theta, phi in radians
                rtp = np.array([[1., np.deg2rad(takeoff), np.deg2rad(90.-azi)]])        
                # to 3D coordinates (x, y, z)
                points = beachball.numpy_rtp2xyz(rtp)        
                # project to 2D with same projection as used in beachball
                x, y = beachball.project(points, projection='lambert').T

                xcoor.append(x)
                ycoor.append(y)
                ax.plot(x, y, '+', ms=10., mew=2.0, mec=color, mfc='none')
                #ax.text(x, y, Phase, color = 'blue', horizontalalignment='center', verticalalignment='center')

            annot = ax.annotate("", xy=(np.mean(xcoor),np.mean(ycoor)), xytext=(-20,-20),textcoords="offset points",
                                        arrowprops=dict(arrowstyle="-", color = color), color = color)
            annot.set_text(Phase)
            annot.set_visible(True)
            
            
            
    def plot_event(self, event_name, config="local", azimuthal_projection=False, iteration_name = "raw", Phases = ['P'], event_info=None, ax=None):
        """
        Plots information about one event on the map.
        """
        import matplotlib.pyplot as plt
        if not ax:
            ax = plt.subplot(111)

        if event_info is None and not self.comm.events.has_event(event_name):
            self.comm.events.update_cache()
            if not self.comm.events.has_event(event_name):
                msg = "Event '%s' not found in project." % event_name
                raise ValueError(msg)

        if config == "teleseismic":
            if azimuthal_projection is True:
                map_object = self.comm.project.domain.plot(Teleseismic=True, azimuthal_projection =True, ax=ax)
            else:
                map_object = self.comm.project.domain.plot(Teleseismic=True, ax=ax)
        else:
            map_object = self.comm.project.domain.plot(ax=ax)

        from lasif import visualization

        # Get the event and extract information from it.
        if event_info is None:
            event_info = self.comm.events.get(event_name)

        # Get a dictionary containing all stations that have data for the
        # current event.
        if event_info is None:
            if iteration_name =='raw':
                stations = self.comm.query.get_all_stations_for_event(event_name)
            else:
                stations = self.comm.query.get_all_stations_for_event_for_iteration(event_name, iteration_name)
        else:
            stations = self.comm.query.get_all_stations()
            
        # Plot the stations if it has some. This will also plot raypaths.
        scatter_stations = visualization.plot_stations_for_event(
            map_object=map_object, station_dict=stations,
            event_info=event_info)

        # Plot the beachball for one event.
        visualization.plot_events(events=[event_info], map_object=map_object)

        # print station names when clicking on it

        fig=plt.gcf()
        annot = ax.annotate("", xy=(0,0), xytext=(-10,10),textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w"),
                            arrowprops=dict(arrowstyle="->"))
        annot.set_visible(False)
        station_names = []
        for station_name in stations:
            station_names.append(station_name)

        def hover_stations(event):        
            vis = annot.get_visible()
            cont, ind = scatter_stations.contains(event)
            if cont:
                pos = scatter_stations.get_offsets()[ind["ind"][0]]
                annot.xy = pos
                text = [station_names[n] for n in ind["ind"]][0]
                annot.set_text(text)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

        fig.canvas.mpl_connect("button_press_event", hover_stations)

        # add beachball with arrival in subaxes
        if config == "local":      
            ax_sub = fig.add_axes([0.05, 0.63, 0.3, 0.25])
            self.comm.visualizations.plot_beachball_arrivals(event_info, stations, Phases=Phases, beachcolor="red", ax=ax_sub)
        return map_object



    def plot_domain(self, plot_simulation_domain=True):
        """
        Plots the simulation domain and the actual physical domain.
        """
        self.comm.project.domain.plot(
            plot_simulation_domain=plot_simulation_domain)

    def plot_raydensity(self, save_plot=True, plot_stations=True):
        """
        Plots the raydensity.
        """
        from lasif import visualization
        import matplotlib.pyplot as plt

        plt.figure(figsize=(20, 21))

        m = self.comm.project.domain.plot()

        event_stations = []
        for event_name, event_info in \
                self.comm.events.get_all_events().items():
            try:
                stations = \
                    self.comm.query.get_all_stations_for_event(event_name)
            except LASIFError:
                stations = {}
            event_stations.append((event_info, stations))

        visualization.plot_raydensity(map_object=m,
                                      station_events=event_stations,
                                      domain=self.comm.project.domain)

        visualization.plot_events(list(self.comm.events.get_all_events().values()),
                                  map_object=m)

        if plot_stations:
            stations = itertools.chain.from_iterable((
                list(_i[1].values()) for _i in event_stations if _i[1]))
            # Remove duplicates
            stations = [(_i["latitude"], _i["longitude"]) for _i in stations]
            stations = set(stations)
            x, y = m([_i[1] for _i in stations],
                     [_i[0] for _i in stations])
            m.scatter(x, y, s=14 ** 2, color="#333333",
                      edgecolor="#111111", alpha=0.6, zorder=200,
                      marker="v")

        plt.tight_layout()

        if save_plot:
            outfile = os.path.join(
                self.comm.project.get_output_folder(
                    type="raydensity_plots", tag="raydensity"),
                "raydensity.png")
            plt.savefig(outfile, dpi=200, transparent=True)
            print("Saved picture at %s" % outfile)


    def plot_raw_waveforms(self, event_name, components = ['E','N','Z'], 
                           scaling = 0.5, Filter=False, freqmin=0.01, freqmax=0.1,
                           plot_arrival=True, Phase='P'):
        """
        Will plot raw seismic waveform gather with data stored in DATA/event_name/raw folder 
        that corresponds to preprocessing parameters of the iteration file
        Parameters
        ----------
        event_name : str
            name of the event to plot
        components: list of str
            list of seismic components to plot
        Filter: Boolean, True or False
            if True, will filter the waveforms 
        freqmin: float
            minimum corner frequency
        freqmax: float
            maximum corner frequency
        plot_window: Boolean, True or False
            if True, will plot the time limits for the pase windowing as defined in the iteration file
        Phase: str
            Name of the seismic phase to window. 

        Raises
        ------
        ValueError
            - if no seismic waveforms stored in raw folder 
            for one of the asked component
            - if no waveforms or the epicentral distance range is zero

        Returns
        -------
        None.

        """  
        from lasif.visualization import plot_waveform_section
        from obspy.geodetics.base import locations2degrees
        from obspy.core import read, Stream
        if plot_arrival:
            from obspy.taup import TauPyModel
            earth_model = TauPyModel("ak135")

        # First check the given components

        # check the number of components
        ncomp = len(components)
        if ncomp>3:
            msg = ("There are more than 3 components given")
            raise LASIFError(msg)

        # sort components so will plot Z data first
        if 'Z' in components:
            components.sort(reverse=False)
            components = np.roll(np.array(components),1).tolist()
        if (('R' in components )or ('T' in components)):
            print("====================================================================================")
            print("!!! No R or T component available for raw data !!!")

        if Filter == "True":
            print(("Will bandpass filter the waveforms between %0.2f and %0.2f Hz, This might take a while ..."%(freqmin,freqmax)))

        # check data component and initiate titles for each component plot
        titles=[]
        for comp in components:
            if comp=='Z':
                titles.append('Vertical')
            elif comp=='E':
                titles.append('East')
            elif comp=='N':
                titles.append('North')
            else: 
                raise ValueError("Invalid data component '%s'. Component should be E, N, R, T or Z." %comp)

        # Get event and iteration infos
        event = self.comm.events.get(event_name)

        # Get station and waveform infos
        station_coordinates = self.comm.query.get_all_stations_for_event(event_name)
        '''
        for station in station_coordinates:
            if station_coordinates[station]["latitude"] is None:
                cha = station.split(".")
                cha = cha[0]+'.'+cha[1]
                coord = self.comm.inventory_db.get_coordinates(cha)
                station_coordinates[station]["latitude"] = coord["latitude"]
                station_coordinates[station]["longitude"] = coord["longitude"]
        '''

        waveforms = self.comm.waveforms.get_metadata_raw(event_name)
        # select waveforms for requested components:
        waveforms = [wav for wav in waveforms 
                     if wav["channel"][-1] in components]
        if not waveforms:
            raise ValueError("No raw available waveforms for '%s'." %event_name)

        # Group by station name.
        def func(x):
            return ".".join(x["channel_id"].split(".")[:2])
        waveforms.sort(key=func)

        # First step is to calculate all epicentral distances.
        for wav in waveforms:
            station = station_coordinates[wav["network"]+'.'+wav["station"]]
            wav["epicentral_distance"] = locations2degrees(
                event["latitude"], event["longitude"], 
                station["latitude"],station["longitude"])
            if plot_arrival :
                tts = earth_model.get_travel_times(source_depth_in_km=event["depth_in_km"],
                                                   distance_in_degree=wav["epicentral_distance"],
                                                   phase_list=[Phase])
                if len(tts):
                    wav["tt_arrival_time"]  = tts[0].time
                else:
                    wav["tt_arrival_time"] = np.nan

        min_epicentral_distance = math.ceil(min(
            _i["epicentral_distance"] for _i in waveforms))
        max_epicentral_distance = math.ceil(max(
            _i["epicentral_distance"] for _i in waveforms))
        epicentral_range = max_epicentral_distance - min_epicentral_distance
        if epicentral_range == 0:
            raise ValueError

        # Plotting
        import matplotlib.pylab as plt
        import matplotlib.gridspec as gridspec   
        fig = plt.figure(figsize=(12, 8))
        spec = gridspec.GridSpec(ncols=ncomp, nrows=1, figure=fig)
        if Filter == "True":
            plt.suptitle("%s - %0.2f-%0.2f Hz"%(event_name,freqmin,freqmax))
        else:
            plt.suptitle(event_name)

        # subplot on each component: read data files + plotting
        sub = 0
        ax1 = fig.add_subplot(spec[0, sub])
        comp = components[sub]
        file_list = [wav["filename"] for wav in waveforms if comp in wav["channel"]]
        if len(file_list)==0:
            if ncomp == 1:
                raise LASIFError("No raw data available for component '%s'" %comp)
            else:
                print(("No raw data available for component '%s'" %comp))
        else:
            offset = [wav["epicentral_distance"] for wav in waveforms if comp in wav["channel"]]
            tt_arrival_times = [wav["tt_arrival_time"] for wav in waveforms if comp in wav["channel"]]
            st = Stream()
            for files, tt_arrival in zip(file_list, tt_arrival_times):
                tr = read(files)
                
                if Filter == "True":
                    tr.detrend("linear")
                    tr.detrend("demean")
                    tr.taper(0.05, type="cosine")
                    tr.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=3,
                              zerophase=False)
                    tr.detrend("linear")
                    tr.detrend("demean")
                    tr.taper(0.05, type="cosine")
                    tr.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=3,
                              zerophase=False)
                tr
                st += tr
            plot_waveform_section(ax1,st, offset, reftime = event["origin_time"], scale=scaling)
            if plot_arrival:
                for arrival, dist in zip(tt_arrival_times, offset):
                    plt.plot(arrival, dist,'r.',markersize=2)

        ax1.set_ylim(min_epicentral_distance-0.1*epicentral_range,
                     max_epicentral_distance+0.1*epicentral_range)
        ax1.set_title(titles[sub])
        ax1.set_xlabel('Time since event (s)')
        ax1.set_ylabel('Epicentral distance (degree)')
        annot1 = ax1.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                              bbox=dict(boxstyle="round", fc="w"),
                              arrowprops=dict(arrowstyle="->"))
        annot1.set_visible(False)            


        if ncomp>1:
            sub += 1
            ax2 = fig.add_subplot(spec[0, sub], sharex = ax1, sharey=ax1)
            comp = components[sub]
            file_list = [wav["filename"] for wav in waveforms if comp in wav["channel"]]
            if len(file_list)==0:
                print(("No raw data available for component '%s'" %comp))
            else:
                offset = [wav["epicentral_distance"] for wav in waveforms if comp in wav["channel"]]
                tt_arrival_times = [wav["tt_arrival_time"] for wav in waveforms if comp in wav["channel"]]
                st = Stream()
                for files, tt_arrival in zip(file_list, tt_arrival_times):
                    tr = read(files)
                    
                    if Filter == "True":
                        tr.detrend("linear")
                        tr.detrend("demean")
                        tr.taper(0.05, type="cosine")
                        tr.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=3,
                                  zerophase=False)
                        tr.detrend("linear")
                        tr.detrend("demean")
                        tr.taper(0.05, type="cosine")
                        tr.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=3,
                                  zerophase=False)
                    st += tr
                plot_waveform_section(ax2,st, offset, reftime = event["origin_time"], scale=scaling)
                if plot_arrival:
                    for arrival, dist in zip(tt_arrival_times, offset):
                        plt.plot(arrival, dist,'r.',markersize=2)

            #ax2.set_ylim(min_epicentral_distance, max_epicentral_distance)
            ax2.set_title(titles[sub])
            ax2.set_xlabel('Time since event (s)')
            annot2 = ax2.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                                  bbox=dict(boxstyle="round", fc="w"),
                                  arrowprops=dict(arrowstyle="->"))
            annot2.set_visible(False) 


            if ncomp>2:
                sub += 1
                ax3 = fig.add_subplot(spec[0, sub], sharex = ax1, sharey=ax1)
                comp = components[sub]
                file_list = [wav["filename"] for wav in waveforms if comp in wav["channel"]]
                if len(file_list)==0:
                    print(("No raw data available for component '%s'" %comp))
                else:
                    offset = [wav["epicentral_distance"] for wav in waveforms if comp in wav["channel"]]
                    tt_arrival_times = [wav["tt_arrival_time"] for wav in waveforms if comp in wav["channel"]]
                    st = Stream()
                    for files, tt_arrival in zip(file_list, tt_arrival_times):
                        tr = read(files)
                        
                        if Filter == "True":
                            tr.detrend("linear")
                            tr.detrend("demean")
                            tr.taper(0.05, type="cosine")
                            tr.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=3,
                                      zerophase=False)
                            tr.detrend("linear")
                            tr.detrend("demean")
                            tr.taper(0.05, type="cosine")
                            tr.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=3,
                                      zerophase=False)
                        st += tr
                    plot_waveform_section(ax3,st, offset, reftime = event["origin_time"], scale=scaling)
                    if plot_arrival:
                        for arrival, dist in zip(tt_arrival_times, offset):
                            plt.plot(arrival, dist,'r.',markersize=2)

                #ax3.set_ylim(min_epicentral_distance, max_epicentral_distance)
                ax3.set_title(titles[sub])
                ax3.set_xlabel('Time since event (s)')
                annot3 = ax3.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                                      bbox=dict(boxstyle="round", fc="w"),
                                      arrowprops=dict(arrowstyle="->"))
                annot3.set_visible(False) 

        def on_plot_hover(event):
            # Iterating over each data member plotted
            if event.inaxes == ax1:
                xpos = ax1.get_xlim()[0]
                for curve in ax1.get_lines():
                    # Searching which data member corresponds to current mouse position
                    if curve.contains(event)[0]:
                        annot1.xy = (xpos,curve.get_data()[1][0])
                        annot1.set_text(curve.get_label())
                        annot1.set_visible(True)
                        fig.canvas.draw_idle()
            elif len(fig.axes)>1 and event.inaxes == ax2:
                xpos = ax2.get_xlim()[0]
                for curve in ax2.get_lines():
                    # Searching which data member corresponds to current mouse position
                    if curve.contains(event)[0]:
                        annot2.xy = (xpos,curve.get_data()[1][0])
                        annot2.set_text(curve.get_label())
                        annot2.set_visible(True)
                        fig.canvas.draw_idle()
            elif len(fig.axes)>2 and event.inaxes == ax3:
                xpos = ax3.get_xlim()[0]
                for curve in ax3.get_lines():
                    # Searching which data member corresponds to current mouse position
                    if curve.contains(event)[0]:
                        annot3.xy = (xpos,curve.get_data()[1][0])
                        annot3.set_text(curve.get_label())
                        annot3.set_visible(True)
                        fig.canvas.draw_idle()

        fig.canvas.mpl_connect('button_press_event', on_plot_hover) 



    def plot_preprocessed_waveforms(self, event_name, iteration_name, components = ['E','N','Z'], 
                                    scaling = 0.5, plot_raw=False, plot_window=True, Phase='P'):
        """
        Will plot seismic waveform gather with data stored in DATA/event_name/preprocessed folder
        that corresponds to preprocessing parameters of the iteration file
        Parameters
        ----------
        event_name : str
            name of the event to plot
        iteration_name: int
            name of the iteration
        components: list of str
            list of seismic components to plot
        plot_raw: Boolean, True or False
            if True, will plot on top of preprocessed waveforms, the raw waveforms
            which did not pass the preprocessing selection process
        plot_window: Boolean, True or False
            if True, will plot the time limits for the pase windowing as defined in the iteration file
        Phase: str
            Name of the seismic phase to window. 

        Raises
        ------
        ValueError
            - if no seismic waveforms stored in preprocessed folder
            for one of the asked component
            - if no waveforms or the epicentral distance range is zero

        Returns
        -------
        None.

        """
        from lasif.visualization import plot_waveform_section
        from obspy.geodetics.base import locations2degrees
        from obspy.core import read, Stream
        if plot_window:
            from obspy.taup import TauPyModel
            earth_model = TauPyModel("ak135")

        # First check the given components

        # check the number of components
        ncomp = len(components)
        if ncomp > 3:
            msg = ("There are more than 3 components given")
            raise LASIFError(msg)

        # sort components so will plot Z data first
        if 'Z' in components:
            components.sort(reverse=False)
            components = np.roll(np.array(components), 1).tolist()
        if (('R' in components)or ('T' in components)) and plot_raw:
            print(
                "====================================================================================")
            print(
                "!!! No R or T component available for raw data, plot_raw data option is disabled !!!")
            plot_raw = False

        # check data component and initiate titles for each component plot
        titles = []
        for comp in components:
            if comp == 'Z':
                titles.append('Vertical')
            elif comp == 'E':
                titles.append('East')
            elif comp == 'N':
                titles.append('North')
            elif comp == 'R':
                titles.append('Radial')
            elif comp == 'T':
                titles.append('Transverse')
            else: 
                raise ValueError("Invalid data component '%s'. Component should be E, N, R, T or Z." %comp)

        # Get event and iteration infos
        event = self.comm.events.get(event_name)
        iteration = self.comm.iterations.get(iteration_name)
        pparam = iteration.get_process_params()
        processing_tag = iteration.processing_tag

        # Get station and waveform infos
        station_coordinates = self.comm.query.get_all_stations_for_event(event_name)
        waveforms = self.comm.waveforms.get_metadata_processed(event_name, processing_tag)
        # select waveforms for requested components:
        waveforms = [wav for wav in waveforms 
                     if wav["channel"][-1] in components]
        # Group by station name.

        def func(x):
            return ".".join(x["channel_id"].split(".")[:2])
        waveforms.sort(key=func)
        if plot_raw:
            raws = self.comm.waveforms.get_metadata_raw(event_name)
            raws.sort(key=func)
            # we keep only raw files that were not processed to save time
            raw_sta = [raw["channel_id"] for raw in raws]
            wav_sta = [wav["channel_id"] for wav in waveforms]
            not_in_wav = list(set(raw_sta) ^ set(wav_sta))
            raws = [raw for raw in raws if raw["channel_id"] in not_in_wav]

        # First step is to calculate all epicentral distances.
        for wav in waveforms:
            station = station_coordinates[wav["network"]+'.'+wav["station"]]
            wav["epicentral_distance"] = locations2degrees(
                event["latitude"], event["longitude"], 
                station["latitude"],station["longitude"])
            if plot_window :
                tts = earth_model.get_travel_times(source_depth_in_km=event["depth_in_km"],
                                                   distance_in_degree=wav["epicentral_distance"],
                                                   phase_list=[Phase])
                if len(tts):
                    wav["tt_arrival_time"]  = tts[0].time
                else:
                    wav["tt_arrival_time"] = np.nan
        if plot_raw == "True" and raws:
            for raw in raws:
                station = station_coordinates[wav["network"]+'.'+wav["station"]]
                raw["epicentral_distance"] = locations2degrees(
                    event["latitude"], event["longitude"], 
                    station["latitude"],station["longitude"])

        if plot_raw == "True" and raws:
            min_epicentral_distance = math.ceil(np.min([min(
                _i["epicentral_distance"] for _i in waveforms),
                min(_i["epicentral_distance"] for _i in raws)]))
            max_epicentral_distance = math.ceil(np.max([max(
                _i["epicentral_distance"] for _i in waveforms),
                max(_i["epicentral_distance"] for _i in raws)]))
        else:
            min_epicentral_distance = math.ceil(min(
                _i["epicentral_distance"] for _i in waveforms))
            max_epicentral_distance = math.ceil(max(
                _i["epicentral_distance"] for _i in waveforms))
        epicentral_range = max_epicentral_distance - min_epicentral_distance
        if epicentral_range == 0:
            raise ValueError

        # Plotting
        import matplotlib.pylab as plt
        import matplotlib.gridspec as gridspec
        fig = plt.figure(figsize=(12, 8))
        spec = gridspec.GridSpec(ncols=ncomp, nrows=1, figure=fig)
        plt.suptitle(event_name)

        # subplot on each component: read data files + plotting
        sub = 0
        ax1 = fig.add_subplot(spec[0, sub])
        comp = components[sub]
        file_list = [wav["filename"] for wav in waveforms if comp in wav["channel"]]
        if len(file_list)==0:
            if ncomp == 1:
                raise LASIFError("No preprocessed data available for component '%s'" %comp)
            else:
                print(("No preprocessed data available for component '%s'" %comp))
        else:
            offset = [wav["epicentral_distance"] for wav in waveforms if comp in wav["channel"]]
            tt_arrival_times = [wav["tt_arrival_time"] for wav in waveforms if comp in wav["channel"]]
            st = Stream()
            #corrected_arrival_times = []
            for files, tt_arrival in zip(file_list, tt_arrival_times):
                tr = read(files)
                #delta = event["origin_time"].datetime - tr[0].stats.starttime.datetime
                #corrected_arrival_times.append(tt_arrival + delta.total_seconds())
                st += tr
            if plot_raw == "True" and raws:
                file_list_raw = [raw["filename"] for raw in raws if comp in raw["channel"]]
                if len(file_list_raw)==0:
                    print(("No raw data available for component '%s'" %comp))
                else:
                    offset_raw = [raw["epicentral_distance"] for raw in raws if comp in raw["channel"]]
                    st_raw = Stream()
                    print(("Bandpass filtering raw data for components %s"%comp))
                    for i, files in enumerate(file_list_raw):
                        tr = read(files)
                        tr.detrend("linear")
                        tr.detrend("demean")
                        tr.taper(0.05, type="cosine")
                        tr.filter(
                            "bandpass",
                            freqmin=pparam["highpass"],
                            freqmax=pparam["lowpass"],
                            corners=3,
                            zerophase=False)
                        tr.detrend("linear")
                        tr.detrend("demean")
                        tr.taper(0.05, type="cosine")
                        tr.filter(
                            "bandpass",
                            freqmin=pparam["highpass"],
                            freqmax=pparam["lowpass"],
                            corners=3,
                            zerophase=False)
                        tr.interpolate(
                            sampling_rate=1.0 /
                            pparam["dt"],
                            method="lanczos",
                            starttime=event["origin_time"],
                            window="blackman",
                            a=12,
                            npts=pparam["npts"])
                        st_raw += tr
                    plot_waveform_section(ax1,st_raw, offset_raw, reftime = event["origin_time"], scale=scaling, 
                                          colors=(.14, .59, .78))

            plot_waveform_section(ax1,st, offset, reftime = event["origin_time"], scale=scaling)
            if plot_window:
                for arrival, dist in zip(tt_arrival_times, offset):
                    plt.plot(arrival-pparam["seconds_prior_arrival"], dist,'r.',markersize=2)
                    plt.plot(arrival-pparam["seconds_prior_arrival"]+pparam["window_length_in_sec"], dist,'r.', markersize=2)

        ax1.set_ylim(min_epicentral_distance-0.1*epicentral_range,
                     max_epicentral_distance+0.1*epicentral_range)
        ax1.set_title(titles[sub])
        ax1.set_xlabel('Time since event (s)')
        ax1.set_ylabel('Epicentral distance (degree)')
        annot1 = ax1.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                              bbox=dict(boxstyle="round", fc="w"),
                              arrowprops=dict(arrowstyle="->"))
        annot1.set_visible(False)            


        if ncomp>1:
            sub += 1
            ax2 = fig.add_subplot(spec[0, sub], sharex = ax1, sharey=ax1)
            comp = components[sub]
            file_list = [wav["filename"] for wav in waveforms if comp in wav["channel"]]
            if len(file_list)==0:
                print(("No preprocessed data available for component '%s'" %comp))
            else:
                offset = [wav["epicentral_distance"] for wav in waveforms if comp in wav["channel"]]
                tt_arrival_times = [wav["tt_arrival_time"] for wav in waveforms if comp in wav["channel"]]
                st = Stream()
                #corrected_arrival_times = []
                for files, tt_arrival in zip(file_list, tt_arrival_times):
                    tr = read(files)
                    #delta = event["origin_time"].datetime - tr[0].stats.starttime.datetime
                    #corrected_arrival_times.append(tt_arrival + delta.total_seconds())
                    st += tr

                if plot_raw == "True" and raws:
                    file_list_raw = [raw["filename"] for raw in raws if comp in raw["channel"]]
                    if len(file_list_raw)==0:
                        print(("No raw data available for component '%s'" %comp))
                    else:
                        offset_raw = [raw["epicentral_distance"] for raw in raws if comp in raw["channel"]]
                        st_raw = Stream()
                        print(("Bandpass filtering raw data for components %s"%comp))
                        for files in file_list_raw:
                            tr = read(files)
                            tr.detrend("linear")
                            tr.detrend("demean")
                            tr.taper(0.05, type="cosine")
                            tr.filter("bandpass", freqmin=pparam["highpass"], freqmax=pparam["lowpass"], corners=3,
                                      zerophase=False)
                            tr.detrend("linear")
                            tr.detrend("demean")
                            tr.taper(0.05, type="cosine")
                            tr.filter("bandpass", freqmin=pparam["highpass"], freqmax=pparam["lowpass"], corners=3,
                                      zerophase=False)
                            tr.interpolate(
                                sampling_rate=1.0 / pparam["dt"],
                                method="lanczos", starttime=event["origin_time"], window="blackman", a=12,
                                npts=pparam["npts"])
                            st_raw += tr
                        plot_waveform_section(ax2,st_raw, offset_raw, reftime = event["origin_time"], scale=scaling, 
                                              colors=(.14, .59, .78))

                plot_waveform_section(ax2,st, offset, reftime = event["origin_time"], scale=scaling)
                if plot_window:
                    for arrival, dist in zip(tt_arrival_times, offset):
                        plt.plot(arrival-pparam["seconds_prior_arrival"], dist,'r.',markersize=2)
                        plt.plot(arrival-pparam["seconds_prior_arrival"]+pparam["window_length_in_sec"], dist,'r.', markersize=2)

            #ax2.set_ylim(min_epicentral_distance, max_epicentral_distance)
            ax2.set_title(titles[sub])
            ax2.set_xlabel('Time since event (s)')
            annot2 = ax2.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                                  bbox=dict(boxstyle="round", fc="w"),
                                  arrowprops=dict(arrowstyle="->"))
            annot2.set_visible(False) 


            if ncomp>2:
                sub += 1
                ax3 = fig.add_subplot(spec[0, sub], sharex = ax1, sharey=ax1)
                comp = components[sub]
                file_list = [wav["filename"] for wav in waveforms if comp in wav["channel"]]
                if len(file_list)==0:
                    print(("No preprocessed data available for component '%s'" %comp))
                else:
                    offset = [wav["epicentral_distance"] for wav in waveforms if comp in wav["channel"]]
                    tt_arrival_times = [wav["tt_arrival_time"] for wav in waveforms if comp in wav["channel"]]
                    st = Stream()
                    #corrected_arrival_times = []
                    for files, tt_arrival in zip(file_list, tt_arrival_times):
                        tr = read(files)
                        #delta = event["origin_time"].datetime - tr[0].stats.starttime.datetime
                        #corrected_arrival_times.append(tt_arrival + delta.total_seconds())
                        st += tr

                    if plot_raw == "True" and raws:
                        file_list_raw = [raw["filename"] for raw in raws if comp in raw["channel"]]
                        if len(file_list_raw)==0:
                            print(("No raw data available for component '%s'" %comp))
                        else:
                            offset_raw = [raw["epicentral_distance"] for raw in raws if comp in raw["channel"]]
                            st_raw = Stream()
                            print(("Bandpass filtering raw data for components %s"%comp))
                            for files in file_list_raw:
                                tr = read(files)
                                tr.detrend("linear")
                                tr.detrend("demean")
                                tr.taper(0.05, type="cosine")
                                tr.filter("bandpass", freqmin=pparam["highpass"], freqmax=pparam["lowpass"], corners=3,
                                          zerophase=False)
                                tr.detrend("linear")
                                tr.detrend("demean")
                                tr.taper(0.05, type="cosine")
                                tr.filter("bandpass", freqmin=pparam["highpass"], freqmax=pparam["lowpass"], corners=3,
                                          zerophase=False)
                                tr.interpolate(
                                    sampling_rate=1.0 / pparam["dt"],
                                    method="lanczos", starttime=event["origin_time"], window="blackman", a=12,
                                    npts=pparam["npts"])
                                st_raw += tr
                            plot_waveform_section(ax3,st_raw, offset_raw, reftime = event["origin_time"], scale=scaling, 
                                                  colors=(.14, .59, .78))

                    plot_waveform_section(ax3,st, offset, reftime = event["origin_time"], scale=scaling)
                    if plot_window:
                        for arrival, dist in zip(tt_arrival_times, offset):
                            plt.plot(arrival-pparam["seconds_prior_arrival"], dist,'r.',markersize=2)
                            plt.plot(arrival-pparam["seconds_prior_arrival"]+pparam["window_length_in_sec"], dist,'r.', markersize=2)

                #ax3.set_ylim(min_epicentral_distance, max_epicentral_distance)
                ax3.set_title(titles[sub])
                ax3.set_xlabel('Time since event (s)')
                annot3 = ax3.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                                      bbox=dict(boxstyle="round", fc="w"),
                                      arrowprops=dict(arrowstyle="->"))
                annot3.set_visible(False) 

        def on_plot_hover(event):
            # Iterating over each data member plotted
            if event.inaxes == ax1:
                xpos = ax1.get_xlim()[0]
                for curve in ax1.get_lines():
                    # Searching which data member corresponds to current mouse position
                    if curve.contains(event)[0]:
                        annot1.xy = (xpos,curve.get_data()[1][0])
                        annot1.set_text(curve.get_label())
                        annot1.set_visible(True)
                        fig.canvas.draw_idle()
            elif len(fig.axes)>1 and event.inaxes == ax2:
                xpos = ax2.get_xlim()[0]
                for curve in ax2.get_lines():
                    # Searching which data member corresponds to current mouse position
                    if curve.contains(event)[0]:
                        annot2.xy = (xpos,curve.get_data()[1][0])
                        annot2.set_text(curve.get_label())
                        annot2.set_visible(True)
                        fig.canvas.draw_idle()
            elif len(fig.axes)>2 and event.inaxes == ax3:
                xpos = ax3.get_xlim()[0]
                for curve in ax3.get_lines():
                    # Searching which data member corresponds to current mouse position
                    if curve.contains(event)[0]:
                        annot3.xy = (xpos,curve.get_data()[1][0])
                        annot3.set_text(curve.get_label())
                        annot3.set_visible(True)
                        fig.canvas.draw_idle()

        fig.canvas.mpl_connect('button_press_event', on_plot_hover) 


    def plot_synthetic_waveforms(self, event_name, iteration_name, components = ['E','N','Z'], 
                                 scaling = 0.5, plot_window=True, Phase='P'):
        """
        Will plot seismic waveform gather with synthetics and preprocessed data
        that corresponds to preprocessing parameters of the iteration file
        Parameters
        ----------
        event_name : str
            name of the event to plot
        iteration_name: int
            name of the iteration
        components: list of str
            list of seismic components to plot
        plot_window: Boolean, True or False
            if True, will plot the time limits for the pase windowing as defined in the iteration file
        Phase: str
            Name of the seismic phase to window. 


        Raises
        ------
        ValueError
            - if no seismic waveforms stored in the synthetic folder
            for one of the asked component
            - if no waveforms or the epicentral distance range is zero

        Returns
        -------
        None.

        """
        from lasif.visualization import plot_waveform_section
        from obspy.geodetics.base import locations2degrees
        from obspy.core import read, Stream
        if plot_window:
            from obspy.taup import TauPyModel
            earth_model = TauPyModel("ak135")

        # First check the given components

        # check the number of components
        ncomp = len(components)
        if ncomp > 3:
            msg = ("There are more than 3 components given")
            raise LASIFError(msg)

        # sort components so will plot Z data first
        if 'Z' in components:
            components.sort(reverse=False)
            components = np.roll(np.array(components),1).tolist()

        # check data component and initiate titles for each component plot
        titles = []
        for comp in components:
            if comp == 'Z':
                titles.append('Vertical')
            elif comp == 'E':
                titles.append('East')
            elif comp == 'N':
                titles.append('North')
            elif comp == 'R':
                titles.append('Radial')
            elif comp == 'T':
                titles.append('Transverse')
            else: 
                raise ValueError("Invalid data component '%s'. Component should be E, N, R, T or Z." %comp)

        # Get event and iteration infos
        event = self.comm.events.get(event_name)
        iteration = self.comm.iterations.get(iteration_name)
        pparam = iteration.get_process_params()
        processing_tag = iteration.processing_tag

        # Get station and waveform infos
        station_coordinates = self.comm.query.get_all_stations_for_event(event_name)
        waveforms = self.comm.waveforms.get_metadata_processed(event_name, processing_tag)
        synthetics = self.comm.waveforms.get_metadata_synthetic(event_name, iteration_name)
        # select waveforms for requested components:
        waveforms = [wav for wav in waveforms 
                     if wav["channel"][-1] in components]
        synthetics = [syn for syn in synthetics 
                      if syn["channel"][-1] in components]
        # Group by station name.

        def func(x):
            return ".".join(x["channel_id"].split(".")[:2])
        waveforms.sort(key=func)
        synthetics.sort(key=func)


        # First step is to calculate all epicentral distances.   
        #for (wav,syn) in zip(waveforms, synthetics):
        for wav in waveforms:
            station = station_coordinates[wav["network"]+'.'+wav["station"]]
            wav["epicentral_distance"] = locations2degrees(
                event["latitude"], event["longitude"], 
                station["latitude"],station["longitude"])
            if plot_window :
                tts = earth_model.get_travel_times(source_depth_in_km=event["depth_in_km"],
                                                   distance_in_degree=wav["epicentral_distance"],
                                                   phase_list=[Phase])
                if len(tts):
                    wav["tt_arrival_time"]  = tts[0].time
                else:
                    wav["tt_arrival_time"] = np.nan

        for syn in synthetics:
            syn["epicentral_distance"] = [wav["epicentral_distance"]
                                          for wav in waveforms 
                                          if syn["station"] in wav["station"]][0]
            #syn["epicentral_distance"] = wav["epicentral_distance"].copy()
        min_epicentral_distance = math.ceil(min(
            _i["epicentral_distance"] for _i in waveforms))
        max_epicentral_distance = math.ceil(max(
            _i["epicentral_distance"] for _i in waveforms))
        epicentral_range = max_epicentral_distance - min_epicentral_distance
        if epicentral_range == 0:
            raise ValueError

        # Plotting
        import matplotlib.pylab as plt
        import matplotlib.gridspec as gridspec
        fig = plt.figure(figsize=(12, 8))
        spec = gridspec.GridSpec(ncols=ncomp, nrows=1, figure=fig)
        plt.suptitle(event_name)

        # subplot on each component: read data files + plotting
        sub = 0
        ax1 = fig.add_subplot(spec[0, sub])
        comp = components[sub]
        file_list = [syn["filename"] for syn in synthetics if comp in syn["channel"]]
        if len(file_list)==0:
            raise LASIFError("No synthetics data available for component '%s'" %comp)        
        offset = [syn["epicentral_distance"] for syn in synthetics if comp in syn["channel"]]
        tt_arrival_times = [wav["tt_arrival_time"] for wav in waveforms if comp in wav["channel"]]
        st = Stream()
        #corrected_arrival_times = []
        for files, tt_arrival in zip(file_list, tt_arrival_times):
            tr = read(files)
            #delta = event["origin_time"].datetime - tr[0].stats.starttime.datetime
            #corrected_arrival_times.append(tt_arrival + delta.total_seconds())
            st += tr     
        file_list_wav = [wav["filename"] for wav in waveforms if comp in wav["channel"]]
        if len(file_list_wav)==0:
            raise LASIFError("No processed data available for component '%s'" %comp)
        offset_wav = [wav["epicentral_distance"] for wav in waveforms if comp in wav["channel"]]
        st_wav = Stream()
        for files in file_list_wav:
            st_wav += read(files)
        plot_waveform_section(ax1,st, offset, reftime = event["origin_time"], scale=scaling,colors='r')
        plot_waveform_section(ax1,st_wav, offset_wav, reftime = event["origin_time"], scale=scaling, 
                              colors='k')  
        if plot_window:
            for arrival, dist in zip(tt_arrival_times, offset_wav):
                plt.plot(arrival-pparam["seconds_prior_arrival"], dist,'b.',markersize=2)
                plt.plot(arrival-pparam["seconds_prior_arrival"]+pparam["window_length_in_sec"], dist,'b.', markersize=2)

        ax1.set_ylim(min_epicentral_distance-0.1*epicentral_range,
                     max_epicentral_distance+0.1*epicentral_range)
        ax1.set_title(titles[sub])
        ax1.set_xlabel('Time since event (s)')
        ax1.set_ylabel('Epicentral distance (degree)')
        annot1 = ax1.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                              bbox=dict(boxstyle="round", fc="w"),
                              arrowprops=dict(arrowstyle="->"))
        annot1.set_visible(False)

        if ncomp>1:
            sub += 1
            ax2 = fig.add_subplot(spec[0, sub], sharex=ax1, sharey=ax1)
            comp = components[sub]
            file_list = [syn["filename"] for syn in synthetics if comp in syn["channel"]]
            if len(file_list)==0:
                raise LASIFError("No synthetics data available for component '%s'" %comp)        
            offset = [syn["epicentral_distance"] for syn in synthetics if comp in syn["channel"]]
            st = Stream()
            #corrected_arrival_times = []
            for files, tt_arrival in zip(file_list, tt_arrival_times):
                tr = read(files)
                #delta = event["origin_time"].datetime - tr[0].stats.starttime.datetime
                #corrected_arrival_times.append(tt_arrival + delta.total_seconds())
                st += tr
            file_list_wav = [wav["filename"] for wav in waveforms if comp in wav["channel"]]
            if len(file_list_wav)==0:
                raise LASIFError("No processed data available for component '%s'" %comp)
            offset_wav = [wav["epicentral_distance"] for wav in waveforms if comp in wav["channel"]]
            tt_arrival_times = [wav["tt_arrival_time"] for wav in waveforms if comp in wav["channel"]]
            st_wav = Stream()
            for files in file_list_wav:
                st_wav += read(files)
            plot_waveform_section(ax2,st, offset, reftime = event["origin_time"], scale=scaling,colors='r')
            plot_waveform_section(ax2,st_wav, offset_wav, reftime = event["origin_time"], scale=scaling, 
                                  colors='k')
            if plot_window:
                for arrival, dist in zip(tt_arrival_times, offset_wav):
                    plt.plot(arrival-pparam["seconds_prior_arrival"], dist,'b.',markersize=2)
                    plt.plot(arrival-pparam["seconds_prior_arrival"]+pparam["window_length_in_sec"], dist,'b.', markersize=2)

            #ax2.set_ylim(min_epicentral_distance, max_epicentral_distance)
            ax2.set_title(titles[sub])
            ax2.set_xlabel('Time since event (s)')
            annot2 = ax2.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                                  bbox=dict(boxstyle="round", fc="w"),
                                  arrowprops=dict(arrowstyle="->"))
            annot2.set_visible(False) 

            if ncomp>2:
                sub += 1
                ax3 = fig.add_subplot(spec[0, sub], sharex=ax1, sharey=ax1)
                comp = components[sub]
                file_list = [syn["filename"] for syn in synthetics if comp in syn["channel"]]
                if len(file_list)==0:
                    raise LASIFError("No synthetics data available for component '%s'" %comp)        
                offset = [syn["epicentral_distance"] for syn in synthetics if comp in syn["channel"]]
                tt_arrival_times = [wav["tt_arrival_time"] for wav in waveforms if comp in wav["channel"]]
                st = Stream()
                #corrected_arrival_times = []
                for files, tt_arrival in zip(file_list, tt_arrival_times):
                    tr = read(files)
                    #delta = event["origin_time"].datetime - tr[0].stats.starttime.datetime
                    #corrected_arrival_times.append(tt_arrival + delta.total_seconds())
                    st += tr     
                file_list_wav = [wav["filename"] for wav in waveforms if comp in wav["channel"]]
                if len(file_list_wav)==0:
                    raise LASIFError("No processed data available for component '%s'" %comp)
                offset_wav = [wav["epicentral_distance"] for wav in waveforms if comp in wav["channel"]]
                st_wav = Stream()
                for files in file_list_wav:
                    st_wav += read(files)
                plot_waveform_section(ax3,st, offset, reftime = event["origin_time"], scale=scaling,colors='r')
                plot_waveform_section(ax3,st_wav, offset_wav, reftime = event["origin_time"], scale=scaling, 
                                      colors='k')
                if plot_window:
                    for arrival, dist in zip(tt_arrival_times, offset_wav):
                        plt.plot(arrival-pparam["seconds_prior_arrival"], dist,'b.',markersize=2)
                        plt.plot(arrival-pparam["seconds_prior_arrival"]+pparam["window_length_in_sec"], dist,'b.', markersize=2)

                #ax3.set_ylim(min_epicentral_distance, max_epicentral_distance)
                ax3.set_title(titles[sub])
                ax3.set_xlabel('Time since event (s)')
                annot3 = ax3.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                                      bbox=dict(boxstyle="round", fc="w"),
                                      arrowprops=dict(arrowstyle="->"))
                annot3.set_visible(False) 

        def on_plot_hover(event):
            # Iterating over each data member plotted
            if event.inaxes == ax1:
                xpos = ax1.get_xlim()[0]
                for curve in ax1.get_lines():
                    # Searching which data member corresponds to current mouse position
                    if curve.contains(event)[0]:
                        annot1.xy = (xpos,curve.get_data()[1][0])
                        annot1.set_text(curve.get_label())
                        annot1.set_visible(True)
                        fig.canvas.draw_idle()
            elif len(fig.axes)>1 and event.inaxes == ax2:
                xpos = ax2.get_xlim()[0]
                for curve in ax2.get_lines():
                    # Searching which data member corresponds to current mouse position
                    if curve.contains(event)[0]:
                        annot2.xy = (xpos,curve.get_data()[1][0])
                        annot2.set_text(curve.get_label())
                        annot2.set_visible(True)
                        fig.canvas.draw_idle()
            elif len(fig.axes)>2 and event.inaxes == ax3:
                xpos = ax3.get_xlim()[0]
                for curve in ax3.get_lines():
                    # Searching which data member corresponds to current mouse position
                    if curve.contains(event)[0]:
                        annot3.xy = (xpos,curve.get_data()[1][0])
                        annot3.set_text(curve.get_label())
                        annot3.set_visible(True)
                        fig.canvas.draw_idle()

        fig.canvas.mpl_connect('button_press_event', on_plot_hover)


    def plot_synthetic_from_stf(self, event_name, iteration_name, components = ['E','N','Z'], 
                                scaling = 0.5, plot_raw = True, plot_window=True, Phase='P'):
        """
        Will plot seismic waveform gather with synthetics computed prior and after the stf deconvolution
        Parameters
        ----------
        event_name : str
            name of the event to plot
        iteration_name: int
            name of the iteration
        components: list of str
            list of seismic components to plot
        plot_window: Boolean, True or False
            if True, will plot the time limits for the pase windowing as defined in the iteration file
        Phase: str
            Name of the seismic phase to window. 


        Raises
        ------
        ValueError
            - if no seismic waveforms stored in the synthetic folder
            for one of the asked component
            - if no stf
            - if no waveforms or the epicentral distance range is zero

        Returns
        -------
        None.

        """
        from lasif.visualization import plot_waveform_section
        from obspy.geodetics.base import locations2degrees
        from obspy.core import read, Stream
        if plot_window:
            from obspy.taup import TauPyModel
            earth_model = TauPyModel("ak135")

        def compute_synthetics_from_stf(src_array, stream_green):
            nfft = stream_green[0].stats.npts
            src_fft = np.fft.fft(src_array, nfft)
            stream_syn = stream_green.copy()
            for i, sy in enumerate(stream_green):
                sy_fft = np.fft.fft(sy.data, nfft)
                cal = sy.copy()
                cal.data = np.real(np.fft.ifft(src_fft * sy_fft))
                stream_syn[i] = cal
            return stream_syn

        # First check the given components

        # check the number of components
        ncomp = len(components)
        if ncomp > 3:
            msg = ("There are more than 3 components given")
            raise LASIFError(msg)

        # sort components so will plot Z data first
        if 'Z' in components:
            components.sort(reverse=False)
            components = np.roll(np.array(components),1).tolist()

        # check data component and initiate titles for each component plot
        titles = []
        for comp in components:
            if comp == 'Z':
                titles.append('Vertical')
            elif comp == 'E':
                titles.append('East')
            elif comp == 'N':
                titles.append('North')
            elif comp == 'R':
                titles.append('Radial')
            elif comp == 'T':
                titles.append('Transverse')
            else: 
                raise ValueError("Invalid data component '%s'. Component should be E, N, R, T or Z." %comp)

        # Get event and iteration infos
        event = self.comm.events.get(event_name)
        iteration = self.comm.iterations.get(iteration_name)
        pparam = iteration.get_process_params()
        processing_tag = iteration.processing_tag

        # Get station and waveform infos
        station_coordinates = self.comm.query.get_all_stations_for_event(event_name)
        waveforms = self.comm.waveforms.get_metadata_processed(event_name, processing_tag)
        greens = self.comm.waveforms.get_metadata_synthetic(event_name, iteration_name)
        # select waveforms for requested components:
        waveforms = [wav for wav in waveforms 
                     if wav["channel"][-1] in components]
        greens = [green for green in greens 
                  if green["channel"][-1] in components]
        # Group by station name.

        def func(x):
            return ".".join(x["channel_id"].split(".")[:2])
        waveforms.sort(key=func)
        greens.sort(key=func)


        # First step is to calculate all epicentral distances.
        stream_syn = Stream()  
        for wav in waveforms:
            station = station_coordinates[wav["network"]+'.'+wav["station"]]
            wav["epicentral_distance"] = locations2degrees(
                event["latitude"], event["longitude"], 
                station["latitude"],station["longitude"])
            if plot_window :
                tts = earth_model.get_travel_times(source_depth_in_km=event["depth_in_km"],
                                                   distance_in_degree=wav["epicentral_distance"],
                                                   phase_list=[Phase])
                if len(tts):
                    wav["tt_arrival_time"]  = tts[0].time
                else:
                    wav["tt_arrival_time"] = np.nan
            try:
                tr = self.comm.waveforms.get_synthetic_from_stf(event_name, wav["channel_id"], iteration_name)
                stream_syn += tr
            except:
                print(('no stf syn for %s'%wav["channel_id"]))
                continue

        for syn in greens:
            syn["epicentral_distance"] = [wav["epicentral_distance"]
                                          for wav in waveforms 
                                          if syn["station"] in wav["station"]][0]
        min_epicentral_distance = math.ceil(min(
            _i["epicentral_distance"] for _i in waveforms))
        max_epicentral_distance = math.ceil(max(
            _i["epicentral_distance"] for _i in waveforms))
        epicentral_range = max_epicentral_distance - min_epicentral_distance
        if epicentral_range == 0:
            raise ValueError



        # Plotting
        import matplotlib.pylab as plt
        import matplotlib.gridspec as gridspec
        fig = plt.figure(figsize=(12, 8))
        spec = gridspec.GridSpec(ncols=ncomp, nrows=1, figure=fig)
        plt.suptitle(event_name)

        # subplot on each component: read data files + plotting
        sub = 0
        ax1 = fig.add_subplot(spec[0, sub])
        comp = components[sub]

        stf = \
            self.comm.waveforms.get_waveform_stf(event_name, iteration_name, component = 'Z')
        if not stf:
            raise LASIFError("No stf data available for component '%s'" %comp)   
        starttime = stf[0].stats.starttime
        endtime = stf[0].stats.endtime

        file_list_green = [syn["filename"] for syn in greens if comp in syn["channel"]]
        if len(file_list_green)==0:
            raise LASIFError("No synthetics data available for component '%s'" %comp)        
        offset_green = [syn["epicentral_distance"] for syn in greens if comp in syn["channel"]]
        st_green = Stream()
        for files, offset in zip(file_list_green, offset_green):
            tr = read(files)
            tr.trim(starttime,endtime)
            st_green += tr  
        st_syn = stream_syn.select(component = comp)
        offset_syn = []
        for tr in st_syn:
            net = tr.stats.network
            stat = tr.stats.station
            offset_syn.append([wav["epicentral_distance"]
                               for wav in waveforms 
                               if net in wav["channel_id"].split('.')[0]
                               and stat in wav["channel_id"].split('.')[1]][0])

        file_list_wav = [wav["filename"] for wav in waveforms if comp in wav["channel"]]
        if len(file_list_wav)==0:
            raise LASIFError("No processed data available for component '%s'" %comp)
        offset_wav = [wav["epicentral_distance"] for wav in waveforms if comp in wav["channel"]]
        tt_arrival_times = [wav["tt_arrival_time"] for wav in waveforms if comp in wav["channel"]]
        st_wav = Stream()
        for files in file_list_wav:
            tr = read(files)
            tr.trim(starttime, endtime)
            st_wav += tr
        #st_syn = compute_synthetics_from_stf(stf[0].data, st_green)
        #offset_syn = offset_green

        if plot_raw == "True":
            plot_waveform_section(ax1,st_green, offset_green, scale=scaling,colors='r')
        plot_waveform_section(ax1,st_wav, offset_wav, scale=scaling, 
                              colors='k', lw =2)
        plot_waveform_section(ax1,st_syn, offset_syn, scale=scaling, 
                              colors='b')
        if plot_window:
            for arrival, dist in zip(tt_arrival_times, offset_wav):
                arrival -= starttime - event["origin_time"]
                plt.plot(arrival-pparam["seconds_prior_arrival"], dist,'r.',markersize=2)
                plt.plot(arrival-pparam["seconds_prior_arrival"]+pparam["window_length_in_sec"], dist,'r.', markersize=2)

        ax1.set_ylim(min_epicentral_distance-0.1*epicentral_range,
                     max_epicentral_distance+0.1*epicentral_range)
        ax1.set_title(titles[sub])
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Epicentral distance (degree)')
        annot1 = ax1.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                              bbox=dict(boxstyle="round", fc="w"),
                              arrowprops=dict(arrowstyle="->"))
        annot1.set_visible(False)

        if ncomp>1:
            sub += 1
            ax2 = fig.add_subplot(spec[0, sub], sharex=ax1, sharey=ax1)
            comp = components[sub]

            '''
            stf = \
                self.comm.waveforms.get_waveform_stf(event_name, iteration_name, component = comp)
            if not stf:
                raise LASIFError(
                    "No stf data available for component '%s'" %
                    comp)
            print(stf[0])
            starttime = stf[0].stats.starttime
            endtime = stf[0].stats.endtime
            '''

            file_list_green = [syn["filename"] for syn in greens if comp in syn["channel"]]
            if len(file_list_green)==0:
                raise LASIFError("No synthetics data available for component '%s'" %comp)        
            offset_green = [syn["epicentral_distance"] for syn in greens if comp in syn["channel"]]
            st_green = Stream()
            for files in file_list_green:
                tr = read(files)
                tr.trim(starttime, endtime)
                st_green += tr
            file_list_wav = [wav["filename"] for wav in waveforms if comp in wav["channel"]]
            if len(file_list_wav)==0:
                raise LASIFError("No processed data available for component '%s'" %comp)
            offset_wav = [wav["epicentral_distance"] for wav in waveforms if comp in wav["channel"]]
            tt_arrival_times = [wav["tt_arrival_time"] for wav in waveforms if comp in wav["channel"]]
            st_wav = Stream()
            for files in file_list_wav:
                tr = read(files)
                tr.trim(starttime, endtime)
                st_wav += tr
            #st_syn = compute_synthetics_from_stf(stf[0].data, st_green)
            st_syn = stream_syn.select(component = comp)
            offset_syn = []
            for tr in st_syn:
                net = tr.stats.network
                stat = tr.stats.station
                offset_syn.append([wav["epicentral_distance"]
                                   for wav in waveforms 
                                   if net in wav["channel_id"].split('.')[0]
                                   and stat in wav["channel_id"].split('.')[1]][0])
            if plot_raw == "True":
                plot_waveform_section(ax2,st_green, offset_green, scale=scaling, colors='r')
            plot_waveform_section(ax2,st_wav, offset_wav, scale=scaling, 
                                  colors='k', lw =2)
            plot_waveform_section(ax2,st_syn, offset_syn, scale=scaling, 
                                  colors='b')
            if plot_window:
                for arrival, dist in zip(tt_arrival_times, offset_wav):
                    arrival -= starttime - event["origin_time"]
                    plt.plot(arrival-pparam["seconds_prior_arrival"], dist,'r.',markersize=2)
                    plt.plot(arrival-pparam["seconds_prior_arrival"]+pparam["window_length_in_sec"], dist,'r.', markersize=2)

            #ax2.set_ylim(min_epicentral_distance, max_epicentral_distance)
            ax2.set_title(titles[sub])
            ax2.set_xlabel('Time (s)')
            annot2 = ax2.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                                  bbox=dict(boxstyle="round", fc="w"),
                                  arrowprops=dict(arrowstyle="->"))
            annot2.set_visible(False) 

            if ncomp>2:
                sub += 1
                ax3 = fig.add_subplot(spec[0, sub], sharex=ax1, sharey=ax1)
                comp = components[sub]

                '''
                stf = \
                    self.comm.waveforms.get_waveform_stf(event_name, iteration_name, component = comp)
                if not stf:
                    raise LASIFError(
                        "No stf data available for component '%s'" %
                        comp)
                print(stf[0])
                starttime = stf[0].stats.starttime
                endtime = stf[0].stats.endtime
                '''

                file_list_green = [syn["filename"] for syn in greens if comp in syn["channel"]]
                if len(file_list_green)==0:
                    raise LASIFError("No synthetics data available for component '%s'" %comp)        
                offset_green = [syn["epicentral_distance"] for syn in greens if comp in syn["channel"]]
                st_green = Stream()
                for files in file_list_green:
                    tr = read(files)
                    tr.trim(starttime, endtime)
                    st_green += tr
                file_list_wav = [wav["filename"] for wav in waveforms if comp in wav["channel"]]
                if len(file_list_wav)==0:
                    raise LASIFError("No processed data available for component '%s'" %comp)
                offset_wav = [wav["epicentral_distance"] for wav in waveforms if comp in wav["channel"]]
                tt_arrival_times = [wav["tt_arrival_time"] for wav in waveforms if comp in wav["channel"]]
                st_wav = Stream()
                for files in file_list_wav:
                    tr = read(files)
                    tr.trim(starttime, endtime)
                    st_wav += tr
                #st_syn = compute_synthetics_from_stf(stf[0].data, st_green)
                st_syn = stream_syn.select(component = comp)
                offset_syn = []
                for tr in st_syn:
                    net = tr.stats.network
                    stat = tr.stats.station
                    offset_syn.append([wav["epicentral_distance"]
                                       for wav in waveforms 
                                       if net in wav["channel_id"].split('.')[0]
                                       and stat in wav["channel_id"].split('.')[1]][0])

                if plot_raw == "True":
                    plot_waveform_section(ax3,st_green, offset_green, scale=scaling,colors='r')
                plot_waveform_section(ax3,st_wav, offset_wav, scale=scaling, 
                                      colors='k', lw =2)
                plot_waveform_section(ax3,st_syn, offset_syn, scale=scaling, 
                                      colors='b')
                if plot_window:
                    for arrival, dist in zip(tt_arrival_times, offset_wav):
                        arrival -= starttime - event["origin_time"]
                        plt.plot(arrival-pparam["seconds_prior_arrival"], dist,'r.',markersize=2)
                        plt.plot(arrival-pparam["seconds_prior_arrival"]+pparam["window_length_in_sec"], dist,'r.', markersize=2)

                #ax3.set_ylim(min_epicentral_distance, max_epicentral_distance)
                ax3.set_title(titles[sub])
                ax3.set_xlabel('Time (s)')
                annot3 = ax3.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                                      bbox=dict(boxstyle="round", fc="w"),
                                      arrowprops=dict(arrowstyle="->"))
                annot3.set_visible(False) 

        def on_plot_hover(event):
            # Iterating over each data member plotted
            if event.inaxes == ax1:
                xpos = ax1.get_xlim()[0]
                for curve in ax1.get_lines():
                    # Searching which data member corresponds to current mouse position
                    if curve.contains(event)[0]:
                        annot1.xy = (xpos,curve.get_data()[1][0])
                        annot1.set_text(curve.get_label())
                        annot1.set_visible(True)
                        fig.canvas.draw_idle()
            elif len(fig.axes)>1 and event.inaxes == ax2:
                xpos = ax2.get_xlim()[0]
                for curve in ax2.get_lines():
                    # Searching which data member corresponds to current mouse position
                    if curve.contains(event)[0]:
                        annot2.xy = (xpos,curve.get_data()[1][0])
                        annot2.set_text(curve.get_label())
                        annot2.set_visible(True)
                        fig.canvas.draw_idle()
            elif len(fig.axes)>2 and event.inaxes == ax3:
                xpos = ax3.get_xlim()[0]
                for curve in ax3.get_lines():
                    # Searching which data member corresponds to current mouse position
                    if curve.contains(event)[0]:
                        annot3.xy = (xpos,curve.get_data()[1][0])
                        annot3.set_text(curve.get_label())
                        annot3.set_visible(True)
                        fig.canvas.draw_idle()

        fig.canvas.mpl_connect('button_press_event', on_plot_hover)


    def plot_windows(self, event, iteration, distance_bins=500,
                     ax=None, show=True):
        """
        Plot all selected windows on a epicentral distance vs duration plot
        with the color encoding the selected channels. This gives a quick
        overview of how well selected the windows for a certain event and
        iteration are.

        :param event: The event.
        :param iteration: The iteration.
        :param distance_bins: The number of bins on the epicentral
            distance axis.
        :param ax: If given, it will be plotted to this ax.
        :param show: If true, ``plt.show()`` will be called before returning.
        :return: The potentially created axes object.
        """
        from obspy.geodetics.base import locations2degrees

        event = self.comm.events.get(event)
        iteration = self.comm.iterations.get(iteration)
        pparam = iteration.get_process_params()
        window_manager = self.comm.windows.get(event, iteration)

        starttime = event["origin_time"]
        duration = (pparam["npts"] - 1) * pparam["dt"]

        # First step is to calculate all epicentral distances.
        stations = copy.deepcopy(self.comm.query.get_all_stations_for_event(
            event["event_name"]))
        for s in list(stations.values()):
            s["epicentral_distance"] = locations2degrees(
                event["latitude"], event["longitude"], s["latitude"],
                s["longitude"])

        # Plot from 0 to however far it goes.
        min_epicentral_distance = 0
        max_epicentral_distance = math.ceil(max(
            _i["epicentral_distance"] for _i in list(stations.values())))
        epicentral_range = max_epicentral_distance - min_epicentral_distance

        if epicentral_range == 0:
            raise ValueError

        # Create the image that will represent the pictures in an epicentral
        # distance plot. By default everything is black.
        #
        # First dimension: Epicentral distance.
        # Second dimension: Time.
        # Third dimension: RGB tuple.
        len_time = 1000
        len_dist = distance_bins
        image = np.zeros((len_dist, len_time, 3), dtype=np.uint8)

        # Helper functions calculating the indices.
        def _time_index(value):
            frac = np.clip((value - starttime) / duration, 0, 1)
            return int(round(frac * (len_time - 1)))

        def _space_index(value):
            frac = np.clip(
                (value - min_epicentral_distance) / epicentral_range, 0, 1)
            return int(round(frac * (len_dist - 1)))

        def _color_index(channel):
            _map = {
                "Z": 2,
                "N": 1,
                "E": 0
            }
            channel = channel[-1].upper()
            if channel not in _map:
                raise ValueError
            return _map[channel]

        for channel in window_manager.list():
            station = ".".join(channel.split(".")[:2])
            for win in window_manager.get(channel):
                image[
                    _space_index(stations[station]["epicentral_distance"]),
                    _time_index(win.starttime):_time_index(win.endtime),
                    _color_index(channel)] = 255

        # From http://colorbrewer2.org/
        color_map = {
            (255, 0, 0): (228, 26, 28),  # red
            (0, 255, 0): (77, 175, 74),  # green
            (0, 0, 255): (55, 126, 184),  # blue
            (255, 0, 255): (152, 78, 163),  # purple
            (0, 255, 255): (255, 127, 0),  # orange
            (255, 255, 0): (255, 255, 51),  # yellow
            (255, 255, 255): (250, 250, 250),  # white
            (0, 0, 0): (50, 50, 50)  # More pleasent gray background
        }

        # Replace colors...fairly complex. Not sure if there is another way...
        red, green, blue = image[:, :, 0], image[:, :, 1], image[:, :, 2]
        for color, replacement in list(color_map.items()):
            image[:, :, :][(red == color[0]) & (green == color[1]) &
                           (blue == color[2])] = replacement

        def _one(i):
            return [_i / 255.0 for _i in i]

        import matplotlib.pylab as plt
        plt.style.use("ggplot")

        artists = [
            plt.Rectangle((0, 1), 1, 1, color=_one(color_map[(0, 0, 255)])),
            plt.Rectangle((0, 1), 1, 1, color=_one(color_map[(0, 255, 0)])),
            plt.Rectangle((0, 1), 1, 1, color=_one(color_map[(255, 0, 0)])),
            plt.Rectangle((0, 1), 1, 1, color=_one(color_map[(0, 255, 255)])),
            plt.Rectangle((0, 1), 1, 1, color=_one(color_map[(255, 0, 255)])),
            plt.Rectangle((0, 1), 1, 1, color=_one(color_map[(255, 255, 0)])),
            plt.Rectangle((0, 1), 1, 1,
                          color=_one(color_map[(255, 255, 255)]))
        ]
        labels = [
            "Z",
            "N",
            "E",
            "Z + N",
            "Z + E",
            "N + E",
            "Z + N + E"
        ]

        if ax is None:
            plt.figure(figsize=(16, 9))
            ax = plt.gca()

        ax.imshow(image, aspect="auto", interpolation="nearest", vmin=0,
                  vmax=255, origin="lower")
        ax.grid()
        ax.set_title("Selected windows for iteration %s and event %s" % (
                     iteration.name, event["event_name"]))

        ax.legend(artists, labels, loc="lower right",
                  title="Selected Components")

        # Set the x-ticks.
        xticks = []
        for time in ax.get_xticks():
            # They are offset by -0.5.
            time += 0.5
            # Convert to actual time
            frac = time / float(len_time)
            time = frac * duration
            xticks.append("%.1f" % time)
        ax.set_xticklabels(xticks)
        ax.set_xlabel("Time since event in seconds")

        yticks = []
        for dist in ax.get_yticks():
            # They are offset by -0.5.
            dist += 0.5
            # Convert to actual epicentral distance.
            frac = dist / float(len_dist)
            dist = min_epicentral_distance + (frac * epicentral_range)
            yticks.append("%.1f" % dist)
        ax.set_yticklabels(yticks)
        ax.set_ylabel("Epicentral distance in degree [Binned in %i distances]"
                      % distance_bins)

        if show:
            plt.tight_layout()
            plt.show()
            plt.close()

        return ax

    def plot_window_statistics(self, iteration, ax=None, show=True,
                               cache=True):
        """
        Plots the statistics of windows for one iteration.

        :param iteration: The chosen iteration.
        :param ax: If given, it will be plotted to this ax.
        :param show: If true, ``plt.show()`` will be called before returning.
        :param cache: Use the cache for the statistics.

        :return: The potentially created axes object.
        """
        # Get the statistics.
        data = self.comm.windows.get_window_statistics(iteration=iteration,
                                                       cache=cache)

        import matplotlib
        import matplotlib.pylab as plt
        import seaborn as sns

        if ax is None:
            plt.figure(figsize=(20, 12))
            ax = plt.gca()

        ax.invert_yaxis()

        pal = sns.color_palette("Set1", n_colors=4)

        total_count = []
        count_z = []
        count_n = []
        count_e = []
        event_names = []

        width = 0.2
        ind = np.arange(len(data))

        cm = matplotlib.cm.RdYlGn

        for _i, event in enumerate(sorted(data.keys())):
            d = data[event]
            event_names.append(event)
            total_count.append(d["total_station_count"])
            count_z.append(d["stations_with_vertical_windows"])
            count_n.append(d["stations_with_north_windows"])
            count_e.append(d["stations_with_east_windows"])

            frac = int(round(100 * d["stations_with_windows"] /
                             float(d["total_station_count"])))

            color = cm(frac / 70.0)

            ax.text(-10, _i, "%i%%" % frac,
                    fontdict=dict(fontsize="x-small", ha='right', va='top'),
                    bbox=dict(boxstyle="round", fc=color, alpha=0.8))

        ax.barh(ind, count_z, width, color=pal[0],
                label="Stations with Vertical Component Windows")
        ax.barh(ind + 1 * width, count_n, width, color=pal[1],
                label="Stations with North Component Windows")
        ax.barh(ind + 2 * width, count_e, width, color=pal[2],
                label="Stations with East Component Windows")
        ax.barh(ind + 3 * width, total_count, width, color='0.4',
                label="Total Stations")

        ax.set_xlabel("Station Count")

        ax.set_yticks(ind + 2 * width, event_names)
        ax.yaxis.set_tick_params(pad=30)
        ax.set_ylim(len(data), -width)

        ax.legend(frameon=True)

        plt.suptitle("Window Statistics for Iteration %i" % 1)

        plt.tight_layout()
        plt.subplots_adjust(top=0.95)

        if show:
            plt.show()

        return ax

    def plot_data_and_synthetics(self, event, iteration, channel_id, ax=None,
                                 show=True):
        """
        Plots the data and corresponding synthetics for a given event,
        iteration, and channel.

        :param event: The event.
        :param iteration: The iteration.
        :param channel_id: The channel id.
        :param ax: If given, it will be plotted to this ax.
        :param show: If true, ``plt.show()`` will be called before returning.
        :return: The potentially created axes object.
        """
        import matplotlib.pylab as plt

        data = self.comm.query.get_matching_waveforms(event, iteration,
                                                      channel_id)
        if ax is None:
            plt.figure(figsize=(15, 3))
            ax = plt.gca()

        iteration = self.comm.iterations.get(iteration)

        ax.plot(data.data[0].times(), data.data[0].data, color="black",
                label="observed")
        ax.plot(data.synthetics[0].times(), data.synthetics[0].data,
                color="red",
                label="synthetic, iteration %s" % str(iteration.name))
        ax.legend()

        ax.set_xlabel("Seconds since event")
        ax.set_ylabel("m/s")
        ax.set_title(channel_id)
        ax.grid()

        if iteration.scale_data_to_synthetics:
            ax.text(0.995, 0.005, "data scaled to synthetics",
                    horizontalalignment="right", verticalalignment="bottom",
                    transform=ax.transAxes, color="0.2")

        if show:
            plt.tight_layout()
            plt.show()
            plt.close()

        return ax

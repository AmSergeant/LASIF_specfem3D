#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The main LASIF console script.

It is important to import necessary things at the method level to make
importing this file as fast as possible. Otherwise using the command line
interface feels sluggish and slow.


All functions starting with "lasif_" will automatically be available as
subcommands to the main "lasif" command. A decorator to determine the category
of a function is provided.

The help for every function can be accessed either via

lasif help CMD_NAME

or

lasif CMD_NAME --help


The former will be converted to the later and each subcommand is responsible
for handling the --help argument.


Each function will be passed a parser and args. It is the function author's
responsibility to add any arguments and call

parser.parse_args(args)

when done. See the existing functions for some examples. This architecture
should scale fairly well and makes it trivial to add new methods.


:copyright:
    Lion Krischer (krischer@geophysik.uni-muenchen.de), 2013
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""
import os
import lasif
from lasif import LASIFError

os.environ["OPENBLAS_NUM_THREADS"] = "1"

import argparse
import collections
import colorama
import difflib
import itertools
import progressbar
import sys
import time
import traceback
import warnings

from mpi4py import MPI

from lasif import LASIFNotFoundError
from lasif.components.project import Project

# Try to disable the ObsPy deprecation warnings. This makes LASIF work with
# the latest ObsPy stable and the master.
try:
    # It only exists for certain ObsPy versions.
    from obspy.core.util.deprecation_helpers import ObsPyDeprecationWarning
except:
    pass
else:
    warnings.filterwarnings("ignore", category=ObsPyDeprecationWarning)


FCT_PREFIX = "lasif_"


# Documentation for the subcommand groups. This will appear in the CLI
# documentation.
COMMAND_GROUP_DOCS = {
    "Data Acquisition": (
        "These functions are used to acquire and archive different types of "
        "data usually from webservices."
    ),
    "Event Management": (
        "Function helping in organzing the earthquakes inside a LASIF "
        "project."
    ),
    "Iteration Management": (
        "Functions dealing with one or more iterations inside a LASIF "
        "project."
    ),
    "Misc": (
        "All functions that do not fit in one of the other categories."
    ),
    "Misc": (
        "All functions that do not fit in one of the other categories."
    ),
    "Plotting": (
        "Functions producing pictures."
    ),
    "Project Management": (
        "Functions dealing with LASIF projects as a whole."
    )
}


def command_group(group_name):
    """
    Decorator to be able to logically group commands.
    """
    def wrapper(func):
        func.group_name = group_name
        return func
    return wrapper


def mpi_enabled(func):
    """
    Decorator to mark function with mpi capabilities.
    """
    func._is_mpi_enabled = True
    return func


class LASIFCommandLineException(Exception):
    pass


def _find_project_comm(folder, read_only_caches):
    """
    Will search upwards from the given folder until a folder containing a
    LASIF root structure is found. The absolute path to the root is returned.
    """
    max_folder_depth = 10
    folder = folder
    for _ in range(max_folder_depth):
        if os.path.exists(os.path.join(folder, "config.xml")):
            return Project(
                os.path.abspath(folder),
                read_only_caches=read_only_caches).get_communicator()
        folder = os.path.join(folder, os.path.pardir)
    msg = "Not inside a LASIF project."
    raise LASIFCommandLineException(msg)


def _find_project_comm_mpi(folder, read_only_caches):
    """
    Parallel version. Will open the caches for rank 0 with write access,
    caches from the other ranks can only read.

    :param folder: The folder were to start the search.
    :param read_only_caches: Read-only caches for rank 0. All others will
        always be read-only.
    """
    if MPI.COMM_WORLD.rank == 0:
        # Rank 0 can write the caches, the others cannot. The
        # "--read_only_caches" flag overwrites this behaviour.
        comm = _find_project_comm(folder, read_only_caches=read_only_caches)

    # Open the caches for the other ranks after rank zero has opened it to
    # allow for the initial caches to be written.
    MPI.COMM_WORLD.barrier()

    if MPI.COMM_WORLD.rank != 0:
        comm = _find_project_comm(folder, read_only_caches=True)

    return comm


def split(container, count):
    """
    Simple and elegant function splitting a container into count
    equal chunks.

    Order is not preserved but for the use case at hand this is
    potentially an advantage as data sitting in the same folder thus
    have a higher at being processed at the same time thus the disc
    head does not have to jump around so much. Of course very
    architecture dependent.
    """
    return [container[_i::count] for _i in range(count)]


@command_group("Plotting")
def lasif_plot_domain(parser, args):
    """
    Plot the project's domain on a map.
    """
    parser.add_argument("--no_simulation_domain",
                        help="Don't plot the simulation domain",
                        action="store_false")
    args = parser.parse_args(args)
    comm = _find_project_comm(".", args.read_only_caches)

    comm.visualizations.plot_domain(
        plot_simulation_domain=args.no_simulation_domain)

    import matplotlib.pyplot as plt
    plt.show()


@command_group("Misc")
def lasif_shell(parser, args):
    """
    Drops you into a shell with an active communicator instance.
    """
    args = parser.parse_args(args)

    comm = _find_project_comm(".", args.read_only_caches)
    print("LASIF shell, 'comm' object is available in the local namespace.\n")
    print(comm)
    from IPython import embed
    embed(display_banner=False)


@command_group("Plotting")
def lasif_plot_raw_waveforms(parser, args):
    """
    Plot seismic waveform gather resulting from data preprocessing steps
    """
    parser.add_argument("event_name", help="name of the event to plot waveform gather")
    parser.add_argument("--components", default="ENZ", help="list of components to plot, examples: ENZ, or Z or RTZ")
    parser.add_argument("--filter", default="False", choices=["False", "True"],
                        help="For bandpass filtering the waveforms, False by default")
    parser.add_argument("--freqmin", default=0.01,
                        help="Float for minimum frequency corner, 0.01 by default")
    parser.add_argument("--freqmax", default=0.1,
                        help="Float for maximum frequency corner, 0.01 by default")
    parser.add_argument("--scaling", default=0.5,
                        help="Float for scaling the waveforms, 0.5 by default")
    parser.add_argument("--plot_arrival", default="True", choices=["False", "True"],
                        help="For additionally plotting the seismic arrivals for a specific phase, False by default")
    parser.add_argument("--phase", default=None,
                        help="Name of the seismic phase you want to mark for the plot_arrival option"
                        "default: <phase_of_interest> in config file")

    args = parser.parse_args(args)
    event_name = args.event_name
    comps = args.components
    Filter = args.filter
    freqmin = float(args.freqmin)
    freqmax = float(args.freqmax)
    scale = float(args.scaling)
    components=[]
    for comp in comps:
        components.append(comp)
    plot_arrival = args.plot_arrival
    phase = args.phase

    comm = _find_project_comm(".", args.read_only_caches)
    comm.visualizations.plot_raw_waveforms(event_name, 
                                           components, scaling=scale, 
                                           Filter=Filter, freqmin=freqmin, freqmax=freqmax, 
                                           plot_arrival=plot_arrival, Phase=phase)

    import matplotlib.pyplot as plt
    plt.show()


@command_group("Plotting")
def lasif_plot_preprocessed_waveforms(parser, args):
    """
    Plot seismic waveform gather resulting from data preprocessing steps
    """
    parser.add_argument("iteration_name", help="name of the iteration")
    parser.add_argument("event_name", help="name of the event to plot waveform gather")
    parser.add_argument("--components", default="ENZ", help="list of components to plot, examples: ENZ, or Z or RTZ")
    parser.add_argument("--raw", default="False", choices=["False", "True"],
                        help="For additionally plotting non-selected raw waveforms on top, False by default")
    parser.add_argument("--scaling", default=0.5,
                        help="Float for scaling the waveforms, 0.5 by default")
    parser.add_argument("--plot_window", default="True", choices=["False", "True"],
                        help="For additionally plotting the phase window defined in the iteration, False by default")
    parser.add_argument("--phase", default=None,
                        help="Name of the seismic phase you want to window for the plot_window option"
                        "default: <phase_of_interest> in config file")

    args = parser.parse_args(args)
    iteration_name = args.iteration_name
    event_name = args.event_name
    comps = args.components
    raw = args.raw
    scale = float(args.scaling)
    components=[]
    for comp in comps:
        components.append(comp)
    plot_window = args.plot_window
    phase = args.phase

    comm = _find_project_comm(".", args.read_only_caches)
    comm.visualizations.plot_preprocessed_waveforms(event_name, iteration_name, 
                                                    components, scaling=scale, plot_raw=raw, 
                                                    plot_window=plot_window, Phase=phase)

    import matplotlib.pyplot as plt
    plt.show()


@command_group("Plotting")
def lasif_plot_synthetic_waveforms(parser, args):
    """
    Plot seismic waveform gather for preprocessed data and synthetic waveforms
    """
    parser.add_argument("iteration_name", help="name of the iteration")
    parser.add_argument("event_name", help="name of the event to plot waveform gather")
    parser.add_argument("--components", default="ENZ", help="list of components to plot, examples: ENZ, or Z or RTZ")
    parser.add_argument("--scaling", default=0.5,
                        help="Float for scaling the waveforms, 0.5 by default")
    parser.add_argument("--plot_window", default="True", choices=["False", "True"],
                        help="For additionally plotting the phase window defined in the iteration, False by default")
    parser.add_argument("--phase", default=None,
                        help="Name of the seismic phase you want to window for the plot_window option"
                        "default: <phase_of_interest> in config file")

    args = parser.parse_args(args)
    iteration_name = args.iteration_name
    event_name = args.event_name
    scale = float(args.scaling)
    comps = args.components
    components=[]
    for comp in comps:
        components.append(comp)
    plot_window = args.plot_window
    phase = args.phase

    comm = _find_project_comm(".", args.read_only_caches)
    comm.visualizations.plot_synthetic_waveforms(event_name, iteration_name, 
                                                 components, scaling = scale, 
                                                 plot_window=plot_window, Phase=phase)

    import matplotlib.pyplot as plt
    plt.show()


@command_group("Plotting")
def lasif_plot_event(parser, args):
    """
    Plot a single event including stations on a map.
    """
    taup_phases = ['P', 'pP', 'sP', 'PcP', 'PP', 'PKiKP', 'sPKiKP', 'S', 'pS', 'SP',
                       'sS', 'PS', 'SKS', 'SKKS', 'ScS', 'SKiKP', 'pSKS', 'sSKS', 'SS',
                       'PKIKKIK', 'SKIKKIKP', 'PKIKKIKS', 'SKIKKIKS', 'PKIKPPKIKP', 'PKPPKP',
                       'PKPPKP', 'SKIKSSKIKS']
    
    parser.add_argument("event_name", help="name of the event to plot")
    parser.add_argument("--config", default="local", choices=["local", "teleseismic"],
                        help="the type of map plot. "
                        "``local``: for map bounded to the domain, "
                        "``teleseismic``: for global map for teleseismic configuration, ")
    parser.add_argument("--azimuthal_proj", default="False", choices=["True", "False"],
                        help="the type of map projection. "
                        "``False``: for classic map projection, by default"
                        "``True``: for Azimuthal Equidistant Projection centered on the domain, ")
    parser.add_argument("--iteration", default="raw", 
                        help="iteration_name. "
                        "None: will plot stations for available waveforms in the data raw folder, "
                        "FLOAT: will plot stations for available waveforms in the data preprocessed folder that corresponds to the iteration processing tag, ")
    parser.add_argument("--phase", default=None, 
                        help="Names of seismic phases you want to plot on the beachball, option enabled only for local configuration. "
                        "default: <phase_of_interest> in config file "
                        "example: P,S,pP"
                        "List of seismic phases:%s"%(', ').join(taup_phases))
    args = parser.parse_args(args)
    event_name = args.event_name
    config=args.config
    azimuthal_proj=args.azimuthal_proj
    iteration_name = args.iteration
    phases = args.phase
    if phases is None:
        Phases = phases
    else:
        phases = phases.split(',')
        Phases=[]
        for phase in phases:
            Phases.append(phase)


    comm = _find_project_comm(".", args.read_only_caches)

    comm.visualizations.plot_event(event_name, config, azimuthal_proj, iteration_name, Phases)

    import matplotlib.pyplot as plt
    plt.show()


@command_group("Plotting")
def lasif_plot_events(parser, args):
    """
    Plot all events.

    type can be one of:
        * ``map`` (default) - a map view of the events
        * ``depth`` - a depth distribution histogram
        * ``time`` - a time distribution histogram
    """
    parser.add_argument("--type", default="map", choices=["map", "depth",
                                                          "time", "azimuth"],
                        help="the type of plot. "
                        "``map``: beachballs on a map, "
                        "``depth``: depth distribution histogram, "
                        "``time``: time distribution histogram,"
                        "``azimuth``: azimuth distribution histogram")
    parser.add_argument("--azimuthal_proj", default=False, choices=["True", "False"],
                        help="the type of map projection. "
                        "``False``: for classic map projection, by default"
                        "``True``: for Azimuthal Equidistant Projection centered on the domain, ")
    parser.add_argument("--iteration", default="raw", 
                        help="iteration_name. "
                        "None: will all found events, "
                        "FLOAT: will plot only events for available waveforms in the data preprocessed folder that corresponds to the iteration processing tag, ")
    args = parser.parse_args(args)
    plot_type = args.type
    azimuthal_proj=args.azimuthal_proj
    iteration_name = args.iteration
    # some trick
    if azimuthal_proj == "True":
        azimuthal_proj = True
        
    comm = _find_project_comm(".", args.read_only_caches)
    comm.visualizations.plot_events(plot_type, azimuthal_projection=azimuthal_proj, iteration_name = iteration_name)

    import matplotlib.pyplot as plt
    plt.show()



@command_group("Plotting")
def lasif_plot_stations(parser,args):
    """
    Plot a single event including stations on a map.
    """
    parser.add_argument("--relief", default=True, choices=["True","False"],
                        help="For plotting relief as a background image. "
                        "True: relief background, by default"
                        "False: no relief background, ")
    parser.add_argument("--color_per_network", default="True", choices=["True","False"],
                        help="For using a color code based on seismic network. "
                        "True: use network color code, by default"
                        "False: no color code, plot in red triangles")
    args = parser.parse_args(args)
    plot_type = args.relief
    color_code = args.color_per_network
    if plot_type == "False":
        plot_type = False
    if color_code == "False":
        color_code = False
    comm = _find_project_comm(".", args.read_only_caches)
    comm.visualizations.plot_stations(plot_type, color_code)

    import matplotlib.pyplot as plt
    plt.show()    


@command_group("Plotting")
def lasif_plot_raydensity(parser, args):
    """
    Plot a binned raycoverage plot for all events.
    """
    parser.add_argument("--plot_stations", help="also plot the stations",
                        action="store_true")
    args = parser.parse_args(args)

    comm = _find_project_comm(".", args.read_only_caches)
    comm.visualizations.plot_raydensity(plot_stations=args.plot_stations)


@command_group("Data Acquisition")
def lasif_add_spud_event(parser, args):
    """
    Add an event from the IRIS SPUD webservice to the project.
    """
    parser.add_argument("url", help="any SPUD momenttensor URL")
    args = parser.parse_args(args)
    url = args.url

    from lasif.scripts.iris2quakeml import iris2quakeml

    comm = _find_project_comm(".", args.read_only_caches)
    iris2quakeml(url, comm.project.paths["events"])



@command_group("Plotting")
def lasif_query_and_plot_available_stations(parser, args):
    """
    Query and plot data availability for stations in the chosen domain 
    and for a specific time range
    """
    parser.add_argument("start_date", type=str,
                        help="start date from which to query stations"
                        "example: 2011-01-01")
    parser.add_argument("end_date", type=str,
                        help="end date from which to query stations")
    

    args = parser.parse_args(args)
    start_date = args.start_date
    end_date = args.end_date
    from obspy.core import UTCDateTime
    s_yyyy, s_mm, s_dd = start_date.split('-') 
    start_date = UTCDateTime(int(s_yyyy), int(s_mm), int(s_dd), 0,0,0)
    e_yyyy, e_mm, e_dd = end_date.split('-') 
    end_date = UTCDateTime(int(e_yyyy), int(e_mm), int(e_dd), 0,0,0)
    
    comm = _find_project_comm(".", args.read_only_caches)
    
    
    inventory = comm.downloads.query_station_inventory(start_date, end_date)
    inventory = inventory.select(channel="*Z")
    number_of_stations = comm.visualizations.plot_station_inventory_in_query(inventory)
    import matplotlib.pyplot as plt
    plt.title("%d stations available between %s and %s" %(number_of_stations, 
                                                          start_date.datetime.strftime("%Y-%m-%d"),
                                                          end_date.datetime.strftime("%Y-%m-%d")))
    
    comm.visualizations.plot_data_availability(inventory, start_date, end_date)

    plt.show()


@command_group("Plotting")
def lasif_plot_data_availability(parser, args):
    """
    Plot data availability for the station channels in the project 
    and indicates the event occurence times
    """
    args = parser.parse_args(args) 
    comm = _find_project_comm(".", args.read_only_caches)
    channels = comm.stations.get_all_channels()
    from obspy.core.inventory import Inventory
    import obspy
    read_filenames = []
    inventory = Inventory()
    print("Reading inventory for all stations. This might take a while ...")
    for channel in channels:
        filename = channel["filename"]
        if filename in read_filenames:
            continue
        else:
            inv = obspy.read_inventory(filename)
            if inv:
                inventory.extend(inv)
    inventory = inventory.select(channel="*Z")
    
        
    import matplotlib.pyplot as plt
    import numpy as np
    from obspy.core import UTCDateTime
    from matplotlib.dates import date2num
    print("Reading event catalog ...")
    events = list(comm.events.get_all_events().values())
    dates = []
    for event in events:
        dates.append(event["origin_time"].datetime)
    start_date = np.min(dates)
    start_date = UTCDateTime(start_date.year,start_date.month,start_date.day)
    end_date = np.max(dates)
    end_date = UTCDateTime(end_date.year,end_date.month,end_date.day)
    ax0, ax1 = comm.visualizations.plot_data_availability(inventory, start_date, end_date)
    
    for date in dates:
        ax0.plot([date2num(date), date2num(date)], [ax0.get_ylim()[0], ax0.get_ylim()[1]], 'k')
        ax1.plot([date2num(date), date2num(date)], [ax1.get_ylim()[0], ax1.get_ylim()[1]], 'k')
        
    plt.show()
   
    
    
@command_group("Data Acquisition")
def lasif_add_gcmt_events(parser, args):
    """
    Selects and adds optimally distributed events from the GCMT catalog.
    """
    parser.add_argument("count", type=int,
                        help="maximum amount of events to add")
    parser.add_argument("min_year", default=None, type=int,
                        help="minimum year from which to add events")
    parser.add_argument("max_year", default=None, type=int,
                        help="maximum year from which to add events")
    parser.add_argument("--select", default=1, type=int,
                        help="Flag for statistical selection of the events based on their spatial and time distibution")
    parser.add_argument("--min_magnitude", default=None, type=int,
                        help="minimum magnitude from which to add events")
    parser.add_argument("--max_magnitude", default=None, type=int,
                        help="maximum magnitude from which to add events")
    
    
    args = parser.parse_args(args)
    from lasif.tools.query_gcmt_catalog import add_new_events
    comm = _find_project_comm(".", args.read_only_caches)

    add_new_events(comm=comm, count=args.count,
                   min_year=args.min_year, max_year=args.max_year,
                   statistical_selection=args.select,
                   min_magnitude=args.min_magnitude, max_magnitude=args.max_magnitude)



@command_group("Project Management")
def lasif_info(parser, args):
    """
    Print a summary of the project.
    """
    args = parser.parse_args(args)

    comm = _find_project_comm(".", args.read_only_caches)
    print((comm.project))


@mpi_enabled
@command_group("Data Acquisition")
def lasif_download_data(parser, args):
    """
    Download waveform and station data for one event.
    """
    parser.add_argument(
        "events", help="One or more events. If none given, all will be done.",
        nargs="*")
    parser.add_argument("--providers", default=None,
                        type=str,
                        help="FDSN providers to query. Will use all known "
                             "ones if not set.")
    parser.add_argument("--networks", default=None,
                        type=str, 
                        help="seismic networks (comma separated) to download in the domain, eg. ``IU,G`` ")

    args = parser.parse_args(args)
    events = args.events if args.events else None
    providers = args.providers
    if providers is not None:
        providers = providers.split(',')
        Providers=[]
        for prov in providers:
            Providers.append(prov)
    else:
        Providers = None
    Networks = args.networks
    
    comm = _find_project_comm(".", args.read_only_caches)
    
    # No need to perform these checks on all ranks.
    exceptions = []
    if MPI.COMM_WORLD.rank == 0:

        # Check if the event ids are valid.
        if  events:
            for event_name in events:

                if not comm.events.has_event(event_name):
                    msg = "Event '%s' not found." % event_name
                    exceptions.append(msg)
                    break

    # Raise any exceptions on all ranks if necessary.
    exceptions = MPI.COMM_WORLD.bcast(exceptions, root=0)
    if exceptions:
        raise LASIFCommandLineException(exceptions[0])
    
    comm.downloads.download_data(events, providers=Providers, networks=Networks)


@command_group("Event Management")
def lasif_list_events(parser, args):
    """
    Print a list of all events in the project.
    """
    parser.add_argument("--details", help="print details of filecounts for "
                                          "all events",
                        action="store_true")
    parser.add_argument("--list", help="Show only a list of events. Good for "
                                       "scripting bash.",
                        action="store_true")
    args = parser.parse_args(args)

    from lasif.tools.prettytable import PrettyTable
    comm = _find_project_comm(".", args.read_only_caches)

    if args.list and args.details:
        raise LASIFCommandLineException("--list and --details cannot both "
                                        "be specified.")

    if args.list is False:
        print(("%i event%s in project:" % (comm.events.count(),
                                          "s" if comm.events.count() != 1 else "")))

    if args.details is True:
        tab = PrettyTable(["Event Name", "Lat/Lng/Depth(km)/Mag",
                           "# raw/preproc/synth"])
        tab.align["Event Name"] = "l"
        for event in comm.events.list():
            ev = comm.events.get(event)
            count = comm.project.get_filecounts_for_event(event)
            tab.add_row([
                event, "%6.1f / %6.1f / %3i / %3.1f" % (
                    ev["latitude"], ev["longitude"], int(ev["depth_in_km"]),
                    ev["magnitude"]),
                "%4i / %5i / %4i" % (
                    count["raw_waveform_file_count"],
                    count["preprocessed_waveform_file_count"],
                    count["synthetic_waveform_file_count"])])
        print(tab)
    elif args.list is True:
        for event in sorted(comm.events.list()):
            print(event)
    else:
        tab = PrettyTable(["Event Name", "Lat/Lng/Depth(km)/Mag"])
        tab.align["Event Name"] = "l"
        for event in comm.events.list():
            ev = comm.events.get(event)
            tab.add_row([
                event, "%6.1f / %6.1f / %3i / %3.1f" % (
                    ev["latitude"], ev["longitude"], int(ev["depth_in_km"]),
                    ev["magnitude"])])
        print(tab)


@command_group("Project Management")
def lasif_build_all_caches(parser, args):
    """
    Build all caches to speed up subsequent operations.

    This is optional and might take a while. Otherwise the caches are built
    on demand which works fine but might impede on some workflows.
    """
    parser.add_argument("--quick", help="Only check caches for folders that "
                                        "do not have a cache.",
                        action="store_true")
    args = parser.parse_args(args)

    comm = _find_project_comm(".", read_only_caches=False)
    comm.project.build_all_caches(quick=args.quick)


@command_group("Project Management")
def lasif_list_models(parser, args):
    """
    Print a list of all models in the project.
    """
    args = parser.parse_args(args)

    comm = _find_project_comm(".", args.read_only_caches)
    models = comm.models.list()
    print(("%i model%s in project:" % (len(models), "s" if len(models) != 1
                                      else "")))
    for model in models:
        print(("\t%s" % model))


@command_group("Project Management")
def lasif_list_kernels(parser, args):
    """
    Print a list of all kernels in this project.
    """
    args = parser.parse_args(args)

    comm = _find_project_comm(".", args.read_only_caches)
    kernels = comm.kernels.list()
    print(("%i kernel%s in project:" % (
        len(kernels), "s" if len(kernels) != 1 else "")))
    for kernel in kernels:
        print(("\tIteration %3s and Event %s" % (kernel["iteration"],
                                                kernel["event"])))


@command_group("Plotting")
def lasif_plot_wavefield(parser, args):
    """
    Plots a SES3D wavefield.
    """
    parser.add_argument("iteration_name", help="name_of_the_iteration")
    parser.add_argument("event_name", help="name of the event")
    args = parser.parse_args(args)

    from lasif import ses3d_models

    comm = _find_project_comm(".", args.read_only_caches)

    event_name = comm.events.get(args.event_name)["event_name"]
    iteration_name = comm.iterations.get(args.iteration_name).long_name

    wavefield_dir = os.path.join(comm.project.paths["wavefields"], event_name,
                                 iteration_name)

    if not os.path.exists(wavefield_dir) or not os.listdir(wavefield_dir):
        msg = "No data available for event and iteration combination."
        raise LASIFCommandLineException(msg)

    handler = ses3d_models.RawSES3DModelHandler(
        wavefield_dir, model_type="wavefield")
    handler.rotation_axis = comm.project.domain["rotation_axis"]
    handler.rotation_angle_in_degree = comm.project.domain["rotation_angle"]

    while True:
        print(handler)
        print("")

        inp = input("Enter 'COMPONENT DEPTH' "
                        "('quit/exit' to exit): ").strip()
        if inp.lower() in ["quit", "q", "exit", "leave"]:
            break
        try:
            component, timestep, depth = inp.split()
        except:
            continue

        component = component = "%s %s" % (component, timestep)

        try:
            handler.parse_component(component)
        except:
            continue
        handler.plot_depth_slice(component, float(depth))


@command_group("Event Management")
def lasif_event_info(parser, args):
    """
    Print information about a single event.
    """
    parser.add_argument("event_name", help="name of the event")
    parser.add_argument("-v", help="Verbose. Print all contained events.",
                        action="store_true")
    args = parser.parse_args(args)
    event_name = args.event_name
    verbose = args.v

    comm = _find_project_comm(".", args.read_only_caches)
    if not comm.events.has_event(event_name):
        msg = "Event '%s' not found in project." % event_name
        raise LASIFCommandLineException(msg)

    event_dict = comm.events.get(event_name)

    print(("Earthquake with %.1f %s at %s" % (
          event_dict["magnitude"], event_dict["magnitude_type"],
          event_dict["region"])))
    print(("\tLatitude: %.3f, Longitude: %.3f, Depth: %.1f km" % (
          event_dict["latitude"], event_dict["longitude"],
          event_dict["depth_in_km"])))
    print(("\t%s UTC" % str(event_dict["origin_time"])))

    try:
        stations = comm.query.get_all_stations_for_event(event_name)
    except LASIFError:
        stations = {}

    if verbose:
        from lasif.utils import table_printer
        print(("\nStation and waveform information available at %i "
              "stations:\n" % len(stations)))
        header = ["id", "latitude", "longitude", "elevation_in_m",
                  "local depth"]
        keys = sorted(stations.keys())
        data = [[
            key, stations[key]["latitude"], stations[key]["longitude"],
            stations[key]["elevation_in_m"], stations[key]["local_depth_in_m"]]
            for key in keys]
        table_printer(header, data)
    else:
        print(("\nStation and waveform information available at %i stations. "
              "Use '-v' to print them." % len(stations)))


@command_group("Plotting")
def lasif_plot_stf(parser, args):
    """
    Plot the source time function for one iteration.
    """
    import lasif.visualization

    parser.add_argument("iteration_name", help="name of the iteration")
    args = parser.parse_args(args)
    iteration_name = args.iteration_name

    comm = _find_project_comm(".", args.read_only_caches)

    iteration = comm.iterations.get(iteration_name)
    pp = iteration.get_process_params()
    freqmin = pp["highpass"]
    freqmax = pp["lowpass"]

    stf = iteration.get_source_time_function()

    # Ignore lots of potential warnings with some plotting functionality.
    lasif.visualization.plot_tf(stf["data"], stf["delta"], freqmin=freqmin,
                                freqmax=freqmax)


@command_group("Iteration Management")
def lasif_generate_all_input_files(parser, args):
    """
    Generates all input files for a certain iteration.

    TYPE denotes the type of simulation to run. Available types are
        * "normal_simulation"
        * "adjoint_forward"
        * "adjoint_reverse"
    """
    parser.add_argument("iteration_name", help="name of the iteration")
    parser.add_argument("--simulation_type",
                        choices=("normal_simulation", "adjoint_forward",
                                 "adjoint_reverse"),
                        default="normal_simulation",
                        help="type of simulation to run")
    args = parser.parse_args(args)
    iteration_name = args.iteration_name
    simulation_type = args.simulation_type

    comm = _find_project_comm(".", args.read_only_caches)
    simulation_type = simulation_type.replace("_", " ")

    it = comm.iterations.get(iteration_name)
    events = sorted(it.events.keys())
    for _i, event in enumerate(events):
        print(("Generating input files for event %i of %i..." % (_i + 1,
                                                                len(events))))
        comm.actions.generate_input_files(iteration_name, event,
                                          simulation_type)


@command_group("Iteration Management")
def lasif_generate_input_files(parser, args):
    """
    Generate the input files for the waveform solver.

    TYPE denotes the type of simulation to run. Available types are
        * "normal_simulation"
        * "adjoint_forward"
        * "adjoint_reverse"
    """
    parser.add_argument("iteration_name", help="name of the iteration")
    parser.add_argument("event_name", help="name of the event")
    parser.add_argument("--simulation_type",
                        choices=("normal_simulation", "adjoint_forward",
                                 "adjoint_reverse"),
                        default="normal_simulation",
                        help="type of simulation to run")
    args = parser.parse_args(args)
    iteration_name = args.iteration_name
    event_name = args.event_name
    simulation_type = args.simulation_type

    comm = _find_project_comm(".", args.read_only_caches)
    simulation_type = simulation_type.replace("_", " ")
    comm.actions.generate_input_files(iteration_name, event_name,
                                      simulation_type)


@command_group("Project Management")
def lasif_init_project(parser, args):
    """
    Create a new project.
    """
    parser.add_argument("folder_path", help="where to create the project")
    args = parser.parse_args(args)
    folder_path = args.folder_path

    if os.path.exists(folder_path):
        msg = "The given FOLDER_PATH already exists. It must not exist yet."
        raise LASIFCommandLineException(msg)
    folder_path = os.path.abspath(folder_path)
    try:
        os.makedirs(folder_path)
    except:
        msg = "Failed creating directory %s. Permissions?" % folder_path
        raise LASIFCommandLineException(msg)

    Project(project_root_path=folder_path,
            init_project=os.path.basename(folder_path))

    print(("Initialized project in: \n\t%s" % folder_path))



@command_group("Iteration Management")
def lasif_finalize_adjoint_sources(parser, args):
    """
    Finalize the adjoint sources.
    """
    parser.add_argument("iteration_name", help="name of the iteration")
    parser.add_argument("event_name", help="name of the event")
    args = parser.parse_args(args)
    iteration_name = args.iteration_name
    event_name = args.event_name

    comm = _find_project_comm(".", args.read_only_caches)
    comm.actions.finalize_adjoint_sources(iteration_name, event_name)


@command_group("Iteration Management")
def lasif_calculate_all_adjoint_sources(parser, args):
    """
    Calculates all adjoint sources for a given iteration and event.
    """
    parser.add_argument("iteration_name", help="name of the iteration")
    parser.add_argument("event_name", help="name of the event")
    args = parser.parse_args(args)
    iteration_name = args.iteration_name
    event_name = args.event_name

    comm = _find_project_comm(".", args.read_only_caches)
    comm.actions.calculate_all_adjoint_sources(iteration_name, event_name)


@mpi_enabled
@command_group("Iteration Management")
def lasif_select_windows(parser, args):
    """
    Autoselect windows for a given event and iteration combination.

    This function works with MPI. Don't use too many cores, I/O quickly
    becomes the limiting factor. It also works without MPI but then only one
    core actually does any work.
    """
    parser.add_argument("iteration_name", help="name of the iteration")
    parser.add_argument("event_name", help="name of the event")
    args = parser.parse_args(args)

    iteration = args.iteration_name
    event = args.event_name

    comm = _find_project_comm_mpi(".", args.read_only_caches)

    comm.actions.select_windows(event, iteration)


@mpi_enabled
@command_group("Iteration Management")
def lasif_select_all_windows(parser, args):
    """
    Autoselect all windows for a given iteration.

    This function works with MPI. Don't use too many cores, I/O quickly
    becomes the limiting factor. It also works without MPI but then only one
    core actually does any work.
    """
    parser.add_argument("iteration_name", help="name of the iteration")
    args = parser.parse_args(args)

    iteration = args.iteration_name

    comm = _find_project_comm_mpi(".", args.read_only_caches)

    events = comm.events.list()

    for _i, event in enumerate(events):
        if MPI.COMM_WORLD.rank == 0:
            print(("\n{green}"
                  "==========================================================="
                  "{reset}".format(green=colorama.Fore.GREEN,
                                   reset=colorama.Style.RESET_ALL)))
            print(("Starting window selection for event %i of %i..." % (
                  _i + 1, len(events))))
            print(("{green}"
                  "==========================================================="
                  "{reset}\n".format(green=colorama.Fore.GREEN,
                                     reset=colorama.Style.RESET_ALL)))
        MPI.COMM_WORLD.barrier()
        comm.actions.select_windows(event, iteration)


@command_group("Iteration Management")
def lasif_launch_misfit_gui(parser, args):
    """
    Launch the misfit GUI.
    """
    args = parser.parse_args(args)

    comm = _find_project_comm(".", args.read_only_caches)

    from lasif.misfit_gui.misfit_gui import launch
    launch(comm)


@command_group("Plotting")
def lasif_launch_model_gui(parser, args):
    """
    Launch the model GUI.
    """
    args = parser.parse_args(args)

    comm = _find_project_comm(".", args.read_only_caches)

    from lasif.ses3d_model_gui.model_gui import launch
    launch(comm)


@command_group("Plotting")
def lasif_plot_model(parser, args):
    """
    Directly plot a model to a file.

    Useful for scripting LASIF.
    """
    parser.add_argument("model_name", help="name of the model")
    parser.add_argument("depth", type=float, help="the depth at which to plot")
    parser.add_argument("component", type=str,
                        help="the component to plot")
    parser.add_argument("filename", type=str,
                        help="Output filename. Use '-' to not write to a "
                             "file but directly show the kernel.")

    args = parser.parse_args(args)

    comm = _find_project_comm(".", args.read_only_caches)

    if args.model_name not in comm.models.list():
        raise LASIFCommandLineException("Model '%s' not known to LASIF." %
                                        args.model_name)

    import matplotlib.pyplot as plt

    plt.figure(figsize=(15, 15))

    model = comm.models.get_model_handler(args.model_name)
    model.parse_component(args.component)

    m = comm.project.domain.plot()
    im = model.plot_depth_slice(component=args.component,
                                depth_in_km=args.depth, m=m)["mesh"]

    # make a colorbar and title
    m.colorbar(im, "right", size="3%", pad='2%')
    plt.title(str(args.depth) + ' km')

    if args.filename == "-":
        plt.show()
    else:
        plt.savefig(args.filename, dpi=100)
    plt.close()


@command_group("Plotting")
def lasif_plot_kernel(parser, args):
    """
    Directly plot a kernel to a file.

    Useful for scripting LASIF. As they are not often in the LASIF project
    this is one of the view commands that will work on data outside of LASIF.
    """
    parser.add_argument("folder", help="The folder containing the gradients.")
    parser.add_argument("depth", type=float,
                        help="The depth at which to plot.")
    parser.add_argument("component", type=str,
                        help="The component to plot.")
    parser.add_argument("filename", type=str,
                        help="Output filename. Use '-' to not write to a "
                             "file but directly show the kernel.")

    args = parser.parse_args(args)

    comm = _find_project_comm(".", args.read_only_caches)

    import matplotlib.pyplot as plt
    from lasif.ses3d_models import RawSES3DModelHandler

    plt.figure(figsize=(15, 15))

    model = RawSES3DModelHandler(
        directory=args.folder, domain=comm.project.domain,
        model_type="kernel")

    model.parse_component(args.component)

    m = comm.project.domain.plot()
    im = model.plot_depth_slice(component=args.component,
                                depth_in_km=args.depth, m=m)["mesh"]

    # make a colorbar and title
    m.colorbar(im, "right", size="3%", pad='2%')
    plt.title(str(args.depth) + ' km')

    if args.filename == "-":
        plt.show()
    else:
        plt.savefig(args.filename, dpi=100)
    plt.close()


@command_group("Iteration Management")
def lasif_create_new_iteration(parser, args):
    """
    Create a new iteration.
    """
    parser.add_argument("iteration_name", help="name of the iteration")
    parser.add_argument("min_period", type=float,
                        help="the minimum period of the iteration")
    parser.add_argument("max_period", type=float,
                        help="the maximum period of the iteration")
    parser.add_argument("solver_name", help="name of the solver",
                        choices=("SES3D_4_1", "SES3D_2_0",
                                 "SPECFEM3D_CARTESIAN",
                                 "SPECFEM3D_GLOBE_CEM"))
    parser.add_argument("--seconds_prior", type=float, default=5.,
                        help="nb of seconds prior the theoretical phase arrival time used to window seismograms for quality control, default 5")
    parser.add_argument("--window_length", type=float, default=50.,
                        help="Time window length in seconds used to window the phase of interest in the seismograms, used for suqlity control, default 50")

    args = parser.parse_args(args)
    iteration_name = args.iteration_name
    solver_name = args.solver_name
    min_period = args.min_period
    max_period = args.max_period
    seconds_prior_arrival = args.seconds_prior
    window_length_in_sec = args.window_length
    if min_period >= max_period:
        msg = "min_period needs to be smaller than max_period."
        raise LASIFCommandLineException(msg)
    if seconds_prior_arrival >= window_length_in_sec:
        msg = "seconds_prior needs to be smaller than window_length."
        raise LASIFCommandLineException(msg)

    comm = _find_project_comm(".", args.read_only_caches)
    #print(comm.query.get_stations_for_all_events())
    comm.iterations.create_new_iteration(
        iteration_name=iteration_name,
        solver_name=solver_name,
        events_dict=comm.query.get_stations_for_all_events(),
        min_period=min_period,
        max_period=max_period,
        seconds_prior_arrival=seconds_prior_arrival,
        window_length_in_sec=window_length_in_sec)


@command_group("Iteration Management")
def lasif_update_iteration(parser, args):
    """
    Update the iteration.
    """
    parser.add_argument("iteration_name", help="name of the current iteration to update")
    parser.add_argument("new_iteration_name", help="name of the new iteration")

    args = parser.parse_args(args)
    iteration_name = args.iteration_name
    new_iteration_name = args.new_iteration_name


    if new_iteration_name == iteration_name:
        msg = "The new iteration name should be different than the current iteration names."
        raise LASIFCommandLineException(msg)


    comm = _find_project_comm(".", args.read_only_caches)
    events_dict = comm.query.get_stations_for_all_processed_events(iteration_name)

    comm.iterations.update_iteration(iteration_name,
                                     new_iteration_name, events_dict)


@command_group("Iteration Management")
def lasif_create_successive_iteration(parser, args):
    """
    Create an iteration based on an existing one.

    It will take all settings in one iteration and transfers them to another
    iteration. Any comments will be deleted.
    """
    parser.add_argument("existing_iteration",
                        help="name of the existing iteration")
    parser.add_argument("new_iteration", help="name of the new iteration")
    args = parser.parse_args(args)
    existing_iteration_name = args.existing_iteration
    new_iteration_name = args.new_iteration

    comm = _find_project_comm(".", args.read_only_caches)

    comm.iterations.create_successive_iteration(
        existing_iteration_name=existing_iteration_name,
        new_iteration_name=new_iteration_name)


@mpi_enabled
@command_group("Iteration Management")
def lasif_compare_misfits(parser, args):
    """
    Compares the misfit between two iterations. Will only consider windows
    that are identical in both iterations as the comparision is otherwise
    meaningless.
    """
    from lasif import LASIFAdjointSourceCalculationError

    parser.add_argument("from_iteration",
                        help="past iteration")
    parser.add_argument("to_iteration", help="current iteration")
    args = parser.parse_args(args)

    comm = _find_project_comm_mpi(".", args.read_only_caches)

    _starting_time = time.time()

    # Read on each core as pickling/broadcasting them prooves to be
    # difficult. Should be possible if this ever becomes a performance issue.
    from_it = comm.iterations.get(args.from_iteration)
    to_it = comm.iterations.get(args.to_iteration)

    if MPI.COMM_WORLD.rank == 0:
        # Get a list of events that are in both, the new and the old iteration.
        events = sorted(set(from_it.events.keys()).intersection(
            set(to_it.events.keys())))
        event_count = len(events)

        # Split into a number of events per MPI process.
        events = split(events, MPI.COMM_WORLD.size)

        print(" => Calculating misfit change from iteration '%s' to " \
            "iteration '%s' ..." % (from_it.name, to_it.name))
        print(" => Launching calculations on %i core(s)\n" % \
            MPI.COMM_WORLD.size)

    else:
        events = None

    # Scatter jobs
    events = MPI.COMM_WORLD.scatter(events, root=0)

    total_misfit_from = 0
    total_misfit_to = 0

    all_events = collections.defaultdict(list)

    # Loop over each event.
    for _i, event in enumerate(events):
        # Get the windows from both.
        window_group_to = comm.windows.get(event, to_it)
        window_group_from = comm.windows.get(event, from_it)

        event_weight = from_it.events[event]["event_weight"]

        # Get a list of channels shared amongst both.
        shared_channels = set(window_group_to.list()).intersection(
            set(window_group_from.list()))

        # On rank 0, show a progressbar because it can take forever.
        if MPI.COMM_WORLD.rank == 0:
            widgets = [
                "Approximately event %i of %i: " % (
                    _i * MPI.COMM_WORLD.size + 1, event_count),
                progressbar.Percentage(),
                progressbar.Bar(), "", progressbar.ETA()]
            pbar = progressbar.ProgressBar(
                widgets=widgets, maxval=len(shared_channels)).start()

        # Loop over each channel.
        for _i, channel in enumerate(shared_channels):
            if MPI.COMM_WORLD.rank == 0:
                pbar.update(_i)
            window_collection_from = window_group_from.get(channel)
            window_collection_to = window_group_to.get(channel)

            station_weight = from_it.events[event]["stations"][
                ".".join(channel.split(".")[:2])]["station_weight"]

            channel_misfit_from = 0
            channel_misfit_to = 0
            total_channel_weight = 0

            # Loop over each window in that channel.
            for win_from in window_collection_from.windows:
                try:
                    idx = window_collection_to.windows.index(win_from)
                    win_to = window_collection_to.windows[idx]
                except ValueError:
                    continue

                try:
                    misfit_from = win_from.misfit_value
                except LASIFAdjointSourceCalculationError:
                    continue
                except LASIFNotFoundError as e:
                    print(str(e))
                    continue

                try:
                    misfit_to = win_to.misfit_value
                except Exception as e:
                    print(e)
                    # Random penalty...but how else to compare?
                    misfit_to = 2.0 * misfit_from

                channel_misfit_from += misfit_from * win_from.weight
                channel_misfit_to += misfit_to * win_from.weight
                total_channel_weight += win_from.weight

            # Rare - but sometimes all windows for a certain channel fail
            # the calculation.
            if total_channel_weight == 0:
                continue

            # Make sure the misfits are consistent with the adjoint source
            # calculations!
            channel_misfit_from *= \
                event_weight * station_weight / total_channel_weight
            channel_misfit_to *= \
                event_weight * station_weight / total_channel_weight

            total_misfit_from += channel_misfit_from
            total_misfit_to += channel_misfit_to

            if (misfit_to - misfit_from) < -1.5:
                print((event, channel, misfit_from - misfit_to))
            all_events[event].append(misfit_to - misfit_from)
        if MPI.COMM_WORLD.rank == 0:
            pbar.finish()

    _all_events = MPI.COMM_WORLD.gather(all_events, root=0)

    total_misfit_from = MPI.COMM_WORLD.reduce(total_misfit_from, root=0)
    total_misfit_to = MPI.COMM_WORLD.reduce(total_misfit_to, root=0)

    # Only rank 0 continues.
    if MPI.COMM_WORLD.rank != 0:
        return

    # Collect in singular dictionary again.
    all_events = {}
    [all_events.update(_i) for _i in _all_events]

    if not all_events:
        raise LASIFCommandLineException("No misfit values could be compared.")

    print("\nTotal misfit in Iteration %s: %g" % (from_it.name,
                                                  total_misfit_from))
    print("Total misfit in Iteration %s: %g" % (to_it.name,
                                                total_misfit_to))

    _ending_time = time.time()

    print("\n => Computation time: %.1f seconds" % (_ending_time -
                                                    _starting_time))

    import matplotlib.pylab as plt
    import numpy as np

    plt.figure(figsize=(20, 3 * len(all_events)))
    plt.suptitle("Misfit change of measurements going from iteration"
                 " '%s' to iteration '%s'" % (from_it.name, to_it.name))
    for i, event_name in enumerate(sorted(all_events.keys())):
        values = np.array(all_events[event_name])
        colors = np.array(["green"] * len(values))
        colors[values > 0] = "red"
        plt.subplot(len(all_events), 1, i + 1)
        plt.bar(np.arange(len(values)), values, color=colors)
        plt.ylabel("difference")
        plt.xlim(0, len(values) - 1)
        plt.xticks([])
        plt.title("%i measurements with identical windows for event '%s'" %
                  (len(values), event_name))

    output_folder = comm.project.get_output_folder(
        type="misfit_comparisons", tag="misfit_comparision")
    filename = os.path.join(output_folder, "misfit_comparision.pdf")
    plt.savefig(filename)
    print("\nSaved figure to '%s'" % os.path.relpath(filename))


@command_group("Iteration Management")
def lasif_migrate_windows(parser, args):
    """
    Migrates windows from one iteration to the next.
    """
    parser.add_argument("from_iteration",
                        help="iteration containing windows")
    parser.add_argument("to_iteration", help="iteration windows will "
                                             "be migrated to")
    args = parser.parse_args(args)
    comm = _find_project_comm(".", args.read_only_caches)

    from_it = comm.iterations.get(args.from_iteration)
    to_it = comm.iterations.get(args.to_iteration)

    print("Migrating windows from iteration '%s' to iteration '%s'..." % (
        from_it.name, to_it.name))

    for event_name, stations in list(to_it.events.items()):
        stations = list(stations["stations"].keys())

        window_group_from = comm.windows.get(event_name, from_it.name)
        window_group_to = comm.windows.get(event_name, to_it.name)
        contents_to = set(window_group_to.list())
        contents_from = set(window_group_from.list())
        contents = contents_from - contents_to

        # Remove all not part of this iterations station.
        filtered_contents = filter(
            lambda x: ".".join(x.split(".")[:2]) in stations,
            contents)

        for channel_id in filtered_contents:
            coll = window_group_from.get(channel_id)
            coll.synthetics_tag = to_it.name
            f = window_group_to._get_window_filename(channel_id)
            coll.filename = f
            coll.write()


@command_group("Iteration Management")
def lasif_list_iterations(parser, args):
    """
    Print a list of all iterations in the project.
    """
    args = parser.parse_args(args)

    comm = _find_project_comm(".", args.read_only_caches)

    it_len = comm.iterations.count()

    print(("%i iteration%s in project:" % (it_len,
                                          "s" if it_len != 1 else "")))
    for iteration in comm.iterations.list():
        print(("\t%s" % iteration))


@command_group("Iteration Management")
def lasif_iteration_info(parser, args):
    """
    Print information about a single iteration.
    """
    parser.add_argument("iteration_name", help="name of the iteration")
    args = parser.parse_args(args)
    iteration_name = args.iteration_name

    comm = _find_project_comm(".", args.read_only_caches)
    if not comm.iterations.has_iteration(iteration_name):
        msg = ("Iteration '%s' not found. Use 'lasif list_iterations' to get "
               "a list of all available iterations.") % iteration_name
        raise LASIFCommandLineException(msg)

    print((comm.iterations.get(iteration_name)))


@command_group("Project Management")
def lasif_remove_empty_coordinate_entries(parser, args):
    """
    Remove all empty coordinate entries in the inventory cache.

    This is useful if you want to try to download coordinates again.
    """
    args = parser.parse_args(args)

    comm = _find_project_comm(".", args.read_only_caches)
    comm.inventory_db.remove_coordinate_less_stations()

    print("SUCCESS")


@mpi_enabled
@command_group("Iteration Management")
def lasif_preprocess_data(parser, args):
    """
    Launch data preprocessing.

    This function works with MPI. Don't use too many cores, I/O quickly
    becomes the limiting factor. It also works without MPI but then only one
    core actually does any work.
    """
    parser.add_argument("iteration_name", help="name of the iteration")
    parser.add_argument(
        "events", help="One or more events. If none given, all will be done.",
        nargs="*")
    parser.add_argument("--snr",
                        help="Relative noise level threshold above which data trace will be disregarded, 0.5 by default")
    parser.add_argument("--components", default="ENZ", help="list of components to process, examples: ENZ, or Z or RTZ")
    parser.add_argument("--svd_selection",default="False", choices=["True", "False"],
                        help="``True``: for teleseismic configuration and waveform selection based on their similarity, "
                        "``False``: preferred for regional waveforms, no selection applied ")
    parser.add_argument("--recompute_files",default="False", choices=["True", "False"],
                        help="``False``: will not process already existed preprocessed files, "
                        "``True``: will recompute already existing preprocessed files, "
                        "          it might be useful if you want to change the snr threshold ")

    args = parser.parse_args(args)
    iteration_name = args.iteration_name
    events = args.events if args.events else None
    svd_selection = args.svd_selection
    recompute_files = args.recompute_files
    noise_threshold = float(args.snr) if args.snr else None
    comps = args.components
    components=[]
    for comp in comps:
        components.append(comp)    

    comm = _find_project_comm_mpi(".", args.read_only_caches)

    # No need to perform these checks on all ranks.
    exceptions = []
    if MPI.COMM_WORLD.rank == 0:
        if not comm.iterations.has_iteration(iteration_name):
            msg = ("Iteration '%s' not found. Use 'lasif list_iterations' to "
                   "get a list of all available iterations.") % iteration_name
            exceptions.append(msg)


        # Check if the event ids are valid.
        if not exceptions and events:
            for event_name in events:

                if not comm.events.has_event(event_name):
                    msg = "Event '%s' not found." % event_name
                    exceptions.append(msg)
                    break

    # Raise any exceptions on all ranks if necessary.
    exceptions = MPI.COMM_WORLD.bcast(exceptions, root=0)
    if exceptions:
        raise LASIFCommandLineException(exceptions[0])


    comm.actions.preprocess_data(iteration_name, components, svd_selection, noise_threshold, events, recompute_files)


@mpi_enabled
@command_group("Iteration Management")
def lasif_deconvolve_stf(parser, args):
    """
    Launch source wavelet estimation.

    This function works with MPI. Don't use too many cores, I/O quickly
    becomes the limiting factor. It also works without MPI but then only one
    core actually does any work.
    """
    parser.add_argument("iteration_name", help="name of the iteration")
    parser.add_argument(
        "events", help="One or more events. If none given, all will be done.",
        nargs="*")
    parser.add_argument("--components", default="ENZ", help="list of components to process, examples: ENZ, or Z or RTZ")

    args = parser.parse_args(args)
    iteration_name = args.iteration_name
    events = args.events if args.events else None
    comps = args.components
    components=[]
    for comp in comps:
        components.append(comp)


    comm = _find_project_comm_mpi(".", args.read_only_caches)

    # No need to perform these checks on all ranks.
    exceptions = []
    if MPI.COMM_WORLD.rank == 0:
        if not comm.iterations.has_iteration(iteration_name):
            msg = ("Iteration '%s' not found. Use 'lasif list_iterations' to "
                   "get a list of all available iterations.") % iteration_name
            exceptions.append(msg)


        # Check if the event ids are valid.
        if not exceptions and events:
            for event_name in events:

                if not comm.events.has_event(event_name):
                    msg = "Event '%s' not found." % event_name
                    exceptions.append(msg)
                    break

    # Raise any exceptions on all ranks if necessary.
    exceptions = MPI.COMM_WORLD.bcast(exceptions, root=0)
    if exceptions:
        raise LASIFCommandLineException(exceptions[0])


    comm.actions.stf_estimate(iteration_name, components, events)

    """
    import obspy
    import numpy as np
    stf = comm.waveforms.get_waveform_stf(event_name, iteration_name, component="Z")
    starttime = stf[0].stats.starttime
    endtime = stf[0].stats.endtime

    st = obspy.Stream()
    offset = np.zeros(len(stations_info),dtype=float)
    tt_arrival = np.zeros(len(stations_info),dtype=float)
    phase_angle= np.zeros(len(stations_info),dtype=float)
    lngs = np.zeros(len(stations_info),dtype=float)
    lats = np.zeros(len(stations_info),dtype=float)
    values = np.zeros(len(stations_info),dtype=float)
    amp_values = values.copy()
    cc = values.copy()
    cc_shift = values.copy()
    labels = []
    for i, station in enumerate(stations_info):
        lngs[i] = stations_info[station]["longitude"]
        lats[i] = stations_info[station]["latitude"]
        values[i] = stations_info[station]["time shift"]
        offset[i] = stations_info[station]["epicentral_distance"]
        tt_arrival[i] = stations_info[station]["tt_arrival"]
        phase_angle[i] = stations_info[station]["phase_angle"]
        amp_values[i] = stations_info[station]["amp anomaly"]
        cc[i] = stations_info[station]["cc"]
        cc_shift[i] = stations_info[station]["cc_shift"]
        labels.append(station)

        st += obspy.read(stations_info[station]["input_file"]).trim(starttime,endtime)


    import numpy as np
    labels=np.array(labels)
    #values -= np.mean(values) #we already get rid of the mean time as we compute it on the averaged synthetics with stf
    # check stations where time shift > a threshold
    # then ritere avec la mediane des times-shitfs definis ds une certaine region autour (rayon de 30km) to check if it is really an outlier dans une region precise

    time_residual_threshold = 1.
    mask = np.ones(len(values), bool)
    indices = [i for i, x in enumerate(np.abs(values)) if x >= time_residual_threshold]
    if indices:
        mask[indices] = False
    vmax = np.max(np.abs(values[mask]))
    print(vmax)
    #vmax = 2*np.median(np.abs(values))

    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(16,7.5))
    ax1 = fig.add_subplot(231)
    map_object = comm.visualizations.plot_stations_for_event(events[0], iteration_name, 0, color = [], ax=ax1)
    x, y = map_object(lngs, lats)
    if indices:
        station_rmv_plot = map_object.scatter(x[indices], y[indices], c=values[indices], s=35, 
                                              marker="o", zorder=5, cmap='coolwarm_r',
                                              vmin = -vmax, vmax=vmax)
        station_rmv_plot._edgecolors = np.array([[0.0, 0.0, 0.0, 1.0]])
        station_rmv_plot._edgecolors_original = "black"
    '''
    sc1 = map_object.scatter(x[mask], y[mask], c=values[mask], s=35, 
                    marker="v", zorder=5, cmap='RdBu',
                    vmin = -vmax, vmax=vmax)
    '''
    sc1 = map_object.scatter(x, y, c=values, s=35, 
                             marker="v", zorder=5, cmap='RdBu',
                             vmin = -vmax, vmax=vmax)
    # Setting the picker overwrites the edgecolor attribute on certain
    # matplotlib and basemap versions. Fix it here.    
    sc1._edgecolors = np.array([[0.0, 0.0, 0.0, 1.0]])
    sc1._edgecolors_original = "black"
    cbar = plt.colorbar(sc1, ax=[ax1], extend = "both", location="top")
    cbar.ax.set_xlabel('Time residuals (s)')
    ax1.set_title('')
    ax1.axes.get_xaxis().set_visible(False)
    ax1.axes.get_yaxis().set_visible(False)
    annot1 = ax1.annotate("", xy=(0,0), xytext=(-10,10),textcoords="offset points",
                          bbox=dict(boxstyle="round", fc="w"),
                          arrowprops=dict(arrowstyle="->"))
    annot1.set_visible(False)    




    ax4 = fig.add_subplot(234)
    plt.hist(values, bins = len(values)/2)
    plt.axvline(x=np.median(np.abs(values)), ymin=0, ymax=10, color='k')
    ax4.set_xlabel('Time residuals (s)')
    ax4.set_ylabel('Count')


    ax2 = fig.add_subplot(232, sharex=ax1, sharey=ax1)
    #amp_values = np.log10(amp_values/np.mean(amp_values))
    #amp_values -= np.mean(amp_values)
    vmax = np.max(np.abs(amp_values))
    map_object = comm.visualizations.plot_stations_for_event(events[0], iteration_name, 0, color = [], ax=ax2)
    x, y = map_object(lngs, lats)
    sc2 = map_object.scatter(x, y, c=amp_values, s=35, 
                             marker="v", zorder=5) #,vmin = -vmax, vmax=vmax)
    # Setting the picker overwrites the edgecolor attribute on certain
    # matplotlib and basemap versions. Fix it here.    
    sc2._edgecolors = np.array([[0.0, 0.0, 0.0, 1.0]])
    sc2._edgecolors_original = "black"
    cbar = plt.colorbar(sc2, ax=[ax2], extend = "both", location="bottom")
    cbar.ax.set_xlabel('Amplitude residuals')
    ax2.set_title('')
    ax2.axes.get_xaxis().set_visible(False)
    ax2.axes.get_yaxis().set_visible(False)
    annot2 = ax2.annotate("", xy=(0,0), xytext=(-10,10),textcoords="offset points",
                          bbox=dict(boxstyle="round", fc="w"),
                          arrowprops=dict(arrowstyle="->"))
    annot2.set_visible(False) 

    ax5 = fig.add_subplot(235)
    plt.hist(amp_values, bins = len(values)/2)
    plt.axvline(x=np.median(np.abs(amp_values)), ymin=0, ymax=10, color='k')
    ax5.set_xlabel('Amplitude residuals')
    ax5.set_ylabel('Count')



    ax3 = fig.add_subplot(233, sharex=ax1, sharey=ax1)
    #amp_values = np.log10(amp_values/np.mean(amp_values))
    #amp_values -= np.mean(amp_values)
    vmax = np.max(np.abs(cc))
    map_object = comm.visualizations.plot_stations_for_event(events[0], iteration_name, 0, color = [], ax=ax3)
    x, y = map_object(lngs, lats)
    sc3 = map_object.scatter(x, y, c=cc, s=35, 
                             marker="v", zorder=5) #,vmin = -vmax, vmax=vmax)
    # Setting the picker overwrites the edgecolor attribute on certain
    # matplotlib and basemap versions. Fix it here.    
    sc3._edgecolors = np.array([[0.0, 0.0, 0.0, 1.0]])
    sc3._edgecolors_original = "black"
    cbar = plt.colorbar(sc3, ax=ax3, extend = "both")
    cbar.ax.set_ylabel('Correlation coefficient')
    ax3.set_title('')
    ax3.set_xticks([])
    ax3.set_yticks([])
    annot3 = ax3.annotate("", xy=(0,0), xytext=(-10,10),textcoords="offset points",
                          bbox=dict(boxstyle="round", fc="w"),
                          arrowprops=dict(arrowstyle="->"))
    annot3.set_visible(False) 

    ax6 = fig.add_subplot(236)
    plt.hist(cc, bins = len(values)/2)
    plt.axvline(x=np.median(np.abs(cc)), ymin=0, ymax=10, color='k')
    ax6.set_xlabel('Correlation coefficient')
    ax6.set_ylabel('Count')

    # see https://matplotlib.org/3.1.1/users/event_handling.html
    # https://stackoverflow.com/questions/7908636/possible-to-make-labels-appear-when-hovering-over-a-point-in-matplotlib




    def hover_stations(event):
        def update_anotation(event, ax, sc, annot, labels, cvalues):
            vis = annot.get_visible()
            cont, ind = sc.contains(event)
            if cont:
                pos = sc.get_offsets()[ind["ind"][0]]
                annot.xy = pos
                text = "%d,%s,%.2f"%(ind["ind"],
                                     [labels[n] for n in ind["ind"]][0],
                                     [cvalues[n] for n in ind["ind"]][0])
                annot.set_text(text)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

        if event.inaxes in [ax1, ax2, ax3]:
            if event.inaxes == ax1:
                update_anotation(event, ax1, sc1, annot1, labels, values)
            elif event.inaxes == ax2:
                update_anotation(event, ax2, sc2, annot2, labels, amp_values)
            elif event.inaxes == ax3:
                update_anotation(event, ax3, sc3, annot3, labels, cc)



    fig.canvas.mpl_connect("button_press_event", hover_stations) #motion_notify_event



    seconds_prior_arrival = 5
    window_length_in_sec = 50
    scaling = 0.5
    fig2 = plt.figure()
    ax = fig2.add_subplot(111)
    from lasif.visualization import plot_waveform_section
    plot_waveform_section(ax, st, offset, scale=scaling,colors='k')
    for arrival, dist, time_shift in zip(tt_arrival, offset, values):
        plt.plot(arrival-seconds_prior_arrival+time_shift, dist,'b.')
        plt.plot(arrival-seconds_prior_arrival+window_length_in_sec+time_shift, dist,'b.')
        plt.plot([arrival, arrival], 
                 [dist-0.25*scaling, dist+0.25*scaling],'k', linewidth = 1)
        #plt.plot([arrival-seconds_prior_arrival, arrival-seconds_prior_arrival], 
        #         [dist-0.5*scaling, dist+0.5*scaling],'k', linewidth = 1)
        #plt.plot([arrival-seconds_prior_arrival+window_length_in_sec, arrival-seconds_prior_arrival+window_length_in_sec], 
        #         [dist-0.5*scaling, dist+0.5*scaling],'k', linewidth = 1)

    '''
    st_select = obspy.Stream()
    offset_select = []
    for tr, dist, time_shift in zip(st, offset, values):
        if np.abs(time_shift)< time_residual_threshold:
            st_select += tr
            offset_select.append(dist)
    
    plot_waveform_section(ax, st_select, offset_select, scale=scaling,colors='r')
    '''

    # synthetics from stf
    stf_folder = comm.waveforms.get_waveform_folder(
        event_name=event_name, data_type="stf",
        tag_or_iteration='ITERATION_%s'%iteration_name)

    st_syn = obspy.Stream()
    for station in stations_info:
        synfile = os.path.join(stf_folder, 
                               stations_info[station]["input_file"].split('/')[-1])
        tr = obspy.read(synfile)
        st_syn += tr
    plot_waveform_section(ax, st_syn, offset, scale=scaling,colors='b')

    annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False) 



    st_not_select = obspy.Stream()
    for mk, tr in zip(mask,st):
        if not mk:
            st_not_select += tr
    res = [i for i, val in enumerate(mask) if not val] 
    plot_waveform_section(ax, st_not_select, offset[res], scale=scaling,colors='r')
    def on_plot_hover(event):
        # Iterating over each data member plotted
        if event.inaxes == ax:
            xpos = ax.get_xlim()[0]
            for curve in ax.get_lines():
                # Searching which data member corresponds to current mouse position
                if curve.contains(event)[0]:
                    annot.xy = (xpos,curve.get_data()[1][0])
                    annot.set_text(curve.get_label())
                    annot.set_visible(True)
                    fig2.canvas.draw_idle()


    fig2.canvas.mpl_connect('button_press_event', on_plot_hover)


    plt.figure()
    vmax = np.max(np.abs(phase_angle))
    map_object = comm.visualizations.plot_stations_for_event(events[0], iteration_name, 0, color = [])
    x, y = map_object(lngs, lats)
    stations_plot = map_object.scatter(x, y, c=phase_angle, s=35, 
                                       marker="v", zorder=5, cmap='RdBu',
                                       vmin = -vmax, vmax=vmax)
    # Setting the picker overwrites the edgecolor attribute on certain
    # matplotlib and basemap versions. Fix it here.    
    stations_plot._edgecolors = np.array([[0.0, 0.0, 0.0, 1.0]])
    stations_plot._edgecolors_original = "black"
    import matplotlib.pyplot as plt
    cbar = plt.colorbar(stations_plot, extend = "both")
    cbar.ax.set_ylabel('Phase angle (rad) at 10 s')

    '''
    plt.figure()
    plt.plot(st_syn[0].times(),st_syn[0].data)
    plt.plot(st[0].times(),st[0].data)
    print(values[0],amp_values[0])
    '''


    # focal mechanism
    fig=plt.figure()
    from obspy.imaging.beachball import beach
    from obspy.geodetics.base import locations2degrees, gps2dist_azimuth
    from obspy.taup import TauPyModel
    earth_model = TauPyModel("ak135")
    event_info = comm.events.get(event_name)
    print(event_info)
    takeoff_angles = []
    azimuths = []
    for stla,stlo in zip(lats,lngs):
        #epicentral_distance = locations2degrees(event_info["latitude"], event_info["latitude"], stla,stlo)
        epicentral_distance, azimuth, baz = gps2dist_azimuth(event_info["latitude"], event_info["longitude"], 
                                                             stla,stlo)
        tt = earth_model.get_travel_times(distance_in_degree=epicentral_distance, 
                                          source_depth_in_km=event_info["depth_in_km"],
                                          phase_list=["P"])
        if tt:
            takeoff_angles.append(tt[0].takeoff_angle)
            azimuths.append(azimuth)

    map_object = comm.visualizations.plot_stations_for_event(events[0], iteration_name, 0, color = [])
    x, y = map_object(lngs[0], lats[0])
    focmec = [event_info["m_rr"], event_info["m_tt"], event_info["m_pp"], event_info["m_rt"],
              event_info["m_rp"], event_info["m_tp"]]
    # Attempt to calculate the best beachball size.
    beachball_size=0.1
    width = max((map_object.xmax - map_object.xmin,
                 map_object.ymax - map_object.ymin)) * beachball_size
    b = beach(focmec, xy=(x, y), width=width, linewidth=1., facecolor="red")

    print(b)
    b.set_zorder(200000000)



    from lasif.tools.pyrocko import beachball
    from lasif.tools.pyrocko import orthodrome
    from lasif.tools.pyrocko import cake
    #axes2=map_object
    axes2 = fig.add_axes(map_object.ax.get_position())#, projection='lambert', label='pol', frameon=False)
    ##axes2.set_aspect('equal')
    ##axes2.set_theta_direction(-1)
    ##axes2.set_theta_zero_location("N")

    axes2.set_axis_off()
    axes2.set_xlim(-1.1, 1.1)
    axes2.set_ylim(-1.1, 1.1)
    axes2.set_aspect('equal')
    map_object.ax.add_collection(b)


    projection='lambert'
    mt = [event_info["m_tt"], event_info["m_pp"], event_info["m_rr"], -event_info["m_tp"],
          event_info["m_rt"], -event_info["m_rp"]]
    beachball.plot_beachball_mpl(
        mt, axes2,
        position=(0, 0),
        size=2.0,
        color_t='red',
        linewidth=1.,
        projection=projection,
        size_units='data')

    for takeoff_angle, azimuth in zip(takeoff_angles,azimuths):
        length_on_equatorial_plane = np.tan(np.deg2rad(takeoff_angle/2));
        x_piecring = length_on_equatorial_plane*np.cos(np.deg2rad(azimuth) - np.pi/2)
        y_piecring = length_on_equatorial_plane*np.sin(np.deg2rad(azimuth) + np.pi/2)
        axes2.plot(x_piecring,y_piecring,'+',color='k',linewidth=1)

        # to spherical coordinates, r, theta, phi in radians
        rtp = np.array([[1., np.deg2rad(takeoff_angle), np.deg2rad(90.-azimuth)]])
        # to 3D coordinates (x, y, z)
        points = beachball.numpy_rtp2xyz(rtp)
        # project to 2D with same projection as used in beachball
        x_arr, y_arr = beachball.project(points, projection='lambert').T

        axes2.plot(x_arr, y_arr,'+', color='r',linewidth=1)


    # earth model and phase for takeoff angle computations
    mod = cake.load_model('ak135-f-average.vf')#'ak135-f-continental.m')
    phases = cake.PhaseDef.classic('P')
    slat = event_info["latitude"]
    slon = event_info["longitude"]
    rdepth = 0.
    sdepth = event_info["depth_in_km"]*1e3
    for rlat, rlon in zip(lats,lngs):
        distance = orthodrome.distance_accurate50m(slat, slon, rlat, rlon)
        rays = mod.arrivals(
            phases=cake.PhaseDef('P'),
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
        x, y = beachball.project(points, projection=projection).T

        axes2.plot(x, y, '+', color='b', ms=10., mew=2.0, mec='black', mfc='none')

    #axes2.set_rmax(np.pi / 2.)
    #axes2.set_yticks([0, np.pi/6., 2.*np.pi/6., np.pi/2.])




    plt.show()
    """


@command_group("Plotting")
def lasif_plot_synthetic_from_stf(parser, args):
    """
    Plot seismic waveform gather for preprocessed data and synthetic waveforms
    """
    parser.add_argument("iteration_name", help="name of the iteration")
    parser.add_argument("event_name", help="name of the event to plot waveform gather")
    parser.add_argument("--components", default="ENZ", help="list of components to plot, examples: ENZ, or Z or RTZ")
    parser.add_argument("--scaling", default=0.5,
                        help="Float for scaling the waveforms, 0.5 by default")
    parser.add_argument("--raw", default="True", choices=["False", "True"],
                        help="For additionally plotting initial synthetics on top, True by default")
    parser.add_argument("--plot_window", default="True", choices=["False", "True"],
                        help="For additionally plotting the phase window defined in the iteration, False by default")
    parser.add_argument("--phase", default="P",
                        help="Name of the seismic phase you want to window for the plot_window option"
                        "default: P")

    args = parser.parse_args(args)
    iteration_name = args.iteration_name
    event_name = args.event_name
    scale = args.scaling
    raw = args.raw
    comps = args.components
    components=[]
    for comp in comps:
        components.append(comp)
    plot_window = args.plot_window
    phase = args.phase

    comm = _find_project_comm(".", args.read_only_caches)
    comm.visualizations.plot_synthetic_from_stf(event_name, iteration_name, 
                                                components, scaling = float(scale), plot_raw=raw, 
                                                plot_window=plot_window, Phase=phase)

    import matplotlib.pyplot as plt
    plt.show()


@command_group("Iteration Management")
def lasif_plot_q_model(parser, args):
    """
    Plots the Q model for a given iteration.
    """
    parser.add_argument("iteration_name", help="name of iteration")
    args = parser.parse_args(args)
    iteration_name = args.iteration_name

    comm = _find_project_comm(".", args.read_only_caches)
    comm.iterations.plot_Q_model(iteration_name)

    import matplotlib.pyplot as plt
    plt.show()


@command_group("Plotting")
def lasif_plot_window_statistics(parser, args):
    """
    Plot the selected windows.
    """
    parser.add_argument("iteration_name", help="name of the iteration")
    args = parser.parse_args(args)

    iteration_name = args.iteration_name

    comm = _find_project_comm(".", args.read_only_caches)

    if args.combine:
        comm.visualizations.plot_window_statistics(
            iteration=iteration_name, ax=None, show=True)


@command_group("Plotting")
def lasif_plot_windows(parser, args):
    """
    Plot the selected windows.
    """
    parser.add_argument("iteration_name", help="name of the iteration")
    parser.add_argument("event_name", help="name of the event")
    parser.add_argument(
        "--combine",
        help="Create a combined plot for all windows of that event.",
        action="store_true")
    parser.add_argument("--distance_bins", type=int,
                        help="The number of bins on the distance axis for "
                             "the combined plot.",
                        default=500)
    args = parser.parse_args(args)

    iteration_name = args.iteration_name
    event_name = args.event_name

    comm = _find_project_comm(".", args.read_only_caches)

    if args.combine:
        comm.visualizations.plot_windows(event=event_name,
                                         iteration=iteration_name, ax=None,
                                         distance_bins=args.distance_bins,
                                         show=True)
    else:
        output_folder = comm.project.get_output_folder(
            type="plotted_windows",
            tag="Iteration_%s__%s" % (event_name, iteration_name))

        window_manager = comm.windows.get(event_name, iteration_name)
        for window_group in window_manager:
            window_group.plot(show=False, filename=os.path.join(output_folder,
                                                                "%s.png" % window_group.channel_id))
            sys.stdout.write(".")
            sys.stdout.flush()
        print("\nDone")

        print(("Done. Written output to folder %s." % output_folder))


@command_group("Project Management")
def lasif_validate_data(parser, args):
    """
    Validate the data currently in the project.

    This commands walks through all available data and checks it for validity.
    It furthermore does some sanity checks to detect common problems. These
    should be fixed.

    By default is only checks some things. A full check is recommended but
    potentially takes a very long time.

    Things the command does:

    Event files:
        * Validate against QuakeML 1.2 scheme.
        * Make sure they contain at least one origin, magnitude and focal
          mechanism object.
        * Check for duplicate ids amongst all QuakeML files.
        * Some simply sanity checks so that the event depth is reasonable and
          the moment tensor values as well. This is rather fragile and mainly
          intended to detect values specified in wrong units.
    """
    parser.add_argument(
        "--station_file_availability",
        help="asserts that all waveform files have corresponding station "
        "files. Very slow.",
        action="store_true")
    parser.add_argument(
        "--raypaths", help="assert that all raypaths are within the "
        "set boundaries. Very slow.", action="store_true")
    parser.add_argument(
        "--waveforms", help="asserts that waveforms for one event have only "
        "a single location and channel type. Fast.", action="store_true")
    parser.add_argument(
        "--station_files", help="check the validity of StationXML files by simulating a correction response"
        "if exceptions araised, will try do download a new file"
        "Slow.", action="store_true")

    parser.add_argument("--full", help="run all validations.",
                        action="store_true")

    args = parser.parse_args(args)
    full_check = args.full
    station_file_availability = args.station_file_availability
    raypaths = args.raypaths
    waveforms = args.waveforms
    stationxml = args.station_files

    # If full check, check everything.
    if full_check:
        station_file_availability = True
        raypaths = True
        waveforms = True
        stationxml = True

    comm = _find_project_comm(".", args.read_only_caches)
    comm.validator.validate_data(
        station_file_availability=station_file_availability,
        raypaths=raypaths, waveforms=waveforms,
        stationxml_file_check=stationxml)
    

@command_group("Project Management")
def lasif_validate_station_files(parser, args):
    """
    Validate the station files currently in the project.
    """
    args = parser.parse_args(args)
    comm = _find_project_comm(".", args.read_only_caches)
    comm.validator.validate_StationXML_files()


@command_group("Iteration Management")
def lasif_iteration_status(parser, args):
    """
    Query the current status of an iteration.
    """
    parser.add_argument("iteration_name", help="name of the iteration")
    args = parser.parse_args(args)
    iteration_name = args.iteration_name

    comm = _find_project_comm(".", args.read_only_caches)
    status = comm.query.get_iteration_status(iteration_name)
    iteration = comm.iterations.get(iteration_name)

    print(("Iteration %s is defined for %i events:" % (iteration_name,
                                                      len(iteration.events))))
    for event in sorted(status.keys()):
        st = status[event]
        print(("\t%s" % event))

        print(("\t\t%.2f %% of the events stations have picked windows" %
              (st["fraction_of_stations_that_have_windows"] * 100)))
        if st["missing_raw"]:
            print(("\t\tLacks raw data for %i stations" %
                  len(st["missing_raw"])))
        if st["missing_processed"]:
            print(("\t\tLacks processed data for %i stations" %
                  len(st["missing_processed"])))
        if st["missing_synthetic"]:
            print(("\t\tLacks synthetic data for %i stations" %
                  len(st["missing_synthetic"])))


def lasif_tutorial(parser, args):
    """
    Open the tutorial in a webbrowser.
    """
    parser.parse_args(args)

    import webbrowser
    webbrowser.open("http://krischer.github.io/LASIF/")


def lasif_calculate_constant_q_model(parser, args):
    """
    Calculate a constant Q model useable by SES3D.
    """
    from lasif.tools import Q_discrete

    parser.add_argument("min_period", type=float,
                        help="minimum period for the constant frequency band")
    parser.add_argument("max_period", type=float,
                        help="maximum period for the constant frequency band")
    args = parser.parse_args(args)

    weights, relaxation_times, = Q_discrete.calculate_Q_model(
        N=3,
        f_min=1.0 / args.max_period,
        f_max=1.0 / args.min_period,
        iterations=10000,
        initial_temperature=0.1,
        cooling_factor=0.9998)


def lasif_debug(parser, args):
    """
    Print information LASIF can gather from a list of files.
    """
    parser.add_argument(
        "files", help="filenames to print debug information about", nargs="+")
    args = parser.parse_args(args)
    comm = _find_project_comm(".", args.read_only_caches)

    for filename in args.files:
        filename = os.path.relpath(filename)
        if not os.path.exists(filename):
            print(("{red}Path '{f}' does not exist.{reset}\n".format(
                f=filename, red=colorama.Fore.RED,
                reset=colorama.Style.RESET_ALL)))
            continue
        print(("{green}Path '{f}':{reset}".format(
            f=filename, green=colorama.Fore.GREEN,
            reset=colorama.Style.RESET_ALL)))

        try:
            info = comm.query.what_is(filename)
        except LASIFError as e:
            info = "Error: %s" % e.message

        print(("\t" + info))
        print("")


@command_group("Misc")
def lasif_serve(parser, args):
    """
    Launches the LASIF webinterface.
    """
    parser.add_argument("--port", default=8008, type=int,
                        help="Port of the webserver.")

    parser.add_argument("--nobrowser", help="Do not open a webbrowser.",
                        action="store_true")
    parser.add_argument("--debug", help="Turn on debugging. Implies "
                                        "'--nobrowser'.",
                        action="store_true")
    parser.add_argument(
        "--open_to_outside",
        help="By default the website can only be opened from the current "
             "computer. Use this argument to access it from any other "
             "computer on the network.",
        action="store_true")
    args = parser.parse_args(args)
    port = args.port
    nobrowser = args.nobrowser
    debug = args.debug
    open_to_outside = args.open_to_outside

    if debug:
        nobrowser = True

    comm = _find_project_comm(".", args.read_only_caches)

    if nobrowser is False:
        import webbrowser
        import threading

        threading.Timer(
            1.0, lambda: webbrowser.open("http://localhost:%i" % port)).start()

    from lasif.webinterface.server import serve
    serve(comm, port=port, debug=debug, open_to_outside=open_to_outside)


def _get_cmd_description(fct):
    """
    Convenience function extracting the first line of a docstring.
    """
    try:
        return fct.__doc__.strip().split("\n")[0].strip()
    except:
        return ""


def _print_generic_help(fcts):
    """
    Small helper function printing a generic help message.
    """
    print((100 * "#"))
    header = ("{default_style}LASIF - Large Scale Seismic "
              "{inverted_style}Inversion"
              "{default_style} Framework{reset_style}  [Version {version}]"
              .format(
                  default_style=colorama.Style.BRIGHT + colorama.Fore.WHITE +
                  colorama.Back.BLACK,
                  inverted_style=colorama.Style.BRIGHT + colorama.Fore.BLACK +
                  colorama.Back.WHITE,
                  reset_style=colorama.Style.RESET_ALL,
                  version=lasif.__version__))
    print(("    " + header))
    print("    http://krischer.github.io/LASIF")
    print((100 * "#"))
    print(("\n{cmd}usage: lasif [--help] COMMAND [ARGS]{reset}\n".format(
        cmd=colorama.Style.BRIGHT + colorama.Fore.RED,
        reset=colorama.Style.RESET_ALL)))

    # Group the functions. Functions with no group will be placed in the group
    # "Misc".
    fct_groups = {}
    for fct_name, fct in fcts.items():
        group_name = fct.group_name if hasattr(fct, "group_name") else "Misc"
        fct_groups.setdefault(group_name, {})
        fct_groups[group_name][fct_name] = fct

    # Print in a grouped manner.
    for group_name in sorted(fct_groups.keys()):
        print(("{0:=>25s} Functions".format(" " + group_name)))
        current_fcts = fct_groups[group_name]
        for name in sorted(current_fcts.keys()):
            print(("%s  %32s: %s%s%s" % (colorama.Fore.YELLOW, name,
                                        colorama.Fore.BLUE,
                                        _get_cmd_description(fcts[name]),
                                        colorama.Style.RESET_ALL)))
    print("\nTo get help for a specific function type")
    print("\tlasif help FUNCTION  or\n\tlasif FUNCTION --help")


def _get_argument_parser(fct):
    """
    Helper function to create a proper argument parser.
    """
    parser = argparse.ArgumentParser(
        prog="lasif %s" % fct.__name__.replace("lasif_", ""),
        description=_get_cmd_description(fct))

    parser.add_argument(
        "--ipdb",
        help="If true, a debugger will be launched upon encountering an "
             "exception. Requires ipdb.",
        action="store_true")

    # Exceptions. If any are missed, its not mission critical but just
    # less nice.
    exceptions = ["lasif_tutorial", "lasif_init_project",
                  "lasif_build_all_caches"]

    if fct.__name__ in exceptions:
        return parser

    # Otherwise add the option to add caches in read-only mode.
    parser.add_argument("--read_only_caches",
                        help="sets all caches to read-only",
                        action="store_true")
    return parser


def _get_functions():
    """
    Get a list of all CLI functions defined in this file.
    """
    # Get all functions in this script starting with "lasif_".
    fcts = {fct_name[len(FCT_PREFIX):]: fct for (fct_name, fct) in
            globals().items()
            if fct_name.startswith(FCT_PREFIX) and hasattr(fct, "__call__")}
    return fcts


def main():
    """
    Main entry point for the LASIF command line interface.

    Essentially just dispatches the different commands to the corresponding
    functions. Also provides some convenience functionality like error catching
    and printing the help.
    """
    fcts = _get_functions()
    # Parse args.
    args = sys.argv[1:]

    if len(args) == 1 and args[0] == "--version":
        print(("LASIF version %s" % lasif.__version__))
        sys.exit(0)

    # Print the generic help/introduction.
    if not args or args == ["help"] or args == ["--help"]:
        _print_generic_help(fcts)
        sys.exit(0)

    # Use lowercase to increase tolerance.
    fct_name = args[0].lower()

    further_args = args[1:]
    # Map "lasif help CMD" to "lasif CMD --help"
    if fct_name == "help":
        if further_args and further_args[0] in fcts:
            fct_name = further_args[0]
            further_args = ["--help"]
        else:
            sys.stderr.write("lasif: Invalid command. See 'lasif --help'.\n")
            sys.exit(1)

    # Unknown function.
    if fct_name not in fcts:
        sys.stderr.write("lasif: '{fct_name}' is not a LASIF command. See "
                         "'lasif --help'.\n".format(fct_name=fct_name))
        # Attempt to fuzzy match commands.
        close_matches = sorted(difflib.get_close_matches(fct_name, list(fcts.keys()),
                                                         n=4))
        if len(close_matches) == 1:
            sys.stderr.write("\nDid you mean this?\n\t{match}\n".format(
                match=close_matches[0]))
        elif close_matches:
            sys.stderr.write(
                "\nDid you mean one of these?\n    {matches}\n".format(
                    matches="\n    ".join(close_matches)))

        sys.exit(1)

    func = fcts[fct_name]

    # Make sure that only MPI enabled functions are called with MPI.
    if MPI.COMM_WORLD.size > 1:
        if not hasattr(func, "_is_mpi_enabled") or \
                func._is_mpi_enabled is not True:
            if MPI.COMM_WORLD.rank != 0:
                return
            sys.stderr.write("'lasif %s' must not be called with MPI.\n" %
                             fct_name)
            return

    # Create a parser and pass it to the single function.
    parser = _get_argument_parser(func)

    # Now actually call the function.
    try:
        func(parser, further_args)
    except LASIFCommandLineException as e:
        print((colorama.Fore.YELLOW + ("Error: %s\n" % str(e)) +
              colorama.Style.RESET_ALL))
        sys.exit(1)
    except Exception as e:
        args = parser.parse_args(further_args)
        # Launch ipdb debugger right at the exception point if desired.
        # Greatly eases debugging things. Requires ipdb to be installed.
        if args.ipdb:
            import ipdb  # NOQA
            _, _, tb = sys.exc_info()
            traceback.print_exc()
            ipdb.post_mortem(tb)
        else:
            print((colorama.Fore.RED))
            traceback.print_exc()
            print((colorama.Style.RESET_ALL))

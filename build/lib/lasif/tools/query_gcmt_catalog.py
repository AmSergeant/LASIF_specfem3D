#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:copyright:
    Lion Krischer (krischer@geophysik.uni-muenchen.de), 2014
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""
import glob
import inspect
import numpy as np
import obspy
from obspy.core.event import Catalog
import os
import random
from scipy.spatial import cKDTree

EARTH_RADIUS = 6371.00

from lasif.utils import get_event_filename


class SphericalNearestNeighbour(object):
    """
    Spherical nearest neighbour queries using scipy's fast
    kd-tree implementation.
    """
    def __init__(self, data):
        cart_data = self.spherical2cartesian(data)
        self.data = data
        self.kd_tree = cKDTree(data=cart_data, leafsize=10)

    def query(self, points, k=10):
        points = self.spherical2cartesian(points)
        d, i = self.kd_tree.query(points, k=k)
        return d, i

    @staticmethod
    def spherical2cartesian(data):
        """
        Converts an array of shape (x, 2) containing latitude/longitude
        pairs into an array of shape (x, 3) containing x/y/z assuming a
        radius of one for points on the surface of a sphere.
        """
        lat = data[:, 0]
        lng = data[:, 1]
        # Convert data from lat/lng to x/y/z, assume radius of 1
        colat = 90 - lat
        cart_data = np.empty((lat.shape[0], 3))

        cart_data[:, 0] = np.sin(np.deg2rad(colat)) * \
            np.cos(np.deg2rad(lng))
        cart_data[:, 1] = np.sin(np.deg2rad(colat)) * \
            np.sin(np.deg2rad(lng))
        cart_data[:, 2] = np.cos(np.deg2rad(colat))

        return cart_data

def download_GCMT_catalog_from_url(year):
    """
    Helper function downloading the GCMT catalog as a ndk file for one year 
    from https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_MONTHLY/

    save monthly ndk files within Lasif to data/GCMT_Catalog

    :param year: The year to download.
    :type year: str
    """
    import requests
    parent_url = 'https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_MONTHLY/'

    month_list = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']
    data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(
        inspect.getfile(inspect.currentframe())))), "data", "GCMT_Catalog")
    for month in month_list:

        ndk_file = '%s%s.ndk'%(month,year[-2:])
        child_url = os.path.join(parent_url, year, ndk_file)

        child_folder = os.path.join(data_dir, year)

        try:
            r = requests.get(child_url)

            if r.status_code ==200:
                print(("Requesting CMT catalog at %s"%child_url))
                if not os.path.exists(child_folder):
                    os.makedirs(child_folder)   
                with open(os.path.join(child_folder, ndk_file), 'w') as f:
                    f.write(r.content)
        except:
            print(("Could not recognize url %s"%child_url))
            continue


def _get_years_from_GCMT_catalog(data_dir = None):
    """
    Helper function reading the GCMT data shipping with LASIF.

    :param data_dir: path to the GCMT catalogs in Lasif
    :type data_dir: str, optional
    """

    if data_dir is None:
        data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(
            inspect.getfile(inspect.currentframe())))), "data", "GCMT_Catalog")
    available_years = [_i for _i in os.listdir(data_dir) if _i.isdigit()]
    available_years.sort()

    return available_years




def _read_GCMT_catalog(min_year=None, max_year=None):
    """
    Helper function reading the GCMT data shipping with LASIF.

    :param min_year: The minimum year to read.
    :type min_year: int, optional
    :param max_year: The maximum year to read.
    :type max_year: int, optional
    """
    # easier tests
    if min_year is None:
        min_year = 0
    else:
        min_year = int(min_year)
    if max_year is None:
        max_year = 3000
    else:
        max_year = int(max_year)

    data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(
        inspect.getfile(inspect.currentframe())))), "data", "GCMT_Catalog")
    available_years = _get_years_from_GCMT_catalog(data_dir)

    years = np.arange(min_year, max_year+1, 1)
    years_to_download = [str(year )for year in years if str(year) not in available_years]
    if years_to_download:
        for year in years_to_download:
            try:
                download_GCMT_catalog_from_url(year)
            except:
                print(("Could not download GMT Catalog for year %d"%year))
                continue
        available_years = _get_years_from_GCMT_catalog(data_dir)
    print(("LASIF currently contains GCMT data from %s to %s/%i." % (
        available_years[0], available_years[-1],
        len(glob.glob(os.path.join(data_dir, available_years[-1], "*.ndk*"))))))

    goal_available_years = \
        [_i for _i in available_years if (min_year <= int(_i) <= max_year)]
    goal_available_years.sort()

    print("Parsing the GCMT catalog. This might take a while...")
    cat = Catalog()
    for year in goal_available_years:
        print(("\tReading year %s ..." % year))
        for filename in glob.glob(os.path.join(data_dir, str(year),
                                               "*.ndk*")):
            cat += obspy.read_events(filename, format="ndk")

    return cat


def add_new_events(comm, count, min_year=None, max_year=None, 
                   statistical_selection=True, min_magnitude=None, max_magnitude=None):

    import warnings
    warnings.filterwarnings("ignore")
    
    # get the configuration in config file 
    # to gather information on requested events
    proj = comm.project
    ds = proj.config["download_settings"]   
    config = ds["configuration"]
    min_dist = ds["minimum_epicentral_distance_in_degree"]
    max_dist = ds["maximum_epicentral_distance_in_degree"]
    min_depth = ds["minimum_depth_in_km"]
    max_depth = ds["maximum_depth_in_km"]
    if min_magnitude is None:
        min_magnitude = ds["minimum_magnitude"]
    if max_magnitude is None:
        max_magnitude = ds["maximum_magnitude"]
    if config == "teleseismic":
        query_inside_domain = 0
    else:
        query_inside_domain = 1
    # The minimum acceptable distance to the next closest event in km 
    # for statistical selection on hypocenters
    threshold_distance_in_km = ds["minimum_adjacent_distance_in_km"]
    if min_year is None:
        min_year = 0
    if max_year is None:
        max_year = 3000

    ####### print to the terminal useful infos ######
    print("--------------------------------------------------------------")
    print(("Will look for %d events in %s configuration:"%(count,config)))
    print(("\t happening bewteen %s and %s"%(min_year,max_year)))
    if config == "teleseismic":
        print("\t outside of the chosen domain")
        print(("\t %0.2f <= epicentral_distance <= %0.2f"%(min_dist, max_dist)))
    else:
        print("\t inside the chosen domain")
    print(("\t %0.2f <= magnitude <= %0.2f"%(min_magnitude, max_magnitude)))
    print(("\t %0.2f <= depth_in_km <= %0.2f"%(min_depth, max_depth)))
    if statistical_selection:
        print(("\t will apply a statistical selection to keep %d events"%count))
        print(("\t\t with adjacent hypocentral distances less than %0.2f km"%threshold_distance_in_km))
        print("\t\t and will make sure the events are more than one day apart from each other")
    print("--------------------------------------------------------------")
    #################################################

    # Get the catalog.
    download_from_iris = False
    Flag_filter = 1
    cat = []
    GCMT_years = _get_years_from_GCMT_catalog()
    if int(min_year) >= int(GCMT_years[0]) and int(max_year) <= int(GCMT_years[-1]):
        # look at events in the GCMT catalog stored in Lasif or download it from Columbia
        cat = _read_GCMT_catalog(min_year=min_year, max_year=max_year)
        print("")
        print(("--> Have found %d events in catalogs, need filter"%len(cat)))
    if not cat:
        download_from_iris = True
        # will look for event catalog at IRIS
        print("")
        print("No event found in stored GCMT catalog")
        print("\t\tWill try at IRIS")
        print("")
        from obspy.clients.fdsn import Client
        from obspy.core import UTCDateTime
        client = Client("IRIS")
        t1 = UTCDateTime(int(min_year),1,1,0,0,0)
        t2 = UTCDateTime(int(max_year),12,31,23,59,59)
        if config == "teleseismic":
            Flag_filter = 0
            central_point = comm.query.center()
            cat = client.get_events(starttime=t1, endtime=t2,
                                    minmagnitude=min_magnitude, maxmagnitude=max_magnitude,
                                    minradius=min_dist,maxradius=max_dist,
                                    latitude=central_point.latitude, longitude=central_point.longitude,
                                    mindepth=min_depth, maxdepth=max_depth, 
                                    orderby="time-asc")
            print(("--> Have found %i valid filtered events with requested event criteria." % len(cat)))
        else:
            cat = client.get_events(starttime=t1, endtime=t2,
                                    minmagnitude=min_magnitude, maxmagnitude=max_magnitude,
                                    mindepth=min_depth, maxdepth=max_depth, 
                                    orderby="time-asc")
        print(("--> Have found in total %d events in catalogs, need filter"%len(cat)))
        
        
            

    # Filter with the magnitudes
    if Flag_filter:
        print("\tFiltering the events based on magnitude")
        cat = cat.filter("magnitude >= %.2f" % min_magnitude,
                         "magnitude <= %.2f" % max_magnitude)

    if not cat:
        print("--> No event left to query, Stopping")
    else:
        # Filtering catalog to only contain events in or outside the domain.
        if query_inside_domain==1:
            print("\tFiltering to only include events inside domain...")
            # Coordinates and the Catalog will have the same order!
            temp_cat = Catalog()
            coordinates = []
            azimuths = []
            from obspy.geodetics import gps2dist_azimuth, locations2degrees
            central_point = comm.query.center()
            for event in cat:
                org = event.preferred_origin() or event.origins[0]
                if not comm.query.point_in_domain(org.latitude, org.longitude):
                    continue
                # filter with hypocentral depth
                event_depth_in_km = event.origins[0]["depth"]*1e-3
                if event_depth_in_km > max_depth or event_depth_in_km < min_depth:
                    continue
                dist_in_deg, azimuth, baz = gps2dist_azimuth(org["latitude"], 
                                                             org["longitude"],
                                                             central_point.latitude,
                                                             central_point.longitude)
                azimuths.append(azimuth)
                temp_cat.events.append(event)
                coordinates.append((org.latitude, org.longitude))      
        else:
            from obspy.geodetics import gps2dist_azimuth, locations2degrees
            central_point = comm.query.center()
            coordinates = []
            azimuths = []
            if Flag_filter:
                print("\tFiltering to only events outside domain within epicentral distance boundaries...")
                temp_cat = Catalog()
                for event in cat:
                    org = event.preferred_origin() or event.origins[0]
                    dist_in_deg = locations2degrees(org["latitude"], 
                                                    org["longitude"],
                                                    central_point.latitude,
                                                    central_point.longitude)
                    _dist, azimuth, _baz = gps2dist_azimuth(org["latitude"], 
                                                            org["longitude"],
                                                            central_point.latitude,
                                                            central_point.longitude)
                    if dist_in_deg<min_dist:
                        continue
                    if dist_in_deg>max_dist:
                        continue
                    # filter with hypocentral depth
                    event_depth_in_km = event.origins[0]["depth"]*1e-3
                    if event_depth_in_km > max_depth or event_depth_in_km < min_depth:
                        continue
                    temp_cat.events.append(event)
                    coordinates.append((org.latitude, org.longitude))
                    azimuths.append(azimuth)
            else: #this selection was already made when calling the catalog from IRIS
                temp_cat = cat
                for event in cat:
                    org = event.preferred_origin() or event.origins[0]
                    dist_in_deg, azimuth, baz = gps2dist_azimuth(org["latitude"], 
                                                                 org["longitude"],
                                                                 central_point.latitude,
                                                                 central_point.longitude)
                    coordinates.append((org.latitude, org.longitude))
                    azimuths.append(azimuth)

        # remove any duplicates from cat
        event_names = []
        for event in temp_cat:
            event_names.append(get_event_filename(event, "GCMT").split(".xml")[0])
        _u, index_unique = np.unique(np.array(event_names), return_index=True)
        cat = Catalog()
        temp_azimuth = []
        for idx in index_unique:
            cat.events.append(temp_cat[idx])
            temp_azimuth.append(azimuths[idx])
        azimuths = temp_azimuth
        
        
        if not cat:
            print("--> No event left to query, Stopping")
        else:

            if Flag_filter:
                print(("--> %i valid filtered events remain with requested event criteria." % len(cat)))

            existing_events = list(comm.events.get_all_events().values())
            # Get the coordinates of all existing events.
            existing_coordinates = [
                (_i["latitude"], _i["longitude"]) for _i in existing_events]
            existing_origin_times = [_i["origin_time"] for _i in existing_events]

            chosen_events = []
            if statistical_selection:
                print(("\tRandom selection of %i events based on statistical hypocentral location distributions"%count))
                # Special case handling in case there are no preexisting events.
                if not existing_coordinates:
                    idx = random.randint(0, len(cat) - 1)

                    chosen_events.append(cat[idx])
                    del cat.events[idx]
                    existing_coordinates.append(coordinates[idx])
                    del coordinates[idx]

                    _t = cat[idx].preferred_origin() or cat[idx].origins[0]
                    existing_origin_times.append(_t.time)

                    count -= 1


                while count:
                    print("\t\tStarting selection process ...")

                    if not coordinates:
                        print("\tNo events left to select from. Stoping here.")
                        break
                    # Build kdtree and query for the point furthest away from any other
                    # point.
                    kdtree = SphericalNearestNeighbour(np.array(existing_coordinates))
                    distances = kdtree.query(np.array(coordinates), k=1)[0]
                    idx = np.argmax(distances)

                    event = cat[idx]
                    coods = coordinates[idx]
                    del cat.events[idx]
                    del coordinates[idx]

                    # Actual distance.
                    distance = EARTH_RADIUS * distances[idx]

                    if distance < threshold_distance_in_km:
                        print(("\tNo events left with distance to the next closest event "
                              "of more then %.1f km. Stoping here." %
                              threshold_distance_in_km))
                        break

                    # Make sure it did not happen within one day of an existing event.
                    # This should also filter out duplicates.
                    _t = event.preferred_origin() or event.origins[0]
                    origin_time = _t.time

                    if min([abs(origin_time - _i) for _i in existing_origin_times]) < \
                            86400:
                        print("\tSelected event temporally to close to existing event. "
                              "Will not be chosen. Skipping to next event.")
                        continue

                    print(("\tSelected event with the next closest event being %.1f km "
                          "away." % distance))

                    chosen_events.append(event)
                    existing_coordinates.append(coods)
                    count -= 1
            else:
                if len(cat)>count:
                    chosen_events = []
                    # get une uniform random selection of count events, which are sorted by origin time
                    from scipy.stats import randint as sp_randint

                    # first check if events from catalog are already in the project database, if yes, delete them
                    if existing_events:
                        temp_cat = Catalog()
                        temp_azimuth = []
                        for event, azimuth in zip(cat, azimuths):
                            gcmt_event_name = get_event_filename(event, "GCMT").split(".xml")[0]
                            exist_event = [True for evt in existing_events if gcmt_event_name in evt["event_name"]]
                            if not exist_event:
                                temp_cat.events.append(event)
                                temp_azimuth.append(azimuth)
                        
                        if temp_cat:
                            print(("--> %d events already stored in the database, disregard them"%(np.abs(len(cat)-len(temp_cat)))))
                            cat = temp_cat
                            azimuths = temp_azimuth
                    
                    
                    # sort by azimuth
                    azimuths = np.array(azimuths)
                    Isort = np.argsort(azimuths)
                    sorted_cat = Catalog()
                    for i in Isort:
                        sorted_cat.events.append(cat[i])
                    cat = sorted_cat

                    # random selection
                    available_count = np.min([len(cat), count])
                    print(("\tRandom selection of %i events from a uniform distribution on event azimuth"%available_count))
                    idx = sp_randint.rvs(0, len(cat), size=available_count, random_state=0)
                    
                    # check repeatitive idx
                    idx_unique,  unique_index, unique_counts = np.unique(np.array(idx), return_index=True, return_counts=True)
                    repeated_idx = [index for index, index_count in zip(idx_unique, unique_counts) if index_count>1]
                    idx_unique = list(idx_unique)
                    if repeated_idx:
                        for item in repeated_idx:
                            new_idx = item + 1
                            while new_idx in idx_unique:
                                new_idx += 1
                                if new_idx >= len(cat):
                                    new_idx = np.random.randint(0,len(cat)-1,1)[0]
                            idx_unique.append(new_idx)
                        idx = idx_unique
                        
                    for index in idx:
                        event = cat[index]
                        chosen_events.append(event)
                else:
                    chosen_events = cat

            print(("--> Selected %i events out of %i initially filtered events." %(len(chosen_events),len(cat))))
            
            # read quakeml cat if downloaded from IRIS
            if download_from_iris is True:
                print("Will download quakeMl from IRIS, this might take a while ...")
                from lasif.scripts.iris2quakeml import sparse_quakeml_from_iris
                import urllib.request, urllib.error, urllib.parse
                temp_cat = Catalog()
                for event in chosen_events:
                    event_id = str(event.resource_id).split('=')[-1]
                    url = "http://ds.iris.edu/spudservice/item?eventid="+event_id
                    r = urllib.request.urlopen(url)
                    if r.code != 200:
                        r.close()
                        continue
                    else:
                        url = [str(line) for line in r.readlines() \
                                       if "http://ds.iris.edu/spudservice/momenttensor/" in str(line)][0].split("\"")[1]
                        url += "/quakeml"
                        try:
                            qml = sparse_quakeml_from_iris(url)
                        except:
                            continue
                    temp_cat.append(qml[0])
                chosen_events = temp_cat
            
            
            # get source mechanism infos for events in catalog
            events = []
            for event in chosen_events:
                origin = event.preferred_origin() or event.origins[0]
                magnitude = event.preferred_magnitude() or event.magnitudes[0]
                focmech = event.focal_mechanisms[0]["moment_tensor"]["tensor"]
                event_info = {"event_name": get_event_filename(event, "GCMT").split('.xml')[0],
                              "latitude": origin["latitude"],
                              "longitude": origin["longitude"],
                              "depth_in_km": origin["depth"]*1e-3,
                              "origin_time": origin["time"],
                              "magnitude": magnitude["mag"],
                              "region": '-',
                              "m_rr": focmech["m_rr"],
                              "m_tt": focmech["m_tt"],
                              "m_pp": focmech["m_pp"],
                              "m_rt": focmech["m_rt"],
                              "m_rp": focmech["m_rp"],
                              "m_tp": focmech["m_tp"]}
                event["event_info"] = event_info
            # if already events in the database, get thes infos too for the events map
            if existing_coordinates:
                for event in existing_events: 
                    events.append(event)
            
            '''
            # query inventory for tations in chosen domain
            t1 = UTCDateTime(int(min_year),1,1,0,0,0)
            t2 = UTCDateTime(int(max_year),12,31,23,59,59)
            inventory = comm.downloads.query_station_inventory(t1, t2)
            '''
            
            # gui for event selection
            from lasif import gui_event_query
            gui_event_query.launch_event_gui(chosen_events, events)
            
                
            print("Updating event cache ...")
            comm.events.update_cache()

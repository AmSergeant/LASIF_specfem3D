#!/usr/bin/env python
# -*- coding: utf-8 -*-


import collections
import colorama
import os
import sys

from .component import Component
from ..domain import GlobalDomain
from .. import LASIFNotFoundError

import warnings


class ValidatorComponent(Component):
    """
    Component responsible for validating data inside a project. Needs access
    to a lot of functionality and should therefore be initialized fairly late.

    :param communicator: The communicator instance.
    :param component_name: The name of this component for the communicator.
    """
    def __init__(self, *args, **kwargs):
        super(ValidatorComponent, self).__init__(*args, **kwargs)
        self._reports = []
        self._total_error_count = 0
        self._inspected_station_files = []

    def _print_ok_message(self):
        """
        Prints a colored OK message when a certain test has been passed.
        """
        ok_string = " %s[%sOK%s]%s" % (
            colorama.Style.BRIGHT, colorama.Style.NORMAL + colorama.Fore.GREEN,
            colorama.Fore.RESET + colorama.Style.BRIGHT,
            colorama.Style.RESET_ALL)
        print(ok_string)

    def _print_fail_message(self):
        """
        Prints a colored fail message when a certain test has been passed.
        """
        fail_string = " %s[%sFAIL%s]%s" % (
            colorama.Style.BRIGHT, colorama.Style.NORMAL + colorama.Fore.RED,
            colorama.Fore.RESET + colorama.Style.BRIGHT,
            colorama.Style.RESET_ALL)
        print(fail_string)

    def _flush_point(self):
        """
        Helper function just flushing a point to stdout to indicate progress.
        """
        sys.stdout.write(".")
        sys.stdout.flush()

    def _add_report(self, message, error_count=1):
        """
        Helper method adding a new error message.
        """
        self._reports.append(message)
        self._total_error_count += error_count
        
    def validate_StationXML_files(self):
        """
        will simulate an instrumental correction, 
        if exception raised, will try to download a new stationwml file 
        from iris and orfeus
        """
        
        from obspy.clients.fdsn import Client
        import obspy
        import numpy as np
        iris_client = Client("IRIS")
        orfeus_client = Client("ORFEUS")
        
        def get_raw_data_filename_for_channel(events, channel):
            """
            return the first raw data file which meets the channel
            """
            for event in events:
                try:
                    wav_cache = self.comm.waveforms.get_waveform_cache(event["event_name"],"raw").get_values()
                except Exception:
                    continue
                wav_for_channelid = [wav for wav in wav_cache if wav["channel_id"] == channel["channel_id"]]
                for value in wav_for_channelid:
                    value["starttime"] = \
                            obspy.UTCDateTime(value["starttime_timestamp"])
                    value["endtime"] = \
                            obspy.UTCDateTime(value["endtime_timestamp"])
                    del value["starttime_timestamp"]
                    del value["endtime_timestamp"]
                wav_filename = [wav["filename"] for wav in wav_for_channelid 
                                if wav["starttime"]>=channel["start_date"] and wav["endtime"]<=channel["end_date"]]
                if wav_filename:
                    return wav_filename[0]
                
        print("Checking all station files ", end=' ')

        all_good = True
        
        channels = self.comm.stations.get_all_channels()
        events = list(self.comm.events.get_all_events().values())
        if not self._inspected_station_files:
            self._inspected_station_files = []
        
        output_units = "VEL"
        pre_filt = (0.004, 0.01, 0.1, 0.5)
        fail_item = 0
        warnings.filterwarnings("ignore")
        for channel in channels:
            station_file = channel["filename"]
            if station_file in self._inspected_station_files:
                continue
            
            if channel["end_date"] is None:
                channel["end_date"] = obspy.UTCDateTime(2030,12,31,0,0,0)
                
            example_data_to_read = get_raw_data_filename_for_channel(events, channel)
            if example_data_to_read:
                st = obspy.read(get_raw_data_filename_for_channel(events, channel))
                tr = st[0].copy()
                        
                if '/StationXML/' in station_file:
                    self._flush_point()
                    self._inspected_station_files.append(station_file)
                    inv = obspy.read_inventory(station_file, format="stationxml")
                    tr.attach_response(inv)
                    try:
                        tr.remove_response(output=output_units, pre_filt=pre_filt,
                                                   zero_mean=False, taper=False)
                        
                        if np.isnan(np.sum(tr.data)):
                            raise Exception
                        print("%s: %sok%s"%(station_file, colorama.Fore.GREEN, colorama.Fore.RESET))
                    except Exception:
                        print("%sProblem with station file %s%s"%(colorama.Fore.YELLOW, station_file, colorama.Fore.RESET))
                        print("\t--> will try to download from IRIS")
                        network, station, location, chan = channel["channel_id"].split('.')
                        try:
                            iris_client.get_stations(network=network, station=station, level="response",
                                                            filename = station_file, format="xml")
                            print("\t%sSuccessfull download%s"%(colorama.Fore.GREEN,colorama.Forc.RESET))
                        except Exception:
                            print("\tNo available station response at IRIS for %s.%s"%(network,station))
                            print("\t--> will try to download from ORFEUS")
                            try:
                                orfeus_client.get_stations(network=network, station=station, level="response",
                                                                filename = station_file, format="xml")
                                print("\t%sSuccessfull download%s"%(colorama.Fore.GREEN, colorama.Fore.RESET))
                            except Exception:
                                print("\tNo available station response at ORFEUS for %s.%s -- %sskipping%s"\
                                      %(network,station,colorama.Fore.RED,colorama.Fore.RESET))
                                fail_item += 1
                                all_good = False
                                self._add_report("Defective station file %s"%station_file)
            
        if all_good:
            self._print_ok_message()
        else:
            print("%s%d/%d station files could not be fixed%s"
                  %(colorama.Fore.RED, fail_item, len(self._inspected_station_files), colorama.Fore.RESET))
            self._print_fail_message()

        
        

    def validate_data(self, station_file_availability=False, raypaths=False,
                      waveforms=False, stationxml_file_check = False):
        """
        Validates all data of the current project.

        This commands walks through all available data and checks it for
        validity.  It furthermore does some sanity checks to detect common
        problems. These should be fixed.

        Event files:
            * Validate against QuakeML 1.2 scheme.
            * Make sure they contain at least one origin, magnitude and focal
              mechanism object.
            * Check for duplicate ids amongst all QuakeML files.
            * Some simply sanity checks so that the event depth is reasonable
              and the moment tensor values as well. This is rather fragile and
              mainly intended to detect values specified in wrong units.
        """
        # get configuration as defined in conf domain 
        # in case of teleseismic configuration, will skip some tests
        proj = self.comm.project
        ds = proj.config["download_settings"]
        config = ds["configuration"]

        # Reset error and report counts.
        self._reports = []
        self._total_error_count = 0

        self._validate_event_files(config)

        # Assert that all waveform files have a corresponding station file.
        if station_file_availability:
            self._validate_station_files_availability()
        else:
            print(("%sSkipping station files availability check.%s" % (
                colorama.Fore.YELLOW, colorama.Fore.RESET)))
        
            

        # Assert that all waveform files have a corresponding station file.
        if waveforms:
            self._validate_waveform_files()
        else:
            print(("%sSkipping waveform file validation.%s" % (
                colorama.Fore.YELLOW, colorama.Fore.RESET)))
            
        # checking station files
        if stationxml_file_check:
            self.validate_StationXML_files()
        else:
            print(("%sSkipping station file validation.%s" % (
                colorama.Fore.YELLOW, colorama.Fore.RESET)))

        # self._validate_coordinate_deduction(ok_string, fail_string,
        # flush_point, add_report)

        files_failing_raypath_test = []
        if raypaths:
            if isinstance(self.comm.project.domain, GlobalDomain):
                print(("%sSkipping raypath checks for global domain...%s" % (
                    colorama.Fore.YELLOW, colorama.Fore.RESET)))
            elif config == "teleseismic":
                print(("%sSkipping raypath checks for teleseismic configuration...%s" % (
                    colorama.Fore.YELLOW, colorama.Fore.RESET)))
            else:
                files_failing_raypath_test = \
                    self.validate_raypaths_in_domain()
        else:
            print(("%sSkipping raypath checks.%s" % (
                colorama.Fore.YELLOW, colorama.Fore.RESET)))
            
        

        # Depending on whether or not the tests passed, report it accordingly.
        if not self._reports:
            print(("\n%sALL CHECKS PASSED%s\n"
                  "The data seems to be valid. If we missed something please "
                  "contact the developers." % (colorama.Fore.GREEN,
                                               colorama.Fore.RESET)))
        else:
            folder = \
                self.comm.project.get_output_folder(
                    type="validation",
                    tag="data_integrity_report")
            filename = os.path.join(folder, "report.txt")
            seperator_string = "\n" + 80 * "=" + "\n" + 80 * "=" + "\n"
            with open(filename, "wt") as fh:
                for report in self._reports:
                    fh.write(report.strip())
                    fh.write(seperator_string)
            print(("\n%sFAILED%s\nEncountered %i errors!\n"
                  "A report has been created at '%s'.\n" %
                  (colorama.Fore.RED, colorama.Fore.RESET,
                   self._total_error_count, os.path.relpath(filename))))
            if files_failing_raypath_test:
                # Put quotes around the filenames
                files_failing_raypath_test = ['"%s"' % _i for _i in
                                              files_failing_raypath_test]
                filename = os.path.join(folder,
                                        "delete_raypath_violating_files.sh")
                with open(filename, "wt") as fh:
                    fh.write("# CHECK THIS FILE BEFORE EXECUTING!!!\n")
                    fh.write("rm ")
                    fh.write("\nrm ".join(files_failing_raypath_test))
                print(("\nSome files failed the raypath in domain checks. A "
                      "script which deletes the violating files has been "
                      "created. Please check and execute it if necessary:\n"
                      "'%s'" % filename))

    def validate_raypaths_in_domain(self):
        """
        Checks that all raypaths are within the specified domain boundaries.

        Returns a list of waveform files violating that assumtion.
        """
        print("Making sure raypaths are within boundaries ", end=' ')

        all_good = True

        # Collect list of files to be deleted.
        files_to_be_deleted = []

        for event_name, event in self.comm.events.get_all_events().items():
            try:
                waveforms = self.comm.waveforms.get_metadata_raw(event_name)
            except LASIFNotFoundError:
                continue
            self._flush_point()
            for station_id, value in \
                    self.comm.query.get_all_stations_for_event(
                        event_name).items():
                network_code, station_code = station_id.split(".")
                # Check if the whole path of the event-station pair is within
                # the domain boundaries.
                if self.is_event_station_raypath_within_boundaries(
                        event_name, value["latitude"], value["longitude"],
                        raypath_steps=12):
                    continue
                # Otherwise get all waveform files for that station.
                waveform_files = [_i["filename"]
                                  for _i in waveforms
                                  if (_i["network"] == network_code) and
                                  (_i["station"] == station_code)]
                if not waveform_files:
                    continue
                all_good = False
                for filename in waveform_files:
                    files_to_be_deleted.append(filename)
                    self._add_report(
                        "WARNING: "
                        "The event-station raypath for the file\n\t'{f}'\n "
                        "does not fully lay within the domain. You might want "
                        "to remove the file or change the domain "
                        "specifications.".format(f=os.path.relpath(filename)))
        if all_good:
            self._print_ok_message()
        else:
            self._print_fail_message()

        return files_to_be_deleted

    def _validate_station_files_availability(self):
        """
        Checks that all waveform files have an associated station file.
        """
        print(("Confirming that station metainformation files exist for "
               "all waveforms "), end=' ')

        all_good = True

        # Loop over all events.
        for event_name in self.comm.events.list():
            self._flush_point()
            # Get all waveform files for the current event.
            try:
                waveform_info = self.comm.waveforms.get_metadata_raw(
                    event_name)
            except LASIFNotFoundError:
                # If there are none, skip.
                continue
            # Now loop over all channels.
            for channel in waveform_info:
                if self.comm.stations.has_channel(channel["channel_id"],
                                                  channel["starttime"]):
                    continue
                self._add_report(
                    "WARNING: "
                    "No station metainformation available for the waveform "
                    "file\n\t'{waveform_file}'\n"
                    "If you have a station file for that channel make sure "
                    "it actually covers the time span of the data.\n"
                    "Otherwise contact the developers...".format(
                        waveform_file=os.path.relpath(channel["filename"])))
                all_good = False
        if all_good:
            self._print_ok_message()
        else:
            self._print_fail_message()

    def _validate_waveform_files(self):
        """
        Makes sure all waveform files are acceptable.

        It checks that:

        * each station only has data from one location for each event.
        """
        print("Checking all waveform files ", end=' ')

        all_good = True

        # Get all iterations.
        iterations = self.comm.iterations.list()
        # Also add the 'raw' iteration -> kind of a hack to ease the later
        # flow. I'm too lazy to implement something fancy now.
        iterations.append('__RAW__')

        # Loop over all events.
        for event_name in self.comm.events.list():
            self._flush_point()
            print(event_name)
            print('Parsing the iterations ...')

            for iteration in iterations:
                if iteration == "__RAW__":
                    try:
                        waveform_info = self.comm.waveforms.get_metadata_raw(
                            event_name)
                    except LASIFNotFoundError:
                        data_path = \
                            self.comm.waveforms.get_waveform_folder(event_name, "raw")
                        print("NO RAW !!!")
                        if os.path.exists(data_path):
                            try:
                                print(('Deleting directory %s'%data_path))
                                os.rmdir(data_path)
                            except OSError as e:
                                print(("Error: %s : %s" % (data_path, e.strerror)))      
                        continue
                else:
                    processing_tag = self.comm.iterations.get(iteration).processing_tag
                    iteration_tag = self.comm.iterations.get(iteration).long_name
                    try:
                        waveform_info = \
                            self.comm.waveforms.get_metadata_processed(
                                event_name=event_name,tag=processing_tag)
                        # Sort by network, station, and component.
                        info = collections.defaultdict(list)
                        for i in waveform_info:
                            info[(i["network"], i["station"],
                                  i["channel"][-1].upper())].append(i["filename"])

                        s_keys = sorted(info.keys())

                        for key in s_keys:
                            value = info[key]
                            if len(value) == 1:
                                continue
                            all_good = False
                            msg = (
                                "Data waveform files for {it_or_raw} for event '{event}' "
                                "and station '{station}' have {count} waveform files "
                                "for component {component}. Please make sure there is "
                                "only 1 file per component! Offending files:"
                                "\n\t{files}").format(
                                it_or_raw="iteration '%s'" % iteration
                                if iteration != "__RAW__" else "the raw waveforms",
                                station=".".join(key[:2]),
                                event=event_name,
                                count=len(value),
                                component=key[-1],
                                files="\n\t".join(["'%s'" % _i for _i in value]))
                            self._add_report(msg)

                    except LASIFNotFoundError:
                        data_path = \
                            self.comm.waveforms.get_waveform_folder(
                                event_name, "processed", processing_tag)
                        print(("NO PROCESSED %s!!!" % iteration))
                        '''
                        if os.path.exists(data_path):
                            try:
                                print(('Deleting directory %s'%data_path))
                                os.rmdir(data_path)
                                os.remove("%_cache.sqlite"%data_path)
                            except OSError as e:
                                print(("Error: %s : %s" % (data_path, e.strerror)))
                        '''
                    try: 
                        synthetic_info = \
                            self.comm.waveforms.get_metadata_synthetic(event_name, 
                                                                       iteration_tag)
                        # Sort by network, station, and component.
                        info = collections.defaultdict(list)
                        for i in synthetic_info:
                            info[(i["network"], i["station"],
                                  i["channel"][-1].upper())].append(i["filename"])

                        s_keys = sorted(info.keys())

                        for key in s_keys:
                            value = info[key]
                            if len(value) == 1:
                                continue
                            all_good = False
                            msg = (
                                "Synthetic waveform files for {it_or_raw} for event '{event}' "
                                "and station '{station}' have {count} waveform files "
                                "for component {component}. Please make sure there is "
                                "only 1 file per component! Offending files:"
                                "\n\t{files}").format(
                                it_or_raw="iteration '%s'" % iteration
                                if iteration != "__RAW__" else "the raw waveforms",
                                station=".".join(key[:2]),
                                event=event_name,
                                count=len(value),
                                component=key[-1],
                                files="\n\t".join(["'%s'" % _i for _i in value]))
                            self._add_report(msg)

                    except LASIFNotFoundError:
                        data_path = \
                            self.comm.waveforms.get_waveform_folder(
                                event_name, "synthetic", 
                                iteration_tag)
                        print(("NO SYNTHETIC %s!!!" % iteration))
                        '''
                        if os.path.exists(data_path):
                            try:
                                print(('Deleting directory %s'%data_path))
                                os.rmdir(data_path)
                                os.remove("%_cache.sqlite"%data_path)
                            except OSError as e:
                                print(("Error: %s : %s" % (data_path, e.strerror)))
                        stf_path = data_path.replace('SYNTHETICS','STF')
                        if os.path.exists(stf_path):
                            try:
                                print(('Deleting directory %s'%stf_path))
                                os.rmdir(stf_path)
                            except OSError as e:
                                print(("Error: %s : %s" % (data_path, e.strerror)))
                        '''

                        continue



        if all_good:
            self._print_ok_message()
        else:
            self._print_fail_message()

    def _validate_event_files(self, config):
        """
        Validates all event files in the currently active project.

        The following tasks are performed:
            * Validate against QuakeML 1.2 scheme.
            * Check for duplicate ids amongst all QuakeML files.
            * Make sure they contain at least one origin, magnitude and focal
              mechanism object.
            * Some simply sanity checks so that the event depth is reasonable
              and the moment tensor values as well. This is rather fragile and
              mainly intended to detect values specified in wrong units.
            * Events that are too close in time. Events that are less then one
              hour apart can in general not be used for adjoint tomography.
              This will naturally also detect duplicate events.
        """
        import collections
        import itertools
        import math
        from obspy import read_events
        from obspy.io.quakeml.core import _validate as validate_quakeml
        from lxml import etree

        print("Validating %i event files ..." % self.comm.events.count())

        # Start with the schema validation.
        print("\tValidating against QuakeML 1.2 schema ", end=' ')
        all_valid = True
        for event in list(self.comm.events.get_all_events().values()):
            filename = event["filename"]
            self._flush_point()
            if validate_quakeml(filename) is not True:
                all_valid = False
                msg = (
                    "ERROR: "
                    "The QuakeML file '{basename}' did not validate against "
                    "the QuakeML 1.2 schema. Unfortunately the error messages "
                    "delivered by lxml are not useful at all. To get useful "
                    "error messages make sure jing is installed "
                    "('brew install jing' (OSX) or "
                    "'sudo apt-get install jing' (Debian/Ubuntu)) and "
                    "execute the following command:\n\n"
                    "\tjing http://quake.ethz.ch/schema/rng/QuakeML-1.2.rng "
                    "{filename}\n\n"
                    "Alternatively you could also use the "
                    "'lasif add_spud_event' command to redownload the event "
                    "if it is in the GCMT "
                    "catalog.\n\n").format(
                    basename=os.path.basename(filename),
                    filename=os.path.relpath(filename))
                self._add_report(msg)
        if all_valid is True:
            self._print_ok_message()
        else:
            self._print_fail_message()

        # Now check for duplicate public IDs.
        print("\tChecking for duplicate public IDs ", end=' ')
        ids = collections.defaultdict(list)
        for event in list(self.comm.events.get_all_events().values()):
            filename = event["filename"]
            self._flush_point()
            # Now walk all files and collect all public ids. Each should be
            # unique!
            with open(filename, "rb") as fh:
                for event, elem in etree.iterparse(fh, events=("start",)):
                    if "publicID" not in list(elem.keys()) or \
                            elem.tag.endswith("eventParameters"):
                        continue
                    ids[elem.get("publicID")].append(filename)
        ids = {key: list(set(value)) for (key, value) in ids.items()
               if len(value) > 1}
        if not ids:
            self._print_ok_message()
        else:
            self._print_fail_message()
            self._add_report(
                "Found the following duplicate publicIDs:\n" +
                "\n".join(["\t%s in files: %s" % (
                    id_string,
                    ", ".join([os.path.basename(i) for i in faulty_files]))
                    for id_string, faulty_files in ids.items()]),
                error_count=len(ids))

        def print_warning(filename, message):
            self._add_report("WARNING: File '{event_name}' "
                             "contains {msg}.\n".format(
                                 event_name=os.path.basename(filename),
                                 msg=message))

        # Performing simple sanity checks.
        print("\tPerforming some basic sanity checks ", end=' ')
        all_good = True
        for event in list(self.comm.events.get_all_events().values()):
            filename = event["filename"]
            self._flush_point()
            cat = read_events(filename)
            filename = os.path.basename(filename)
            # Check that all files contain exactly one event!
            if len(cat) != 1:
                all_good = False
                print_warning(filename, "%i events instead of only one." %
                              len(cat))
            event = cat[0]

            # Sanity checks related to the origin.
            if not event.origins:
                all_good = False
                print_warning(filename, "no origin")
                continue
            origin = event.preferred_origin() or event.origins[0]
            if (origin.depth % 100.0):
                all_good = False
                print_warning(
                    filename, "a depth of %.1f meters. This kind of accuracy "
                              "seems unrealistic. The depth in the QuakeML "
                              "file has to be specified in meters. Checking "
                              "all other QuakeML files for the correct units "
                              "might be a good idea"
                    % origin.depth)
            if (origin.depth > (800.0 * 1000.0)):
                all_good = False
                print_warning(filename, "a depth of more than 800 km. This is"
                                        " likely wrong.")

            # Sanity checks related to the magnitude.
            if not event.magnitudes:
                all_good = False
                print_warning(filename, "no magnitude")
                continue

            # Sanity checks related to the focal mechanism.
            if not event.focal_mechanisms:
                all_good = False
                print_warning(filename, "no focal mechanism")
                continue

            focmec = event.preferred_focal_mechanism() or \
                event.focal_mechanisms[0]
            if not hasattr(focmec, "moment_tensor") or \
                    not focmec.moment_tensor:
                all_good = False
                print_warning(filename, "no moment tensor")
                continue

            mt = focmec.moment_tensor
            if not hasattr(mt, "tensor") or \
                    not mt.tensor:
                all_good = False
                print_warning(filename, "no actual moment tensor")
                continue
            tensor = mt.tensor

            # Convert the moment tensor to a magnitude and see if it is
            # reasonable.
            mag_in_file = event.preferred_magnitude() or event.magnitudes[0]
            mag_in_file = mag_in_file.mag
            M_0 = 1.0 / math.sqrt(2.0) * math.sqrt(
                tensor.m_rr ** 2 + tensor.m_tt ** 2 + tensor.m_pp ** 2)
            magnitude = 2.0 / 3.0 * math.log10(M_0) - 6.0
            # Use some buffer to account for different magnitudes.
            if not (mag_in_file - 1.0) < magnitude < (mag_in_file + 1.0):
                all_good = False
                print_warning(
                    filename, "a moment tensor that would result in a moment "
                              "magnitude of %.2f. The magnitude specified in "
                              "the file is %.2f. Please check that all "
                              "components of the tensor are in Newton * meter"
                    % (magnitude, mag_in_file))

        if all_good is True:
            self._print_ok_message()
        else:
            self._print_fail_message()

        # Collect event times
        event_infos = list(self.comm.events.get_all_events().values())

        # Now check the time distribution of events.
        print("\tChecking for duplicates and events too close in time %s" % \
              (self.comm.events.count() * "."), end=' ')
        all_good = True
        # Sort the events by time.
        event_infos = sorted(event_infos, key=lambda x: x["origin_time"])
        # Loop over adjacent indices.
        a, b = itertools.tee(event_infos)
        next(b, None)
        for event_1, event_2 in zip(a, b):
            time_diff = abs(event_2["origin_time"] - event_1["origin_time"])
            # If time difference is under one hour, it could be either a
            # duplicate event or interfering events.
            if time_diff <= 3600.0:
                all_good = False
                self._add_report(
                    "WARNING: "
                    "The time difference between events '{file_1}' and "
                    "'{file_2}' is only {diff:.1f} minutes. This could "
                    "be either due to a duplicate event or events that have "
                    "interfering waveforms.\n".format(
                        file_1=event_1["filename"],
                        file_2=event_2["filename"],
                        diff=time_diff / 60.0))
        if all_good is True:
            self._print_ok_message()
        else:
            self._print_fail_message()

        # Check that all events fall within the chosen boundaries.
        if config == "teleseismic":
            print((\
                "\tSkipping the check that all %s events are in chosen domain for teleseismic configuration"%\
                self.comm.events.count()))
            all_good = True
        else:
            print("\tAssure all events are in chosen domain %s" % \
                  (self.comm.events.count() * "."), end=' ')
            all_good = True
            domain = self.comm.project.domain
            for event in event_infos:
                if domain.point_in_domain(latitude=event["latitude"],
                                          longitude=event["longitude"]):
                    continue
                all_good = False
                self._add_report(
                    "\nWARNING: "
                    "Event '{filename}' is out of bounds of the chosen domain."
                    "\n".format(filename=event["filename"]))
        if all_good is True:
            self._print_ok_message()
        else:
            self._print_fail_message()

    def is_event_station_raypath_within_boundaries(
            self, event_name, station_latitude, station_longitude,
            raypath_steps=25):
        """
        Checks if the full station-event raypath is within the project's domain
        boundaries.

        Returns True if this is the case, False if not.

        :type event_latitude: float
        :param event_latitude: The event latitude.
        :type event_longitude: float
        :param event_longitude: The event longitude.
        :type station_latitude: float
        :param station_latitude: The station latitude.
        :type station_longitude: float
        :param station_longitude: The station longitude.
        :type raypath_steps: int
        :param raypath_steps: The number of discrete points along the raypath
            that will be checked. Optional.
        """
        from lasif.utils import greatcircle_points, Point
        import lasif.domain

        ev = self.comm.events.get(event_name)

        domain = self.comm.project.domain

        # Shortcircuit.
        if isinstance(domain, lasif.domain.GlobalDomain):
            return True

        for point in greatcircle_points(
                Point(station_latitude, station_longitude),
                Point(ev["latitude"], ev["longitude"]),
                max_npts=raypath_steps):

            if not domain.point_in_domain(latitude=point.lat,
                                          longitude=point.lng):
                return False
        return True

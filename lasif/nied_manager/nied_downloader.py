#!/usr/bin/env python
# -*- coding: utf-8 -*-

# NIED version function of MassDownloader (obspy.client.fdsn.mass_downloader)
# used at
# lasif.components.downloads DownloadsComponent.donwload_data_for_one_evnt

import os, datetime, getpass, glob, shutil
import obspy
from obspy.clients.fdsn.mass_downloader import GlobalDomain
from obspy import UTCDateTime

from HinetPy import Client, win32


class NiedDownloader():
    """
    download manager for seismic data provided by NIED
    To use this class, login account for accessing NIED server is required.
    win32tools is also required for converting win32 format to SAC PZ.

    network name for NIED should be "NIED.0101" (NIED. + station code)
    """

    def __init__(self,starttime=None,endtime=None,network=None,
                      domain=None):
        self.starttime    = starttime
        self.endtime      = endtime
        self.network      = network.split(",")
        self.domain       = domain
        self.outdir       = "./test_nied_lasif" # debug

        # time values in JST
        self.starttime_JST = self.starttime + datetime.timedelta(hours=9)
        self.endtime_JST   = self.endtime + datetime.timedelta(hours=9)

        # authentication for NIED server
        user = input("Username: ")
        pswd = getpass.getpass("Password: ")
        self.client = Client(user,pswd)

        # select stations from domain boundary
        for network in self.network:
            network_code = network.lstrip("NIED.")
            print("downloading NIED network-{} data".format(network_code))
            if self.domain.__class__.__name__ == "SphericalSectionDomain":
                params = self.domain.get_query_parameters()
                self.client.select_stations(network_code,
                                                latitude =params["latitude"],
                                                longitude=params["longitude"],
                                                maxradius=params["maxradius"])
            else:
                self.client.select_stations(network_code, minlatitude=self.domain.min_latitude,
                                                maxlatitude =self.domain.max_latitude,
                                                minlongitude=self.domain.min_longitude,
                                                maxlongitude=self.domain.max_longitude)


    def download(self, outdir=None, stationdir=None):
        self.outdir    =outdir
        self.stationdir=stationdir

        # download waveform
        #for network in self.network:
        #    network_code = network.lstrip("NIED.")
        #    output_dir   = os.path.join(self.outdir,network_code)
        #    span         = int((self.endtime_JST-self.starttime_JST)/60.0) # in minutes
        #    data, ctable = self.client.get_continuous_waveform(network_code,self.starttime_JST.datetime,span,outdir=output_dir)

        ## convert waveform file format
        #self.convert_format()

        # make symlinks of SACPZ files in Station/SACPZ
        sacs = glob.glob(self.outdir+"/*.SAC")
        pzs  = glob.glob(self.outdir+"/*.SAC_PZ")

        for sac in sacs:
            try: # make symlink if not exits
                # here event name is added at the head of filename
                # to be found by station_cache registration process.
                link_name=self.outdir.split("/")[-2]+"_"+os.path.basename(sac)
                print(link_name)
                os.symlink(sac, os.path.join(stationdir,link_name))
            except:
                pass
        # pz files are replaced into the station files
        for pz in pzs:
            try:
                shutil.move(pz, os.path.join(stationdir,os.path.basename(pz)))
            except:
                pass

    def convert_format(self):
        network_dirs = [ f.path for f in os.scandir(self.outdir) if f.is_dir() ]

        for network_dir in network_dirs:
            datas   = sorted(glob.glob(network_dir+"/*.cnt"))
            ctables = sorted(glob.glob(network_dir+"/*.ch"))

            for j in range(len(datas)):
                data   = datas[j]
                ctable = ctables[j]

                # extract wave
                win32.extract_sac(data,ctable,outdir=self.outdir)
                # extract instrumental responses
                win32.extract_pz(ctable,outdir=self.outdir)

            # convert SACPZ to mseed for LASIF's requirement
            #sacs = glob.glob(self.outdir+"/*.SAC")
            #for sac in sacs:
            #    try:
            #        st = obspy.read(sac)
            #        st.write(sac.rstrip(".SAC")+".mseed",format="MSEED")
            #    except:
            #        print("{} skipped. some problem on this file.".format(sac))

            # erase the win32 format datafile
            shutil.rmtree(network_dir)
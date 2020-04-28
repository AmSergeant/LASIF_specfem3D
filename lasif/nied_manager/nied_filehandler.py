#!/usr/bin/env python
# -*- coding: utf-8 -*-
import obspy, datetime

def simple_sacpz_parser(filename):
    """
    return a list of channels' dict for a given filename

    one dict must includes:
        channel_id,start_date,end_date,
        latitude,longitude,elevation_in_m,local_depth_in_m

    The datetime value in NIED SAC files are all in JST time (UTC + 9).
    """

    channel = {}
    st = obspy.read(filename)
    stat = st[0].stats
    sac = st[0].stats.sac
    channel["channel_id"]       = "{}.{}.{}.{}".format(stat["network"],stat["station"],stat["location"],stat["channel"])
    channel["start_date"]       = stat["starttime"]+datetime.timedelta(hours=-9)
    channel["end_date"]         = stat["endtime"]  +datetime.timedelta(hours=-9)
    channel["latitude"]         = sac.stla
    channel["longitude"]        = sac.stlo
    channel["elevation_in_m"]   = sac.stel
    channel["local_depth_in_m"] = -1.*sac.stel

    return [channel]
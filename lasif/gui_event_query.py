#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 18:19:38 2020

Project specific function for a gui app to select the evets in the query 
and store them in the project.

:copyright:
    Amandine Sergeant (sergeant@lma.cnrs-mrs.fr), 2019
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""

# useful tutorials:
# https://www.learnpyqt.com/courses/start/creating-your-first-window/
# check here: https://stackoverflow.com/questions/8963082/how-to-update-a-matplotlib-figure-from-a-function-in-a-qt-gui
# http://zetcode.com/gui/pyqt5/firstprograms/
# https://pythonspot.com/pyqt5-matplotlib/


import os
import sys

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from PyQt5 import QtCore, QtWidgets, QtGui
from PyQt5.QtCore import pyqtSlot

#from PyQt5.QtWidgets import *
#from PyQt5.QtCore import *
#from PyQt5.QtGui import *

from lasif.components.project import Project
from lasif.utils import get_event_filename
import lasif.visualization

from obspy.core.event import Catalog
from obspy.imaging.beachball import beach
from obspy.geodetics.base import gps2dist_azimuth, locations2degrees
from obspy.taup import TauPyModel

import numpy as np


class PlotCanvas(FigureCanvas):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.figure = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.figure.add_subplot(111)

        FigureCanvas.__init__(self, self.figure)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                QtWidgets.QSizePolicy.Expanding,
                QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        
        
        
class MainWindow(QtWidgets.QMainWindow):
    
    def __init__(self, comm, chosen_events, events):
        super().__init__()
        self.left = 10
        self.top = 10
        self.title = 'Select suitable events in the filtered catalog'
        self.width = 800
        self.height = 600
        
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
          
        self.init_figure(comm)
        self.init_buttons()
        self.init_shortcuts()
        self.ui(comm, chosen_events, events)
        
        
    def init_figure(self, comm):
        # Create the maptlotlib FigureCanvas object, 
        # which defines a single set of axes as self.axes.
        self.canvas = PlotCanvas(self, width=5, height=4)
        self.setCentralWidget(self.canvas)
               
        # Define axes
        # for global map
        self.map_ax = self.canvas.axes
        # for the current beachball
        self.beachball_ax = self.canvas.figure.add_axes([0.048, 0.65, 0.3, 0.25])
        # for seismic arrivals on current beachball
        self.beachball_subax = self.canvas.figure.add_axes(self.beachball_ax.get_position(), 
                                  projection='polar', label='pol', frameon=False)
        self.beachball_subax.set_aspect('equal')
        self.beachball_subax.set_theta_direction(-1)
        self.beachball_subax.set_theta_zero_location("N")
        
        
    def init_buttons(self):
        # select all events and quit button
        self.quit_btn = QtWidgets.QPushButton('Select all remaining events\n and Quit', self)
        self.quit_btn.setStyleSheet('QPushButton {color: red;}')
        self.quit_btn.resize(self.quit_btn.sizeHint())
        self.quit_btn.move(580, 500)
        
        # go to next event button
        self.next_btn = QtWidgets.QPushButton('Go to next event\n or press \'right arrow\'', self)
        self.next_btn.resize(self.next_btn.sizeHint())
        self.next_btn.move(580, 550)
        
        # select/disregard the event buttons
        self.select_btn = QtWidgets.QCheckBox(self)
        self.select_btn.setText('select the event\n or press \'s\'')
        self.select_btn.resize(self.select_btn.sizeHint())
        self.select_btn.move(430, 550)
        self.unselect_btn = QtWidgets.QCheckBox(self)
        self.unselect_btn.setText('disregard the event\n or press \'q\'')
        self.unselect_btn.resize(self.unselect_btn.sizeHint())
        self.unselect_btn.move(280, 550)
        self.unselect_btn.setChecked(False)
        self.select_btn.setChecked(True)
        
    def init_shortcuts(self):
        # shortcut "right arrow" for "go to next event" button
        self.next_shortcut = QtWidgets.QShortcut(QtGui.QKeySequence("right"), self.next_btn)
        self.next_shortcut.setEnabled(True)
        
        # shortcut "q" for unselect/disregard the event
        self.unselect_shortcut = QtWidgets.QShortcut(QtGui.QKeySequence("q"), self.unselect_btn)
        self.unselect_shortcut.setEnabled(True)
        # shortcut "s" for select the event
        self.select_shortcut = QtWidgets.QShortcut(QtGui.QKeySequence("s"), self.select_btn)
        self.select_shortcut.setEnabled(True)
        
        
                  

    def ui(self, comm, chosen_events, events):
        # Get infos from Lasif project
        self.comm = comm
        self._gather_infos_from_lasif_comm()
        # earth model for taup
        self.earth_model = TauPyModel("ak135")
        # store the nb of events already existing in the project database
        self.nb_of_events_in_database = len(list(comm.events.get_all_events().values()))
        
        
        # store the events to inspect
        self.events = chosen_events
        self.nb_of_events_to_inspect = len(self.events)
        
        # initalize the final event catalog that will be selected
        self.final_event_cat = []
        
        # plot event global map, draw just once
        self.basemap = comm.visualizations.plot_events(events=events,
            ax=self.map_ax, azimuthal_projection=self.azimuthal_proj)
               
        # launch with first event
        self.index_current_event = 0
        self.current_event = self.events[self.index_current_event]
        # initiate the current beachball patch
        self.current_mt_patches = []
        # do plots
        self._update_event_map() 
        
        
        # connect the "go to next event" button to action
        self.next_btn.clicked.connect(self._go_to_next_event)
        self.next_shortcut.activated.connect(self._go_to_next_event)
        # connect the "select all events and quit" button to action
        self.quit_btn.clicked.connect(self._select_all_remaining_events_and_quit)

        # check status of "select event" buttons
        self.unselect_btn.stateChanged.connect(self.unselect_checkBoxChangedAction)
        self.select_btn.stateChanged.connect(self.select_checkBoxChangedAction)
        self.select_shortcut.activated.connect(self.click_on_select_button)
        self.unselect_shortcut.activated.connect(self.click_on_unselect_button)
        
        # show
        self.show() 
       
        
    def _draw(self):
        self.canvas.draw()
        
        
    def add_selected_phases_on_beachball(self):
        # compute ray parameters and plot them on beachball
        def compute_ray_parameters(self):
            azimuths = []
            takeoff = []
            for point in self.domain_points:
                station = self.domain_points[point]
                epicentral_distance, azimuth, baz = gps2dist_azimuth(
                    self.event_info["latitude"], self.event_info["longitude"], 
                    station["latitude"],station["longitude"])
                dist_in_deg = locations2degrees(
                    self.event_info["latitude"], self.event_info["longitude"], 
                    station["latitude"],station["longitude"])
                tts = self.earth_model.get_travel_times(source_depth_in_km=self.event_info["depth_in_km"],
                                                        distance_in_degree=dist_in_deg,
                                                        phase_list=self.Phase)
                if tts:
                    azimuths.append(baz)
                    takeoff.append(tts[0].takeoff_angle)
            return azimuths, takeoff
        
        # plot
        azimuths, takeoff = compute_ray_parameters(self)
        for azimuth, takeoff_angle in zip(azimuths,takeoff):
            self.beachball_subax.plot(azimuth, takeoff_angle, '+', ms=10., mew=2.0, mec='black') #withdash=True
            self.beachball_subax.set_rmax(np.pi / 2.)
            self.beachball_subax.set_yticks([0, np.pi/6., 2.*np.pi/6., np.pi/2.])
            self.beachball_subax.set_yticklabels([])
            self.beachball_subax.set_xticklabels([])
        self.beachball_subax.set_axis_off()
            
        
    def _update_event_map(self):
        # clear beachball axes
        self.beachball_ax.cla()
        self.beachball_subax.cla()
        # clear current event mt on global map
        for i in self.current_mt_patches:
            i.remove()
        
        self.event_info = self.current_event["event_info"]
        
        # draw current event mt on global map
        self.current_mt_patches = lasif.visualization.plot_events(
            events=[self.event_info], map_object=self.basemap,
            beachball_size=self.beachballsize, color="blue")
        
        # draw current event mt and arrivals on beachball axes
        focmec = [self.event_info["m_rr"], self.event_info["m_tt"], 
                  self.event_info["m_pp"], self.event_info["m_rt"],
                  self.event_info["m_rp"], self.event_info["m_tp"]]
        mybeachball = beach(focmec, linewidth=1, xy=(0.5,0.5), width=1, facecolor='b')
        self.beachball_ax.add_collection(mybeachball)
        self.beachball_ax.set_aspect('equal')
        self.beachball_ax.set_axis_off()
        self.add_selected_phases_on_beachball()
        
    
        # titles
        self.canvas.figure.suptitle("Event %i/%i to inspect\n %s" \
                                     %(self.index_current_event+1, 
                                       self.nb_of_events_to_inspect,
                                       self.event_info["event_name"]), fontsize=10)
        self.map_ax.set_title("%i events in the project database"
                                 %self.nb_of_events_in_database, 
                                 fontsize=8, fontweight = "bold")
        self._draw()
        
        
       
    def _gather_infos_from_lasif_comm(self):
        
        # get project infos
        proj = self.comm.project
        
        # phase of interest to plot
        self.Phase = [proj.config["download_settings"]["phase_of_interest"]]
        
        # get the domain corners and center to plot seismic arrivals on beachball
        self.domain_points = self.comm.query.domain_corners()
        central_point = self.comm.query.center()
        self.domain_points["center"] = \
            {"latitude": central_point.latitude,
             "longitude": central_point.longitude}
        
        self.azimuthal_proj = False
        self.beachballsize = 0.02
        if proj.config["download_settings"]["configuration"] == "teleseismic":
            self.azimuthal_proj = True
            self.beachballsize = 0.04
            
        
     
    @pyqtSlot()
    def _check_and_store_event_status(self):
        
        # if select button is True: store the event and upgrade nb of saved events
        if self.select_btn.isChecked() == True :
            self.unselect_btn.setChecked(False)
            
            print("Will save the current event: %s"%self.event_info["event_name"])
            self.final_event_cat.append(self.current_event)
            self.nb_of_events_in_database += 1
            
            # upgrade the label for save events
            #self._label_for_events_in_database()
        
        elif self.unselect_btn.isChecked() == True :
            self.select_btn.setChecked(False)
            print("Will remove the current event: %s"%self.event_info["event_name"])
            
            
            
    @pyqtSlot()    
    def _update_current_event(self):
        index_next_event = self.index_current_event +1
        if index_next_event == self.nb_of_events_to_inspect :
            # we reached the end, now quit
            print("You have inspected all events, now quitting")
            print(self.final_event_cat)
            self._write_catalog_and_quit_the_app()
        else:
            self.index_current_event = index_next_event
            self.current_event = self.events[self.index_current_event]
            
            
    
    @pyqtSlot()
    def _go_to_next_event(self):
        self._check_and_store_event_status()
        self._update_current_event()
        self.close()
        self._update_event_map()
        # set the select buttons back to defaults values
        self.unselect_btn.setChecked(False)
        self.select_btn.setChecked(True) 
        self.show()
            
        
    @pyqtSlot()
    def _select_all_remaining_events_and_quit(self):
        msg = ("--> Will select all remaining events and quit the GUI")
        print(msg)
        if self.unselect_btn.isChecked() == True :
            first_event_to_select = self.index_current_event + 1
        else:
            first_event_to_select = self.index_current_event 
            
        for index_event in np.arange(first_event_to_select, 
                                      self.nb_of_events_to_inspect,1):
            self.final_event_cat.append(self.events[index_event])
            
        self._write_catalog_and_quit_the_app()
    
    
    @pyqtSlot()
    def _write_catalog_and_quit_the_app(self):
        
        print("You have finally selected %i events"%len(self.final_event_cat))
        folder = self.comm.project.paths["events"]
        for item in self.final_event_cat:
            filename = os.path.join(folder, get_event_filename(item, "GCMT"))
            Catalog(events=[item]).write(filename, format="quakeml",
                                          validate=True)
            print(("Written %s" % (os.path.relpath(filename))))
        
        # Quit the gui
        QtCore.QCoreApplication.instance().quit()
        
        
    
    
    def select_checkBoxChangedAction(self, state):
        if (QtCore.Qt.Checked == state):
                self.unselect_btn.setChecked(False)
                
    def unselect_checkBoxChangedAction(self, state):
        if (QtCore.Qt.Checked == state):
                self.select_btn.setChecked(False)
                
    def click_on_unselect_button(self):
        self.select_btn.setChecked(False)
        self.unselect_btn.setChecked(True)
        
    def click_on_select_button(self):
        self.unselect_btn.setChecked(False)
        self.select_btn.setChecked(True)
                
                
        
    
        
        
        

def launch_event_gui(chosen_events, events):
    
    comm = Project('.').get_communicator()
        
    # Launch and open the window.
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow(comm, chosen_events, events)
    #app.exec_()
    # Move window to center of screen.
    window.move(
        app.desktop().screen().rect().center() - window.rect().center())
    # Show and bring window to foreground.
    window.show()
    window.raise_()
    sys.exit(app.exec_())



    
    
    







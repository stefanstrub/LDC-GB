# -*- coding: utf-8 -*-

"""
Graphs for including GW source in the default LISANode instrument. 

"""
from .naming import lisa_links, lisa_indices, indexed, fromto
from .propag import SingleLinkPropagation, LaserLinks
from .instrument import LISA, Spacecraft
from .tdi import LISAWithTDI, TDI
from .config import LISA_PUBLISH_MEASUREMENTS
from .config import LISA_PUBLISH_BEATNOTE_FREQUENCIES
from .config import LISA_MEASUREMENT_FS, LISA_PHYSICS_FS, LISA_TDI_FS
from ..compiler import Graph

from ldc.lisa.projection import from_file
import h5py
import os
import yaml
from scipy.constants import c as clight


DEFAULT_CONFIG = "config.yml"

def gw_infos(config):
    """ Return GW infos from config file
    """
    d = yaml.load(open(config, "r"), Loader=yaml.BaseLoader)
    d["source_signal"] = d["gwstrain_file"]
    yArm, source_names, links, t_min, t_max, dt = from_file(d["source_signal"],
                                                            nodata=True)
    d["upsampling"] = int(LISA_PHYSICS_FS/(1./dt))
    d["links"] = links
    d["dt"] = dt
    return d

def orbits_infos(config):
    """ Return Orbits infos from config file
    """
    d = yaml.load(open(config, "r"), Loader=yaml.BaseLoader)
    d["t_min"] = float(d["t_min"])
    d["t_max"] = float(d["t_max"])
    d["dt"] = float(d["dt"]) 
    d["nominal_arm_length"] = float(d["nominal_arm_length"]) 
    d["init_rotation"] = float(d["initial_rotation"]) 
    d["init_position"] = float(d["initial_position"])
    d["travel_time_order"] = int(d["travel_time_order"]) 
    return d

def noise_infos(config):
    """ Return Noise infos from config file
    """
    d = yaml.load(open(config, "r"), Loader=yaml.BaseLoader)
    d["on_off"] = int(d["noise"])
    return d
    
class LISAWithGWAndTDI(LISAWithTDI):
    def __init__(self, **kwargs):
        Graph.__init__(self, "LISAWithGWAndTDI")

        if "config" in kwargs:
            self.config = kwargs["config"]
        else:
            self.config = DEFAULT_CONFIG
            
        self.add(LISAWithGW, name="lisa", measurements=True, 
                 beatnotes=True, config=self.config)
        self.add(TDI, name="tdi")
        
        self.connect_lisa_to_tdi()
        self.publish_all_params()
        self.publish_laser_offsets()
        self.publish_tdi_variables()

        tt_mean_value = 8 #sec
        tdi_order = min(((tt_mean_value * 2) - 1) * LISA_MEASUREMENT_FS, 31)
        self.nodes['tdi'].params = {'eta_order': tdi_order}

        O = orbits_infos(self.config)
        self.params = {"armlength":O["nominal_arm_length"],
                       "init_rotation":O["init_rotation"],
                       "init_position":O["init_position"],
                       "ldc_tt_order":O["travel_time_order"]}
        self.duration = O["t_max"] + 10
        
    def publish_tdi_variables(self):
        """Publish TDI variables."""
        
        for t in ["X","Y","Z"]:
            avg =  "avg_" + t
            self.add("EllipticFilter<double>", avg)
            dt_target = LISA_TDI_FS
            upsampling = int(LISA_PHYSICS_FS/(1./dt_target))
            self.nodes[avg].params = {
                'passband_freq': 0.4 * (1/3.45), 
                'stopband_freq': 0.48 * (1/3.45), 
                'minimum_passband_gain': 0.1,
                'minimum_attenuation': 100}
            decim = "decimation_" + t
            self.add("Decimation<double>", decim)
            self.nodes[decim].downsampling = upsampling
            self.nodes[decim].lag = 10 # group delay found empirically
            self.connect("tdi."+t, avg + ".input")
            self.connect(avg+".result", decim + ".input")
            self.publish_output(decim + ".result", t)
            
class LISAWithGW(LISA):
    """LISA constellation. 

    See LISA object documentation in LISANode for more details.
    """
    def __init__(self,
                 measurements=LISA_PUBLISH_MEASUREMENTS,
                 beatnotes=LISA_PUBLISH_BEATNOTE_FREQUENCIES,
                 config=DEFAULT_CONFIG):
        """Initialize the LISA instrument.

        Args:
            measurement: whether interferometric measurements should be published
            travel_times: whether ranging information should be published
            beatnotes: whether beatnote frequencies should be published
        """
        # pylint: disable=R0914,R0915,R0912
        #super().__init__("LISAWithGW")
        Graph.__init__(self, "LISAWithGW")

        # Store arguments
        self.measurements = measurements
        self.beatnotes = beatnotes
        self.config = config
        
        # Add three spacecraft
        self.build_spacecraft(1)
        self.build_spacecraft(2)
        self.build_spacecraft(3)

        # Propagate signals
        SIGNALS = {"laser", "laser_frequency", "dt"}
        self.propagate_signals()#signals)
        SingleLinkPropagation.SIGNALS = SIGNALS
        self.publish_orbit_parameters()

        N = noise_infos(config)
        dnoise = {"lasernoise_on_off":0, "accelnoise_on_off":0, "readoutnoise_on_off":0,
                  "opticalnoise_on_off":0, "usonoise_on_off":0, "modulationnoise_on_off":0,
                  "rangingnoise_on_off":0}

        if N["on_off"]:
            dnoise["accelnoise_on_off"] = 1
            dnoise["readoutnoise_on_off"] = 1
            dnoise["opticalnoise_on_off"] = 1

            if 'opticalnoise' in N.keys():
                dnoise["opticalnoise_sci_gain"] = float(N["opticalnoise"])/clight
            if 'readoutnoise' in N.keys():
                dnoise["readoutnoise_sci_gain"] = float(N["readoutnoise"])/clight

        for i in range(1,4):
            self.nodes["sc{}".format(i)].params = dnoise

        O = orbits_infos(config)
        self.params = {"armlength":O["nominal_arm_length"],
                       "init_rotation":O["init_rotation"],
                       "init_position":O["init_position"],
                       "ldc_tt_order":O["travel_time_order"]}


    def propagate_signals(self):
        """Add and connect a LaserLinks node to the list of signals."""
        # Add a LaserLinks node

        GW = gw_infos(self.config)
        self.add(LaserLinks, name="laser_links", orbit_type="ldc",
                 gw_input=GW["source_signal"]) #signals=signals
        self.nodes['laser_links'].params = {'order': 31}
        # Publish the interpolation order as a parameter
        self.publish_param("laser_links.order", "order")
        # Connect all signals for each link
        for i, j in lisa_links():
            self.propagate_signals_along_link(i, j)




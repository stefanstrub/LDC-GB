# -*- coding: utf-8 -*-

"""
Graphs for including GW source in the default LISANode instrument. 

"""
from .tdi import LISAWithTDI, TDI
from .instrument import LISA
from .config import LISA_MEASUREMENT_FS, LISA_PHYSICS_FS
from .config import TDI_PUBLISH_INTERVAR, LISA_TDI_FS
from .config import LISA_DEBUG_CLOCK_OFFSETS
from ..compiler import Graph

class LISAWithGWAndTDI(LISAWithTDI):

    def __init__(self):
        Graph.__init__(self, "LISAWithGWAndTDI") 
        self.add(LISA, name="lisa")
        self.add(TDI, name="tdi")

        self.connect_lisa_to_tdi()

        self.publish_all_params()
        self.publish_laser_offsets()
        self.publish_tdi_variables()
        self.publish_pprs()
        self.publish_beatnotes()
        self.publish_timer_deviations()
        self.publish_dws()

        if LISA_DEBUG_CLOCK_OFFSETS:
            self.publish_debug_clock_offsets()
        if TDI_PUBLISH_INTERVAR:
            self.publish_intermediary_variables()
        
    def publish_tdi_variables(self):
        """Publish TDI variables."""

        for t in ["X","Y","Z"]:
            self.publish_output("tdi."+t , t)

        
            
    # def publish_tdi_variables(self):
    #     """Publish TDI variables."""
        
    #     for t in ["X","Y","Z"]:
    #         avg =  "avg_" + t
    #         self.add("EllipticFilter<double>", avg)
    #         dt_target = 1/LISA_TDI_FS
    #         upsampling = int(LISA_PHYSICS_FS/(1./dt_target))
    #         self.nodes[avg].params = {
    #             'passband_freq': 0.4 * (1/3.45), 
    #             'stopband_freq': 0.48 * (1/3.45), 
    #             'minimum_passband_gain': 0.1,
    #             'minimum_attenuation': 100}
    #         decim = "decimation_" + t
    #         self.add("Decimation<double>", decim)
    #         self.nodes[decim].downsampling = upsampling
    #         self.nodes[decim].lag = 10 # group delay found empirically
    #         self.connect("tdi."+t, avg + ".input")
    #         self.connect(avg+".result", decim + ".input")
    #         self.publish_output(decim + ".result", t)

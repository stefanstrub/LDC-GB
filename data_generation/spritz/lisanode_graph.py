# -*- coding: utf-8 -*-

"""
Graphs for including GW source in the default LISANode instrument. 

"""
from .tdi import LISAWithTDI, TDI
from .instrument import LISA
from .config import LISA_MEASUREMENT_FS, LISA_PHYSICS_FS
from .config import TDI_PUBLISH_INTERVAR
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

        
            

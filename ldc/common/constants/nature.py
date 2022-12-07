""" 
Physical constants
"""

import lisaconstants
from .supported import supported


class Nature:
    """Temporary class to ensure transition with existing code.

    """
    def __getattr__(self, parname):
        if parname in supported.keys():
            #print(f"Warning: in a near future this class won\'t be supported anymore. Please use lisaconstants.{supported[parname]} instead")
            return getattr(lisaconstants, supported[parname])
        elif parname in self.__dict__.keys():
            # this will throw the right exception if we don't have this attribute
            return self.__dict__[parname] 
        raise AttributeError



    

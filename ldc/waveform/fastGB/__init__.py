from .fastgb_tools import get_default_orbits, get_buffersize 
from fastGB import pyGB as FastGB

try:
    from .fastgb_wrapper import VFastGB    # vectorized
except:
    pass

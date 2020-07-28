"""Provides a suite of I/O routine to load and save parameters in
yml file format.
"""

#import yaml
from astropy import units
from astropy.io.misc import yaml

# def decode(value):
#     """ Convert to numerical value if possible.

#     >>> decode('1.0')
#     1.0
#     """
#     if isinstance(value, list):
#         value = [decode(v) for v in value]
#     elif isinstance(value, dict):
#         value = str(value)
#     else:
#         try:
#             value = units.Quantity(value)
#         except:
#             print(value)
#             try:
#                 value = float(value) if "." in value else int(value)
#             except:
#                 pass
#     return value

def save_config(filename, cfg, name="config"):
    """ Write config to yml file

    >>> save_config("test.yml", {'author':'me', 'date':'today'})
    >>> print(load_config("test.yml"))
    {'author': 'me', 'date': 'today'}
    """
    if name in cfg.keys():
        cfg = cfg[name]
    yaml.dump(cfg, open(filename, "a"), default_flow_style=False)


def load_config(filename, name="config"):
    """ Load config from yml file
    """
    cfg = yaml.load(open(filename, "r"))#, Loader=yaml.BaseLoader)
    if name in cfg.keys():
        cfg = cfg[name]
    #for k, v in cfg.items():
    #    cfg[k] = decode(v)
    return cfg


if __name__ == "__main__":
    import doctest
    doctest.testmod()

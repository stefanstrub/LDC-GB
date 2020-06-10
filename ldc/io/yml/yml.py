""" Yml generic I/O 
"""
import yaml

def decode(value):
    """ Convert to numerical value if possible. 
    """
    try:
        value = float(value) if "." in value else int(value)
    except:
        pass
    return value


class YML:
    """Provides a suite of I/O routine to load and save parameters in
    yml file format.
    """

    def __init__(self, filename):
        """ Initialize YML I/O
        """
        self.filename = filename
        
    def save_config(self, cfg, name="config"):
        """ Write config to yml file
        """
        if name in cfg.keys():
            cfg = cfg[name]
        yaml.dump(cfg, open(self.filename, "a"), default_flow_style=False)
        

    def load_config(self, name="config"):
        """ Load config from yml file
        """
        cfg = yaml.load(open(self.filename, "r"), Loader=yaml.BaseLoader)
        if name in cfg.keys():
            cfg = cfg[name]
        for k,v in cfg.items():
            cfg[k] = decode(v)
        return cfg



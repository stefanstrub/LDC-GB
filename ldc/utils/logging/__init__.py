import logging
import sys

def init_logger(filename=None, level=logging.DEBUG, name='lisa'):
    """ Initialize a logger.
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    if (len(logger.handlers) < 2):
        formatter = logging.Formatter("%(asctime)s - %(name)s - "
                                      "%(levelname)s - %(message)s")
        if filename:
            rfhandler = logging.FileHandler(filename)
            logger.addHandler(rfhandler)
            rfhandler.setFormatter(formatter)
        if level:
            shandler = logging.StreamHandler(sys.stdout)
            shandler.setLevel(level)
            shandler.setFormatter(formatter)
            logger.addHandler(shandler)
    return logger

def close_logger(logger):
    """ Close a logger
    """
    for h in logger.handlers:
        logger.removeHandler(h)
        h.flush()
        h.close()

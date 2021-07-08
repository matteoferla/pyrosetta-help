import logging, sys, io, pyrosetta
import re
from typing import *

def configure_logger() -> logging.Logger:
    """
    The function `get_logger`, simply adds a stringIO handler to the log and captures the log,
    thus making it easier to use.
    The function `get_log_entries`, spits out entries of a given level.

    :return: logger
    """
    pyrosetta.logging_support.set_logging_sink()
    logger = logging.getLogger("rosetta")
    logger.setLevel(logging.INFO) # default = logging.WARNING
    stringio = io.StringIO()
    handler = logging.StreamHandler(stringio)
    handler.setLevel(logging.INFO)
    handler.set_name('stringio')
    handler.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s - %(message)s'))
    logger.addHandler(handler)
    return logger


def get_log_entries(levelname: Union[str, int]=logging.INFO, query: Optional[str]=None) -> List[Dict[str, Any]]:
    """
    Get a list of all entries in log at a given level.
    levelname can be either an int (``logging.INFO`` etc. are numbers multiples of 10 in increasing severity)
    or a string of the level.
    Note that it is very crude: if INFO is requested, ERROR is not shown!

    :param levelname: int for the level number or str of the name
    :return: List of str
    """
    if isinstance(levelname, int):
        # logging.INFO is actually an int, not an enum
        levelnumber = logging.getLevelName(levelname)
        entries = get_all_log_entries()
        if query is None:
            return [e for e in entries if e['level'] == levelnumber]
        else:
            return [e for e in entries if str(query) in e['text']]

def get_all_log_entries() -> List[Dict[str, Any]]:
    stringio = logging.getLogger("rosetta").handlers[0].stream
    cleaned = []
    previous = None
    for entry in stringio.getvalue().split('\n'):
        rex = re.match(f'\[(.*)\] (.*) - (.*)', entry)
        if rex:
            if previous:
                cleaned.append(previous)
            previous = dict(datetime=rex.group(1),
                            level=logging.getLevelName(rex.group(2)),
                            text=rex.group(3))
        else:
            previous['text'] += '\n' + entry
    if previous:
        cleaned.append(previous)
    return cleaned
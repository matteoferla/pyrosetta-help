import logging, sys, io, pyrosetta
import re
from typing import Union

def get_logger():
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

def get_log_entries(levelname: Union[str, int]=logging.INFO):
    if isinstance(levelname, int):
        # logging.INFO is actually an int, not an enum
        levelname = logging.getLevelName(levelname)
    stringio = logging.getLogger("rosetta").handlers[0].stream
    return re.findall(f'(\[.*\] {levelname} - [\w\W]*)', stringio.getvalue())
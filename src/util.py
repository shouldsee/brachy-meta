import path
import os

import pymisca.header
from filelock import FileLock
import json
from collections import OrderedDict as _dict
import collections
import time
SRC_DIR = os.path.dirname(__file__)

STATIC_DIR = 'src/static'
INDEX_FILE = os.path.realpath("OUTPUT/index.json.list")

# import pymisca.ext as pyext
from pymisca.events import LinkEvent,CopyEvent
def _file_static_copy( IN_DIR, IN_BASENAME, OUT_DIR, OUT_BASENAME):
    if IN_BASENAME is None:
        IN_BASENAME = IN_DIR
    if OUT_BASENAME is None:
        OUT_BASENAME = '.'.join(IN_BASENAME.lstrip('/').split('/'))  #### caution_1
    return CopyEvent(
        os.path.join(IN_DIR,IN_BASENAME),
        os.path.join(OUT_DIR,OUT_BASENAME),
        force=1)

ON_SERVER = 1
# ON_SERVER = 0


def _get_output_file(FNAME,INDEX_FILE=INDEX_FILE):
    d = [
        ("TS",time.time()),
        ("OUTPUT_FILE",FNAME,),
        ("RUNTIME_FILE", pymisca.header.name__lookup('__file__',level=-1)),
    ]
    d = _dict(d)
    print (json.dumps(d,indent=4))
    with FileLock(  INDEX_FILE +'.lock'):
        with open(INDEX_FILE, "a") as f:
#             f.write(json.dumps(d)+'\n')
            json.dump(d, f)
            f.write('\n')
    return FNAME
get_output_file = _get_output_file

def _get_file(FNAME):
    if ON_SERVER:
#         FNAME = '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/%s'%
        FNAME = FNAME
    else:
        FNAME = path.Path(SRC_DIR)/ "static" / FNAME.lstrip('/').replace('/','.') ##### caution_1
    assert os.path.exists(FNAME),FNAME
    ###
    return FNAME

def _get_middle_file(FNAME):
#     FNAME = 'OUTPUT/'+FNAME
    assert os.path.exists(FNAME),FNAME
    return FNAME
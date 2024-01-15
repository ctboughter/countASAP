# Yet another wild way to do this, but way easier than
# trying to configure things within each app/cli to find
# the test data... Plus it lets users play with it a bit more
import countASAP
import os
import shutil

def getit():
    # This basically finds out *where* the countASAP package is installed
    # then you can get the path, and copy the examples out of it.
    aims_dir = countASAP.__file__[:-11]
    startDir = os.getcwd()
    shutil.copytree(aims_dir+'ex_inputs',startDir+'/ex_inputs')
__author__ = 'joohwi'

import os
from subprocess import *

def run_process(cmdline, logfile = "", append=False):
    cmd = cmdline.split(" ")
    p = Popen(cmd, stdout=PIPE, shell=False)
    output = p.communicate()
    if logfile != "":
        f = open(logfile, "a" if append else "w")
        f.writelines(output)
        f.close()
    return output[0]
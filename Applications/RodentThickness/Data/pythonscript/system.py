__author__ = 'joohwi'

import os
from subprocess import *

def run_process(cmdline):
    cmd = cmdline.split(" ")
    p = Popen(cmd, stdout=PIPE, shell=False)
    output = p.communicate()
    return output[0]
import system
import sys
import logging

system.setup_logger(logfilename="test.log")
retcode = system.run_process("/devel/linux/SPHARM-PDM/spharm-pdm-build-joowhi/bin/GenParaMeshCLP %s %s %s" % (sys.argv[1], sys.argv[2], sys.argv[3]), verbose=True)
print "Return Code = ", retcode

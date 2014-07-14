__author__ = 'joohwi'

import os
import os.path
from subprocess import *
import shlex
import time
import sys
import logging
import fcntl


log_commands = False

# setup the global variable log_commands
# if log_commands is True, it will not execute the command but only logs it
def set_log_commands(setting):
  global log_commands
  log_commands = setting

def setup_logger(logfilename=None, logformat=None, loglevel=logging.INFO):
  if logformat is None:
    logformat="%(levelname)s %(asctime)s %(message)s"
  logging.basicConfig(filename=logfilename, format=logformat, level=loglevel)


def set_nonblocking(instream):
    fd = instream.fileno()
    fl = fcntl.fcntl(fd, fcntl.F_GETFL)
    fcntl.fcntl(fd, fcntl.F_SETFL, fl | os.O_NONBLOCK)
    return fl


def safe_readline(instream):
  try:
    return instream.readline().strip()
  except IOError:
    return ""


def find_program(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"').strip("'")
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def run_process(cmdline, loggerName="", verbose=False):
    # when logging is already set, skip the basic configuration
    logger = logging.getLogger(loggerName)
    loggingCMDS = logging.INFO+1
    logging.addLevelName(loggingCMDS, "CMDS")

    cmd = shlex.split(cmdline)
    exename = cmd[0]
    if exename[0] != os.sep:
      exename = find_program(exename)  
      if exename is None:
        logger.warning("can't find %s" % (cmd[0]))
        return -1
      cmd[0] = exename
      cmdline = " ".join(cmd)

    # if log_commands is True
    # just log the command line and return
    global log_commands
    if log_commands:
      logger.log(loggingCMDS, cmdline)
      return 0

    p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=False)

    # log the command line executing
    logger.log(loggingCMDS, cmdline)

    # make stdin a non-blocking file
    fl_stderr = set_nonblocking(p.stderr)
    fl_stdout = set_nonblocking(p.stdout)

    try:
      while p.poll() is None:
        # non-blocking stdout readline()
        stdline = safe_readline(p.stdout)
        while stdline:
          logger.info(stdline)
          if verbose:
            print stdline
          stdline = safe_readline(p.stdout)

        # non-blocking stderr readline()
        errline = safe_readline(p.stderr)
        while errline:
          logger.warning(errline)
          if verbose:
            print errline
          errline = safe_readline(p.stderr)
        time.sleep(0.5)

      stdline = " "
      while stdline:
        stdline = safe_readline(p.stdout)
        if stdline:
          logger.info(stdline)
          if verbose:
            print stdline

      errline = " "
      while errline:
        errline = safe_readline(p.stderr)
        if errline:
          logger.info(errline)
          if verbose:
            print errline

    except KeyboardInterrupt as ke:
      logger.info("Keyboard Interrupt!")
      raise ke

    logger.info("returncode[%s] = %d" % (cmd[0], 0 if p.returncode is None else p.returncode))
    return p.returncode


def rename_logfile(logfile):
  if not os.path.exists(logfile):
    return
  (filename, fileext) = os.path.splitext(logfile)
  timestamp = time.strftime("%Y%m%d_%H%M%S")
  os.rename(logfile, filename + "." +  timestamp + fileext)
  return

def run_forever():
  j = 0
  try:
    while True:
      j += 1
      sys.stdout.write(str(j))
      sys.stdout.write("\n")
      sys.stdout.flush()
      time.sleep(1)
  except KeyboardInterrupt as ke:
     pass


if __name__ == "__main__":
  run_forever()


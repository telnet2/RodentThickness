__author__ = 'joohwi'

import os
import os.path
from subprocess import *
import shlex
import time
import sys
import logging
import fcntl

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

check_commands = False

# setup the global variable check_commands
# if check_commands is True, it will not execute the command but only logs it
def set_check_commands(setting):
    global check_commands
    check_commands = setting


def setup_logger(logfilename=None, logformat=None, loglevel=logging.INFO):
    if logformat is None:
        logformat = "%(levelname)s %(asctime)s %(message)s"
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


def find_programs_in_path(program_list):
    found_programs = {p: find_program(p) for p in program_list}
    not_found_programs = [p for p in found_programs.keys() if found_programs[p] is None]
    return found_programs, not_found_programs

def run_process(cmdline, loggerName="", verbose=False):
    # when logging is already set, skip the basic configuration
    logger = logging.getLogger(loggerName)
    loggingCMDS = logging.INFO + 1
    logging.addLevelName(loggingCMDS, "CMDS")

    # if verbose print the cmdline
    if verbose:
        print bcolors.OKBLUE + cmdline + bcolors.ENDC

    cmd = shlex.split(cmdline)
    exename = cmd[0]
    if exename[0] != os.sep:
        exename = find_program(exename)
        if exename is None:
            raise RuntimeError("can't find %s" % (cmd[0]))
        cmd[0] = exename
        cmdline = " ".join(cmd)

    # if check_commands is True
    # just log the command line and return
    global check_commands
    if check_commands:
        logger.log(loggingCMDS, cmdline)
        return 0


    # create the process
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
                logger.warning(errline)
                if verbose:
                    print bcolors.WARNING + errline + bcolors.ENDC

    except KeyboardInterrupt as ke:
        logger.info("Keyboard Interrupt!")
        raise ke

    logger.info("returncode[%s] = %d" % (cmd[0], 0 if p.returncode is None else p.returncode))
    return p.returncode


def run_process_args(args):
    run_process(args[0], args[1], args[2])


def scatter_array(arr,n):
    return [arr[a:b] for a,b in zip(range(0,len(arr),n),range(n,len(arr)+n,n))]


def init_worker():
    import signal
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def run_process_pool(cmdlines, pool_size=2, loggerName="", verbose=False):
    from multiprocessing import Pool
    try:
        cmds = scatter_array(cmdlines, pool_size)
        print cmds
        for cmd in cmds:
            print cmd
            p = Pool(pool_size,init_worker)
            arg1 = [loggerName for n in range(len(cmd))] 
            arg2 = [verbose for n in range(len(cmd))] 
            args_list = zip(cmd,arg1,arg2)
            for args in args_list:
                p.apply_async(run_process_args, args)
            p.close()
            p.join()
    except KeyboardInterrupt as ke:
        print "KeyboardInterrupt, terminating workers ...",
        pool.terminate()
        pool.join()
        print " done" 
        raise ke
        

def rename_logfile(logfile):
    if not os.path.exists(logfile):
        return
    (filename, fileext) = os.path.splitext(logfile)
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    os.rename(logfile, filename + "." + timestamp + fileext)
    return


def configure_paths(dirs, programs_list, verbose=True, quit_on_error=True):
    os.environ["PATH"] = sys.path[0] + os.pathsep \
        + ((dirs + os.pathsep) if dirs else "") + os.environ["PATH"]

    found_programs, not_found_programs = find_programs_in_path(programs_list)
    config = {k+"Path": found_programs[k] for k in found_programs}

    if (verbose):
      print "\n".join(["%-40s: %s" % (k[:-4], config[k]) for k in config])
      if not_found_programs:
          print "\n%-10s: %s" % ("Error", "Some programs are not found!")
          print "%-10s: %s" % ("Todo", "Check your PATH environment variable " + \
                                "whether those programs are found!")
          if (quit_on_error):
            sys.exit(0)

    config["?"] = not_found_programs
    return config


    """
    config = {}
    reg = re.compile("[ ]*[Ss][Ee][Tt][ ]*\([ ]*([A-Za-z0-9_]*)[ ]+(.*)\)")
    with open(bmsfile) as f:
        bmslines = f.readlines()
    for line in bmslines:
        bmsmatch = reg.match(line)
        name = bmsmatch.group(1)
        path = bmsmatch.group(2).strip()
        config.update({name:path})
    
    for k in config:
        print os.path.basename(config[k]) 
    """


def is_file_exist(thefile, override=False):
    if override:
        return False
    if not os.path.exists(thefile):
        return False
    if os.stat(thefile).st_size == 0:
        return False
    return True


def is_file_newer(thefile, depfile, override=False):
    # thefile must not exist
    if override:
        return True
    if not os.path.exists(thefile):
        return True
    if os.stat(thefile).st_size == 0:
        return True
    if not os.path.exists(depfile):
        return False
    if os.stat(depfile).st_size == 0:
        return False
    #print thefile, os.path.getmtime(thefile), depfile, os.path.getmtime(depfile)
    if os.path.getmtime(depfile) > os.path.getmtime(thefile):
        return True
    return False




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
    #programs = ["ls", "ccmake", "ssh", "hello"]
    #print find_programs_in_path(programs)
    setup_logger()
    run_process_pool(["ls -alR", "ls -alR"], 2)


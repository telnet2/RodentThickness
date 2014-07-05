#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "unistd.h"
#include "errno.h"

const char* ENV_PYTHONHOME = "/tools/Python/Python-2.7.7";

const char* ENV_PATH = "/tools/bin_linux64:/tools/Rscript/R-2.15.2/bin:/tools/Rscript/R-2.15.2/lib64/bin:/NIRAL/tools/Qt/Qt4.7.4-linux64_THL/bin:/NIRAL/tools/Qt/Qt4.7.4-linux64_THL/lib:/tools/Python/Python-2.7.7/bin:/NIRAL/tools/CMake/cmake-2.8.10.2-Linux-i386/bin:/tools/BatchMake/BatchMake_1.3_linux64_3.20/bin:/usr/bin:/home/joohwi/bin:/usr/local/etc:/usr/local/bin:/usr/bin:/bin:/etc:/usr/bin/X11:/usr/X11/bin:/usr/etc:/usr/sbin:/sbin:/tools/Slicer3/Slicer3-3.6.3-2011-03-04-linux-x86_64:/tools/Slicer3/Slicer3-3.6.3-2011-03-04-linux-x86_64/bin:/tools/Slicer3/Slicer3-3.6.3-2011-03-04-linux-x86_64/lib/Slicer3/Plugins:/tools/Slicer3/Slicer3-3.6.3-2011-03-04-linux-x86_64/lib/MRML";

const char* ENV_LD_LIBRARY_PATH = "/NIRAL/tools/VTK/VTK_5.10.1_dyn_rel/lib/vtk-5.10:/NIRAL/tools/Qt/Qt4.7.4-linux64_THL/lib:/tools/Python/Python-2.7.7/lib:/tools/BatchMake/BatchMake_1.3_linux64_3.20/lib:/usr/lib:/usr/local/lib/X11:/usr/local/lib:/opt/local/lib64:/opt/local/lib:/usr/ucblib:/home/joohwi/lib/linux64:.:/tools/lib_linux64:/tools/lib_linux:/tools/Slicer3/Slicer3-3.6.3-2011-03-04-linux-x86_64/lib/FreeSurfer:/tools/Slicer3/Slicer3-3.6.3-2011-03-04-linux-x86_64/lib/KWWidgets:/tools/Slicer3/Slicer3-3.6.3-2011-03-04-linux-x86_64/lib/MRML:/tools/Slicer3/Slicer3-3.6.3-2011-03-04-linux-x86_64/lib/MRMLCLI:/tools/Slicer3/Slicer3-3.6.3-2011-03-04-linux-x86_64/lib/Slicer3:/tools/Slicer3/Slicer3-3.6.3-2011-03-04-linux-x86_64/lib/Teem-1.11.0:/tools/Slicer3/Slicer3-3.6.3-2011-03-04-linux-x86_64/lib/vtkTeem:/tools/Slicer3/Slicer3-3.6.3-2011-03-04-linux-x86_64/lib/Slicer3/Plugins:/tools/Slicer3/Slicer3-3.6.3-2011-03-04-linux-x86_64/lib/LoadableModule:/tools/Slicer3/Slicer3-3.6.3-2011-03-04-linux-x86_64/lib/ITKCommandIO:/tools/Slicer3/Slicer3-3.6.3-2011-03-04-linux-x86_64/lib/MRMLIDImageIO:/tools/Slicer3/Slicer3-3.6.3-2011-03-04-linux-x86_64/lib/MGHImageIO:/tools/Slicer3/Slicer3-3.6.3-2011-03-04-linux-x86_64/lib/IGT:/tools/Slicer3/Slicer3-3.6.3-2011-03-04-linux-x86_64/lib/RemoteIO:/tools/Slicer3/Slicer3-3.6.3-2011-03-04-linux-x86_64/lib/InsightToolkit:/tools/Slicer3/Slicer3-3.6.3-2011-03-04-linux-x86_64/lib/BatchMake:/tools/Slicer3/Slicer3-3.6.3-2011-03-04-linux-x86_64/lib/bmModuleDescriptionParser:/tools/Slicer3/Slicer3-3.6.3-2011-03-04-linux-x86_64/lib/ModuleDescriptionParser";

char buf[1024];
char exebuf[512];
char exebuf2[512];

void readpath() {
  FILE* file = fopen("PATH.env", "r");
  fgets(buf, sizeof(buf), file);
  int buflen = strlen(buf);
  buf[buflen-1] = 0;
  fclose(file);
}

void extractpath(char* path) {
  int l = strlen(path);
  for (int i = l - 1; i >= 0; i--) {
    if (path[i] == '/') {
      path[i] = 0;
      return;
    }
  }
}

int main(int argc, char* argv[]) {
  int err = 0;
  setenv("PATH", ENV_PATH, 1);
  setenv("LD_LIBRARY_PATH", ENV_LD_LIBRARY_PATH, 1);
  setenv("PYTHONHOME", ENV_PYTHONHOME, 1);
  setenv("PYTHONPATH", NULL, 1);
//  readlink("/proc/self/exe", exebuf, sizeof(exebuf));
 // extractpath(exebuf);
  sprintf(exebuf2, "%s/bin/python", ENV_PYTHONHOME);
  err = execvp(exebuf2, argv);
  return err;
}

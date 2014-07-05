#include <iostream>
#include <itksys/SystemTools.hxx> 

using namespace std;
using namespace itksys;

int main(int argc, char* argv[]) {
  mode_t ITKmode_F_OK = 0;
  cout << "Permission:ITKmode_F_OK => " << itksys::SystemTools::GetPermissions(argv[1], ITKmode_F_OK) << endl;
  return 0;
}

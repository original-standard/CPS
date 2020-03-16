// Initialize the Nuc3pt argument structure
//#include <alg/NprDer_arg.h>
#include "NprDer_arg.h"
#include <util/gjp.h>
#include <util/smalloc.h>
CPS_START_NAMESPACE
NprDerArg::NprDerArg(): 
  num_masses(0),
  t_src(0)
{
  cname = "NprDerArg";
  max_num_iter=10000;
  stop_rsd=1e-8;
  p2max = 3;
  mass[0] = 0.02;
  num_source = 1;
  source = new int*[num_source];
  for(int i = 0;i < num_source;i++)
    source[i] = new int[4] {6*i,6*i,6*i,16*i};
  mom = 6;
}
int NprDerArg::NumMasses(){ return num_masses ; }

Float NprDerArg::Mass(int m){ return mass[m] ; }

void NprDerArg::check_args()
{
  char *fname = "check_args" ;
  char *msg = "src location error\n";
}
CPS_END_NAMESPACE

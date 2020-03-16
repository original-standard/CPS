#include <config.h>
//-------------------------------------------------------------------
//
//  Includes delta S = 1 and 2 contractions, with non-degnerate mass
//  capabilities.  Propagators are stored by shifting around the
//  machine to avoid recalculating propagators.
//
//  Random number stream is presumed loaded into memory.  These
//  values are then restored when the generator is initialized
//  by copying them into the appropriate location in the code load.
//  The initial load to memory is unchanged.
//
//  For K->pipi decay I=3/2
//-------------------------------------------------------------------


#include <stdio.h>
#include <stdlib.h>

#include <comms/sysfunc_cps.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/qcdio.h>
#include <alg/alg_fix_gauge.h>
#include <alg/alg_ghb.h>
#include <alg/alg_plaq.h>
#include <alg/do_arg.h>
#include <alg/common_arg.h>
#include <alg/ghb_arg.h>
#include <alg/no_arg.h>
#include <alg/fix_gauge_arg.h>
#include <util/ReadLatticePar.h>
#include <util/WriteLatticePar.h>
#include "alg_NprDer.h"
#include <util/command_line.h>

//GlobalJobParameter GJP;
//Verbose VRB;
//Error ERR;

USING_NAMESPACE_CPS
//-------------------------------------------------
//  Scan values from the command line at position n
//-------------------------------------------------
 
//extern int GetInt( int n, char *argv[], FILE * fp );
//extern int GetHex( int n, char *argv[], FILE * fp );
//extern float GetFloat( int n, char *argv[], FILE * fp );

//#include "get_arg.CC"

//----------------------------------------------------------------------
// main code for production NprDer job
//----------------------------------------------------------------------

main(int argc,char *argv[])
{
 
  char *cname = argv[0] ;
  char *fname = "main()" ;

  char *filename;
  filename = (char*) smalloc( 128*sizeof(char) );

  CommandLine::is(argc,argv);

  //--------------------------------------
  // Initialize all Global Job Parameters
  //--------------------------------------

  DoArg do_arg;

  //-------------------------------
  // Initialize argument structures
  //-------------------------------

  CommonArg common_arg;
  CommonArg common_arg_ghb;
  CommonArg common_arg_plaq;
  NprDerArg NprDer_arg;
  FixGaugeArg fix_arg; 
  NoArg plaq_arg;



  //------------------------------------------------------------------
  //  Check argument count
  //------------------------------------------------------------------
 
  if (argc < 6)
    ERR.General(cname, fname, "Should be 7 args, not %d.\n", argc+1) ;
 
  // working directory
  const char* dir(CommandLine::arg());
  strcpy(filename,dir);
  strcat(filename,"/do_arg.default");
 
  do_arg.Decode(filename,"do_arg");
  strcpy(filename, dir);
  strcat(filename,"/do_arg.dat");
  do_arg.Encode(filename,"do_arg");

  int conf_flag=0;
  if( do_arg.start_conf_kind == START_CONF_FILE ){
    conf_flag=1;
    do_arg.start_conf_kind = START_CONF_ORD;
  }

  //-------------------------
  // Initialize the GJP class
  //-------------------------
  
  Start(&argc,&argv);
  GJP.Initialize(do_arg);

  const char* conf_file(CommandLine::arg());
  if( conf_flag ) do_arg.start_conf_kind = START_CONF_FILE;

  //-------------------------------------------------------------------
  //  Output file holding WME results.  Also gets a copy of arguments
  //-------------------------------------------------------------------

  int loop_s = CommandLine::arg_as_int();
  int loop_e = CommandLine::arg_as_int();
  int loop_skp = CommandLine::arg_as_int();

  //-----------------------------------------------------------------
  // Gauge fixing args
  //-----------------------------------------------------------------

  fix_arg.fix_gauge_kind = FIX_GAUGE_LANDAU;
  fix_arg.hyperplane_start = 0;
  fix_arg.hyperplane_step = 1;
  fix_arg.hyperplane_num = do_arg.t_nodes * do_arg.t_node_sites;
  fix_arg.max_iter_num = 8000;
  fix_arg.stop_cond = 1e-8;

  //------------------
  // Set verbose level
  //------------------

  VRB.Level(6);



  //--------------------------------------------------------------------
  //  Important:  create lattice, then restore random stream
  //  Follow with plaquette calculation, gauge fixing and NprDer code.
  //--------------------------------------------------------------------




  //--------------------------------------------------------------------
  //
  // conf. loop start
  //
  //--------------------------------------------------------------------
  for(int loop=loop_s; loop<=loop_e; loop+=loop_skp){
    GwilsonFdwf lat;
  //-----------------------------------------------------------------
  // NprDer args
  //-----------------------------------------------------------------
    strcpy(filename,dir);
    strcat(filename,"/npr_der.default");
    
    NprDer_arg.Decode(filename,"NprDer_arg");
    strcpy(filename, dir);
    strcat(filename,"/npr_der.dat");
    sprintf(filename,"%s.%d",filename,loop);
    NprDer_arg.Encode(filename,"NprDer_arg");
    
    
    sprintf(filename,"%s/NPR_Der.dat.%d",dir,loop);
    common_arg.results = filename;
    
    FILE *fp; //open output file
    if( (fp = Fopen((char *)common_arg.results, "a")) == NULL ) {
      ERR.FileA(cname, fname, (char *)common_arg.results );
    }
    
#if TARGET == QCDOC
    Fprintf( fp, "Machine size in nodes:  (X,Y,Z,T) = (%d,%d,%d,%d)\n",
	     SizeX(), SizeY(), SizeZ(), SizeT() );
#else
    Fprintf( fp, "Scalar Machine\n" ) ;
#endif
    Fprintf( fp, "\nCommand line arguments\n");
    Fprintf( fp, "\tprogram name:\t%s\n", argv[0] );
    

    char lat_fname[200];
    if ( do_arg.start_conf_kind == START_CONF_FILE ) {
      sprintf(lat_fname,"%s.%d",conf_file,loop) ;
      ReadLatticeParallel read;
      printf("Before Readlattice from %s\n",lat_fname);
      read.read( lat, lat_fname );
      printf("End Readlattice from %s\n",lat_fname);
    }
    
    Fclose(fp);//output file
    

    AlgPlaq plaq(lat,&common_arg,&plaq_arg);
    plaq.run();
    

    AlgFixGauge fix_gauge(lat,&common_arg,&fix_arg);
    printf("Before fix_gauge\n");
    fix_gauge.run();
    printf("After fix_gauge\n");
    //  fix_gauge.free();
    
    AlgNprDer NprDer(lat,&common_arg,&NprDer_arg);
    NprDer.run();
    
  /*
  if ( do_arg.start_conf_kind == START_CONF_FILE ) {
    sfree((Matrix *)do_arg.start_conf_load_addr);
  }
  */

    fix_gauge.free();


  } // loop end
  sfree(filename);
  //sync();
}

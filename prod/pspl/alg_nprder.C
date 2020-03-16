#include <config.h>
//------------------------------------------------------------------
//
// alg_kpipi.C
//
// AlgNprDer is derived from Alg and is relevant to  
// three point correlation functions with Wilson-type 
// fermions (i.e. Wilson, Clover, Dwf). 
// The type of fermion is determined by the argument to the 
// constructor but it should not be a non Wilson type lattice.
//
//
//------------------------------------------------------------------
// Modification History.
// 09/04/03
// Chateau got this code from Noaki.
// Chateau Modified that for two-pion spectroscopy and K->pipi decay
// two-pion spectroscopy is necessary for \pi\pi state normalizaiton
//
//------------------------------------------------------------------

#include <stdlib.h>	// exit()
#include <stdio.h>

#include <comms/sysfunc_cps.h>

//#include <alg/alg_NprDer.h>
#include "alg_NprDer.h"
#include <comms/scu.h>
#include <util/qcdio.h>
#include <alg/alg_rnd_gauge.h>
#include "meson.h"

CPS_START_NAMESPACE

int siteOffset(const int lcl[], const int lcl_sites[]);
void NprDerStatus( int i );

//------------------------------------------------------------------
// Constructor 
//------------------------------------------------------------------
AlgNprDer::AlgNprDer(Lattice& latt, 
		     CommonArg *c_arg,
		     NprDerArg *arg) : 
		     Alg(latt, c_arg) 
{
  cname = "AlgNprDer";
  char *fname = "AlgNprDer(L&,CommonArg*,NprDerArg*)";
  VRB.Func(cname,fname);

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_NprDer_arg = arg;
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgNprDer::~AlgNprDer() {
  char *fname = "~AlgNprDer()";
  VRB.Func(cname,fname);

  //???
}

static  int NprDer_status_count;
//-------------------------------------------------------------------
// run_cg_contract_xyz() : do the cg and contraction in memory  
//------------------------------------------------------------------
void AlgNprDer::run() 
{
  char *fname = "run_cg_contract()";
  VRB.Func(cname,fname);

  // Set the Lattice pointer
  //------------------------
  Lattice& lat = AlgLattice();

  // set prop args
  //--------------
  QPropWArg prop_arg;
  //  prop_arg.cg = alg_NprDer_arg->cg;
  prop_arg.cg.stop_rsd = alg_NprDer_arg->stop_rsd;
  prop_arg.cg.max_num_iter = alg_NprDer_arg->max_num_iter;


  // gauge fix the props
  prop_arg.gauge_fix_src=0;
  prop_arg.gauge_fix_snk=0;
  //prop_arg.GaugeFixSrc=0;
  //prop_arg.GaugeFixSnk=0;

  // set p2 max
  //--------------

  char test_name[80];

  // Initialize boundary conditions 
  GJP.Xbc(BND_CND_PRD);
  GJP.Ybc(BND_CND_PRD);
  GJP.Zbc(BND_CND_PRD);
  GJP.Tbc(BND_CND_PRD);

  // compute mid_point propagator
  prop_arg.store_midprop = 1;
  // Save props explicitly after averaging
  //prop_arg.SaveProp = 0;
  prop_arg.save_ls_prop = 2;
  //prop_arg.do_half_fermion = 0;

  int time_size = GJP.Tnodes()*GJP.TnodeSites();
  int spatial_size = GJP.Xnodes()*GJP.XnodeSites();
  int ls_glb = GJP.Snodes()*GJP.SnodeSites();

  int p[4];

  //----------------------------------------------------------------------------------
  //
  // start including momentum loop 
  //
  //----------------------------------------------------------------------------------
  int mu,ms,sK,psgn; 
  int cmu, cms, cmuu;

  QPropW* prop;
  QPropW* cprop;

//-----------------------------------------------------------------------
// gauge transformation to Landau gauge
//-----------------------------------------------------------------------
  AlgRotateGauge rot( lat, common_arg );
  printf("Performed random gauge transform\n") ; 
  rot.run();


//-----------------------------------------------------------------------
// start u quark loop
//-----------------------------------------------------------------------
  for( mu = 0 ; mu < 1 ; mu++ ){
  int t_p = alg_NprDer_arg->t_src;

  //---------------------------------------------------------------------------------------
  // for light mass at t_p with p = 0 and rand 0
  //------------------------------------------------------------------------------------------------
    prop_arg.cg.mass = alg_NprDer_arg->mass[mu];
    GJP.Tbc(BND_CND_PRD);

    if (prop_arg.file == NULL) 
      prop_arg.file = (char*)smalloc(100*sizeof(char));
    sprintf(prop_arg.file,"not_read");


    for(int i = 0;i < alg_NprDer_arg->num_source;i++)
{
      {
	int x[4];
	int ** &source = alg_NprDer_arg->source;
	prop_arg.x = source[i][0];
	prop_arg.y = source[i][1];
	prop_arg.z = source[i][2];
	prop_arg.t = source[i][3];
	for(int j = 0; j < 4;j++)
	  x[j] = source[i][j];

	prop_arg.gauge_fix_snk=0;
	prop = new QPropWPointSrc( lat, &prop_arg, common_arg );
	Meson m("PI");
        m.setGamma(-5);
	m.setMass(prop_arg.cg.mass,prop_arg.cg.mass);
	m.Zero();
	m.calcMeson(*prop,*prop);
	if(common_arg->results != 0){
	  FILE *fp;
	  if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
	    ERR.FileA(cname,fname, (char *)common_arg->results);
	  }
	  m.Print(fp);
	  Fclose(fp);
	}
	m.setSrc("Axial");
	m.Zero();
	m.calcAxial(*prop,*prop);
	if(common_arg->results != 0){
	  FILE *fp;
	  if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
	    ERR.FileA(cname,fname, (char *)common_arg->results);
	  }
	  m.Print(fp);
	  Fclose(fp);
	}
	sprintf(test_name,"Fourier ");
	Fourier_transform(*prop,x,test_name);
	

	prop->DeleteQPropLs();
	delete prop;

      }

      for(int dir = 0;dir < 4;dir++)
      {
	int x[4];
	int ** &source = alg_NprDer_arg->source;
	prop_arg.x = source[i][0];
	prop_arg.y = source[i][1];
	prop_arg.z = source[i][2];
	prop_arg.t = source[i][3];
	for(int j = 0; j < 4;j++)
	  x[j] = source[i][j];

	prop_arg.gauge_fix_snk=0;
	
	prop = new QPropWPSPLTSrc( lat, &prop_arg, common_arg ,dir);
	sprintf(test_name,"Fourier split_dir=%d ",dir);
	Fourier_transform(*prop,x,test_name);
	

	prop->DeleteQPropLs();
	delete prop;

      }
      //    cprop = new QPropW(*prop);
}

    /*
    if(GJP.Snodes()==2) {
      //cprop->CopyQPropLs(*prop);
      //cprop->SwapQPropLs();
    }


    if(GJP.Snodes()==2) cprop->DeleteQPropLs();
    delete cprop;
    */


  } // light quark loop End


  sfree(prop_arg.file);
}

void AlgNprDer::Fourier_transform(QPropW& prop,int * source, char* memo)
{
  char *fname = "Fourier_transform()";
  VRB.Clock(cname,fname);

  int time_size=GJP.Tnodes()*GJP.TnodeSites();
  int spatial_size=GJP.Xnodes()*GJP.XnodeSites();
  int time_size_h = time_size/2;
  int spatial_size_h = spatial_size/2;

  Float pi=4*atan(1.0);

  Float pt_u = 2 * pi / time_size;
  Float pt_u_2 = pt_u * pt_u;
  Float px_u = 2 * pi / spatial_size;
  Float px_u_2 = px_u * px_u;

  mom = alg_NprDer_arg->mom;
  int count = (mom * 2 + 1) * (mom * 2 + 1) * (mom * 2 + 1) * (mom * 2 + 1);

  WilsonMatrix* wmat=(WilsonMatrix*)smalloc(count*sizeof(WilsonMatrix));
  if(wmat == 0) ERR.Pointer(cname,fname,"wmat");
  VRB.Smalloc(cname,fname,"wmat",wmat,count*sizeof(WilsonMatrix));

 
  /*
  for(int t = 0; x < time_size;t++)
    for(int x = 0; x < spatial_size;x++)
      for(int y = 0; y < spatial_size;y++)
	for(int z = 0; z < spatial_size;z++){
	  Float temp = t * t * pt_u_2 + x * x * px_u_2 + y * y * px_u_2 + z * z * px_u_2;
	  count 
	}
  */
  int p[count][4];
  for(int ij=0;ij<count;ij++) 
    {
      int l = ij;
      wmat[ij]=0.0;
      for(int i = 0;i < 4;i++){
	p[ij][i] = l % (mom * 2 + 1) - mom;
	l /= (mom * 2 + 1);
      }
    }
  Site site;
  for(site.Begin();site.End();site.nextSite()){

    int x(site.physX()-source[0]), y(site.physY()-source[1]), z(site.physZ()-source[2]), t(site.physT()-source[3]);

    for(int ipn=0; ipn<count; ipn++){
      Float theta = px_u * ( x*p[ipn][0] + y*p[ipn][1] + z*p[ipn][2] ) + pt_u * t*p[ipn][3] ;
      Rcomplex cc(cos(theta), -sin(theta)); // exp(-ipx)
      wmat[ipn]+=prop[site.Index()]*cc;
    }
  }

  for(int ij=0;ij<count;ij++) {
    slice_sum((Float*)&wmat[ij], 288, 99);
  }

  if(common_arg->results != 0){
    FILE *fp;
    if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    Fprintf(fp,"dirac_sink color_sink dirac_source color_source\n");

    for(int ipn=0; ipn<count; ipn++){
      Fprintf(fp,"%s px %d py %d pz %d pt %d x %d y %d z %d t %d\n",memo,p[ipn][0],p[ipn][1],p[ipn][2],p[ipn][3],source[0],source[1],source[2],source[3]);
      for(int dirac_sink=0;dirac_sink<4;dirac_sink++)
        for(int color_sink=0;color_sink<3;color_sink++)
          for(int dirac_source=0;dirac_source<4;dirac_source++)
            for(int color_source=0;color_source<3;color_source++){
              Rcomplex cc(wmat[ipn](dirac_sink,color_sink,dirac_source,color_source));
              Fprintf(fp,"%d %d %d %d %-25.15e %-25.15e\n",dirac_sink,color_sink,dirac_source,color_source,cc.real(),cc.imag());
            }
    }
sync();
    Fclose(fp);
  }else {
    printf("no specific file\n");
  }
  sfree(wmat);
  
}

CPS_END_NAMESPACE

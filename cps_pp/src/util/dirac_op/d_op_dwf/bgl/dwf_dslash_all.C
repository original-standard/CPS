#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008/04/21 14:19:18 $
//  $Header: /space/cvs/cps/cps++/src/util/dirac_op/d_op_dwf/bgl/dwf_dslash_all.C,v 1.4 2008/04/21 14:19:18 chulwoo Exp $
//  $Id: dwf_dslash_all.C,v 1.4 2008/04/21 14:19:18 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: dwf_dslash_all.C,v $
//  $Revision: 1.4 $
//  $Source: /space/cvs/cps/cps++/src/util/dirac_op/d_op_dwf/bgl/dwf_dslash_all.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// dwf_dslash_4.C
//
// dwf_dslash_4 is the derivative part of the space-time part of
// the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice
// cb = 0/1 <--> even/odd checkerboard of in field.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<config.h>
#include<util/dwf.h>
#include<util/wilson.h>
#include<util/gjp.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
CPS_START_NAMESPACE
#ifdef PARALLEL
CPS_END_NAMESPACE
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE
#endif
void dwf_dslash_all(Vector *out, 
		  Matrix *gauge_field, 
		  Vector *in, 
		  Float mass,
		  int cb, 
		  int dag, 
		  Dwf *dwf_lib_arg)
{
  int i;
  int ls;
  IFloat *frm_in;
  IFloat *frm_out;
  IFloat *g_field;
  Wilson *wilson_p;
  int size_cb[2];
  int parity;

  //----------------------------------------------------------------
  // Initializations
  //----------------------------------------------------------------
  ls = dwf_lib_arg->ls;
  frm_in = (IFloat *) in;
  frm_out = (IFloat *) out;
  g_field = (IFloat *) gauge_field;
  wilson_p = dwf_lib_arg->wilson_p;
  size_cb[0] = 24*wilson_p->vol[0];
  size_cb[1] = 24*wilson_p->vol[1];
  
  //----------------------------------------------------------------
  // Apply 4-dimensional Dslash
  //----------------------------------------------------------------
#if TARGET == BGL
  int vec_len=1;
#else
  int vec_len=2;
#endif
  dwf_dslash_5_plus_start(out,in,mass,dag,dwf_lib_arg);
  for(i=0; i<ls; i+= vec_len){

    // parity of 4-D checkerboard
    //------------------------------------------------------------
    parity = (i + cb) % 2;

    // Apply on 4-dim "parity" checkerboard part
    //------------------------------------------------------------
    if(vec_len==1)
      wilson_dslash(frm_out, g_field, frm_in, parity, dag, wilson_p);
    else
      wilson_dslash_two(frm_out, frm_out+size_cb[parity], g_field, frm_in, frm_in+size_cb[parity], parity, 1-parity,dag, wilson_p);
#if 1
    for(int j=0; j<vec_len; j++ ) {
      dwf_dslash_5_plus_slice(out,in,mass,dag,dwf_lib_arg,i+j);
    }
#endif
    frm_in = frm_in + vec_len*size_cb[parity];
    frm_out = frm_out + vec_len*size_cb[parity];
  }

#if 0
  for(i=0; i<ls; i+=1 ) {
    dwf_dslash_5_plus_slice(out,in,mass,dag,dwf_lib_arg,i);
  }
#endif


}



CPS_END_NAMESPACE

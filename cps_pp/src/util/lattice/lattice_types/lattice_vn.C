#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of the GpowerRectFnone class.

  $Id: lattice_vn.C,v 1.4 2004/08/18 11:58:06 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004/08/18 11:58:06 $
//  $Header: /space/cvs/cps/cps++/src/util/lattice/lattice_types/lattice_vn.C,v 1.4 2004/08/18 11:58:06 zs Exp $
//  $Id: lattice_vn.C,v 1.4 2004/08/18 11:58:06 zs Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $Revision: 1.4 $
//  $Source: /space/cvs/cps/cps++/src/util/lattice/lattice_types/lattice_vn.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// lattice_vn.C
//
// This class has double inheritance. The virtual
// base class is Lattice with two types of derived classes. One type
// is relevant to the gauge part and has a name that begins with 
// "G". The other type is relevant to the fermion part and has
// a name that begins with "F". The class below inherits from one 
// GpowerRect and Fnone.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// PowerRect gauge action -- no fermion action
//------------------------------------------------------------------
GpowerRectFnone::GpowerRectFnone()
{
  cname = "GpowerRectFnone";
  char *fname = "GpowerRectFnone()";
  VRB.Func(cname,fname);

  //???
}

GpowerRectFnone::~GpowerRectFnone()
{
  char *fname = "~GpowerRectFnone()";
  VRB.Func(cname,fname);

  //???
}



CPS_END_NAMESPACE

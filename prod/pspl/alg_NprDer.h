//------------------------------------------------------------------
//
// alg_NprDer.h
//
// Header file for all alg classes relevant to Wilson-type fermion
// three point functions. The type of glue or fermion is given as
// an argument of type Lattice& to the constructor. If the 
// fermion type is not F_CLASS_WILSON or F_CLASS_CLOVER or 
// F_CLASS_DWF the constructors exit with a general error.
//
// Modified by chateau 
// based on Jun's threept
//
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_KPIPI_H
#define INCLUDED_ALG_KPIPI_H

#include <alg/common_arg.h>
#include <alg/alg_base.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
#include <comms/glb.h>
#include <alg/qpropw.h>
//#include "qpropw.h"
#include "NprDer_arg.h"

CPS_START_NAMESPACE

//------------------------------------------------------------------
//
// AlgNprDer is derived from Alg and is relevant to  
// meson three point functions with Wilson type 
// fermions (i.e. Wilson, Clover, Dwf). 
// The type of fermion is determined by the argument to the 
// constructor. If the fermion type is not F_CLASS_WILSON or 
// F_CLASS_CLOVER or F_CLASS_DWF the constructors exit with a 
// general error.
//
//------------------------------------------------------------------

class AlgNprDer : public Alg
{
 private:
    char* cname;

    NprDerArg* alg_NprDer_arg;
        // The argument structure for the
        // three point calculation
    int f_size;
        // Node checkerboard size of the fermion field
    int mom;
 public:
    AlgNprDer(Lattice & latt, CommonArg* c_arg, NprDerArg* arg);

    virtual ~AlgNprDer();

    void run();

    void Fourier_transform(QPropW&,int *,char *);

};

CPS_END_NAMESPACE

#endif




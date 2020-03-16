/*  NprDer_arg.h */

/*  The structure type NprDerArg holds the parameters specific to
    meson three point functions for Wilson type fermion's  */

#ifndef INCLUDED_NPRDER_ARG_H
#include <config.h>
#include <util/vml/types.h>
#include <util/vml/vml.h>
#include <util/enum.h>
#include <util/defines.h>
#define INCLUDED_NPRDER_ARG_H

//#include <alg/cg_arg.h>
//#include <util/vector.h>

CPS_START_NAMESPACE

class VML;
class NprDerArg {
public:
         bool Encode(char *filename,char *instance);
         bool Decode(char *filename,char *instance);
         bool Vml(VML *vmls,char *instance);
        char *cname;
	int num_masses;	// number of masses to do
	Float mass[10];	// the list of light masses
	int max_num_iter;     // max iteration number
	Float stop_rsd;       // stopping condition
	int t_src;	        // source times for quarks
	Float p2max;          // max. p2 for Fourier transfomation
	int mom; // max mom
	int num_source;
	int **source;
           NprDerArg (  ) ;
           void check_args (  ) ;
           int NumMasses (  ) ;
           void NumMasses (  int n ) ;
           Float Mass (  int m ) ;
};

/* the xdr functions */

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__STDC__) || defined(__cplusplus)
extern  bool_t vml_NprDerArg (VML *, char *instance, NprDerArg*);

#else /* K&R C */
extern  bool_t vml_NprDerArg (VML *, char *instance, NprDerArg*);

#endif /* K&R C */

#ifdef __cplusplus
}
#endif
CPS_END_NAMESPACE

#endif /* !INCLUDED_NPRDER_ARG_H */

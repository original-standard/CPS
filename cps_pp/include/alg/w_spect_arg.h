/*
 * Please do not edit this file.
 * It was generated by CPC make system.
 */

#ifndef _W_SPECT_ARG_H_RPCGEN
#define _W_SPECT_ARG_H_RPCGEN

#include <rpc/rpc.h>

#include <util/vml/vml.h>
#include <config.h>
#include <util/enum.h>
#include <util/defines.h>

#ifdef __cplusplus
extern "C" {
#endif

CPS_START_NAMESPACE

class WspectOutput {
public:
	 WspectOutput(char *filename);
	 void Encode(char *filename,char *instance);
	 void Decode(char *filename,char *instance);
	 void Vml(VML *vmls,char *instance);
	WbaryonFold fold;
	char *cg;
	char *cg2;
	char *pbp;
	char *mid_point;
	char *a0_p;
	char *a1;
	char *b1;
	char *pion;
	char *pion_prime;
	char *rho;
	char *rho_prime;
	char *a0;
	char *a0_prime;
	char *a1_x;
	char *a1_y;
	char *a1_z;
	char *b1_x;
	char *b1_y;
	char *b1_z;
	char *rho_x;
	char *rho_y;
	char *rho_z;
	char *rho_x_prime;
	char *rho_y_prime;
	char *rho_z_prime;
	char *nucleon;
	char *nucleon_prime;
	char *delta_x;
	char *delta_y;
	char *delta_z;
	char *delta_t;
	   WspectOutput (  ) ;
};

class WspectArg {
public:
	 WspectArg(char *filename);
	 void Encode(char *filename,char *instance);
	 void Decode(char *filename,char *instance);
	 void Vml(VML *vmls,char *instance);
	int prop_dir;
	int num_mom;
	SourceKind source_kind;
	int src_box_b[4];
	int src_box_e[4];
	Float g_epsi;
	int g_n;
	int g_center[4];
	Float rescale_factor;
	int aots_num;
	int aots_start;
	int aots_step;
	int baryons_on;
	int normal_mesons_on;
	int extended_mesons_on;
	int extended_mesonsBE_on;
	int extended_mesons_op_groupId;
	int extended_mesons_first_dev_on;
	int extended_mesons_second_sym_dev_on;
	int extended_mesons_second_antisym_dev_on;
	int extended_mesons_second_diag_dev_on;
	int fuzzing_on;
	int sink_fuzzing_only;
	int fuzzing_level;
	int fuzzing_c_num;
	Float fuzzing_c[MAX_FUZZING_C_NUM];
	int fuzzing_hits;
	int extended_mesonsBE_op_groupId;
	int extended_mesonsBE_Elec_on;
	int extended_mesonsBE_Magn_on;
	int BEfuzzing_on;
	int BEfuzzing_level;
	int BEfuzzing_c_num;
	Float BEfuzzing_c[MAX_FUZZING_C_NUM];
	int BEfuzzing_hits;
	   WspectArg (  ) ;
};

/* the xdr functions */

#if defined(__STDC__) || defined(__cplusplus)
extern  bool_t vml_WspectOutput (VML *, char *instance, WspectOutput*);
extern  bool_t vml_WspectArg (VML *, char *instance, WspectArg*);

#else /* K&R C */
extern  bool_t vml_WspectOutput (VML *, char *instance, WspectOutput*);
extern  bool_t vml_WspectArg (VML *, char *instance, WspectArg*);

#endif /* K&R C */

#ifdef __cplusplus
}
#endif
CPS_END_NAMESPACE

#endif /* !_W_SPECT_ARG_H_RPCGEN */

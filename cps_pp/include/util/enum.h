/*
 * Please do not edit this file.
 * It was generated by CPC make system.
 */

#ifndef _ENUM_H_RPCGEN
#define _ENUM_H_RPCGEN

#include <config.h>
#include <util/vml/types.h>
#include <util/vml/vml.h>
#include <util/enum.h>
#include <util/defines.h>
CPS_START_NAMESPACE

typedef double Float;

typedef double IFloat;

enum DirType {
	DIR_X = 0,
	DIR_Y = 1,
	DIR_Z = 2,
	DIR_T = 3,
	DIR_S = 4,
};
typedef enum DirType DirType;

enum FclassType {
	F_CLASS_NONE = 0,
	F_CLASS_STAG = 1,
	F_CLASS_WILSON = 2,
	F_CLASS_CLOVER = 3,
	F_CLASS_DWF = 4,
	F_CLASS_ASQTAD = 5,
	F_CLASS_P4 = 6,
};
typedef enum FclassType FclassType;

enum GclassType {
	G_CLASS_NONE = 0,
	G_CLASS_WILSON = 1,
	G_CLASS_POWER_PLAQ = 2,
	G_CLASS_IMPR_RECT = 3,
	G_CLASS_POWER_RECT = 4,
	G_CLASS_IMPR_OLSYM = 5,
};
typedef enum GclassType GclassType;

enum StrOrdType {
	CANONICAL = 0,
	STAG = 1,
	WILSON = 2,
	G_WILSON_HB = 3,
};
typedef enum StrOrdType StrOrdType;

enum CnvFrmType {
	CNV_FRM_NO = 0,
	CNV_FRM_YES = 1,
};
typedef enum CnvFrmType CnvFrmType;

enum FermionFieldDimension {
	FOUR_D = 0,
	FIVE_D = 1,
};
typedef enum FermionFieldDimension FermionFieldDimension;

enum PreserveType {
	PRESERVE_NO = 0,
	PRESERVE_YES = 1,
};
typedef enum PreserveType PreserveType;

enum StartConfType {
	START_CONF_ORD = 0,
	START_CONF_DISORD = 1,
	START_CONF_FILE = 2,
	START_CONF_LOAD = 3,
	START_CONF_MEM = 4,
};
typedef enum StartConfType StartConfType;

enum StartSeedType {
	START_SEED_FIXED = 0,
	START_SEED_FIXED_UNIFORM = 1,
	START_SEED = 2,
	START_SEED_UNIFORM = 3,
	START_SEED_INPUT = 4,
	START_SEED_INPUT_UNIFORM = 5,
	START_SEED_INPUT_NODE = 6,
};
typedef enum StartSeedType StartSeedType;

enum ChkbType {
	CHKB_EVEN = 0,
	CHKB_ODD = 1,
};
typedef enum ChkbType ChkbType;

enum DagType {
	DAG_NO = 0,
	DAG_YES = 1,
};
typedef enum DagType DagType;

enum BndCndType {
	BND_CND_PRD = 0,
	BND_CND_APRD = 1,
};
typedef enum BndCndType BndCndType;

enum FixGaugeType {
	FIX_GAUGE_NONE = -2,
	FIX_GAUGE_LANDAU = -1,
	FIX_GAUGE_COULOMB_X = 0,
	FIX_GAUGE_COULOMB_Y = 1,
	FIX_GAUGE_COULOMB_Z = 2,
	FIX_GAUGE_COULOMB_T = 3,
};
typedef enum FixGaugeType FixGaugeType;

enum SprojType {
	SPROJ_XM = 0,
	SPROJ_YM = 1,
	SPROJ_ZM = 2,
	SPROJ_TM = 3,
	SPROJ_XP = 4,
	SPROJ_YP = 5,
	SPROJ_ZP = 6,
	SPROJ_TP = 7,
};
typedef enum SprojType SprojType;

enum SigmaprojType {
	SIGMAPROJ_XY = 0,
	SIGMAPROJ_XZ = 1,
	SIGMAPROJ_XT = 2,
	SIGMAPROJ_YZ = 3,
	SIGMAPROJ_YT = 4,
	SIGMAPROJ_YX = 5,
	SIGMAPROJ_ZT = 6,
	SIGMAPROJ_ZX = 7,
	SIGMAPROJ_ZY = 8,
	SIGMAPROJ_TX = 9,
	SIGMAPROJ_TY = 10,
	SIGMAPROJ_TZ = 11,
};
typedef enum SigmaprojType SigmaprojType;

enum RitzMatType {
	NONE = 0,
	MAT_HERM = 1,
	MATPC_HERM = 2,
	MATPCDAG_MATPC = 3,
	NEG_MATPCDAG_MATPC = 4,
	MATDAG_MAT = 5,
	NEG_MATDAG_MAT = 6,
	MATDAG_MAT_NORM = 7,
	NEG_MATDAG_MAT_NORM = 8,
};
typedef enum RitzMatType RitzMatType;

enum RatApproxType {
	CONSTANT = 0,
	DYNAMIC = 1,
};
typedef enum RatApproxType RatApproxType;

enum MultiShiftSolveType {
	SINGLE = 0,
	MULTI = 1,
};
typedef enum MultiShiftSolveType MultiShiftSolveType;

enum WbaryonFold {
	BARYON_FOLD = 0,
	BARYON_RAW = 1,
	BARYON_PAST = 2,
};
typedef enum WbaryonFold WbaryonFold;

enum SourceKind {
	POINT_W = 0,
	WALL_W = 0 + 1,
	BOX_W = 0 + 2,
	JACOBI_W = 0 + 3,
	MAX_NUM_SINK = 0 + 4,
};
typedef enum SourceKind SourceKind;

enum MomentumKind {
	MOM_000 = 0,
	MOM_001 = 0 + 1,
	MOM_002 = 0 + 2,
	MOM_011 = 0 + 3,
	MOM_022 = 0 + 4,
	MOM_111 = 0 + 5,
	MOM_222 = 0 + 6,
	MAX_NUM_MOMENTA = 0 + 7,
};
typedef enum MomentumKind MomentumKind;

enum DEVOperatorKind {
	UNIT = 0,
	DEV1 = 0 + 1,
	DEV2 = 0 + 2,
	DEV3 = 0 + 3,
	DEV1DEV2 = 0 + 4,
	DEV2DEV1 = 0 + 5,
	DEV2DEV3 = 0 + 6,
	DEV3DEV2 = 0 + 7,
	DEV1DEV3 = 0 + 8,
	DEV3DEV1 = 0 + 9,
	DEV1DEV1 = 0 + 10,
	DEV2DEV2 = 0 + 11,
	DEV3DEV3 = 0 + 12,
	DEV_OP_NUM = 0 + 13,
	SUM_F = 0 + 14,
	SUM_S_ANTISYM = 0 + 15,
	SUM_S_SYM = 0 + 16,
	SUM_S_DIAG = 0 + 17,
	SUM_F_S_ANTISYM = 0 + 18,
	SUM_S_SYM_DIAG = 0 + 19,
	SUM_UNIT_F_S_ANTISYM = 0 + 20,
	END_SUM_OP = 0 + 21,
	BEGIN_BE_OP = 0 + 22,
	FB1_OP = 0,
	FB2_OP = 0 + 1,
	FB3_OP = 0 + 2,
	FE1_OP = 0 + 3,
	FE2_OP = 0 + 4,
	FE3_OP = 0 + 5,
	FUNIT_OP = 0 + 6,
	SUM_MAGN_OP = 0 + 7,
	SUM_ELEC_OP = 0 + 8,
	SUM_MAGN_ELEC_OP = 0 + 9,
	END_BE_OP = 0 + 10,
};
typedef enum DEVOperatorKind DEVOperatorKind;

enum WMesonOpKind {
	MO_a0xP_x = 0,
	MO_a0xP_y = 1,
	MO_a0xP_z = 2,
	MO_pionxP_x = 3,
	MO_pionxP_y = 4,
	MO_pionxP_z = 5,
	MO_a0_primexP_x = 6,
	MO_a0_primexP_y = 7,
	MO_a0_primexP_z = 8,
	MO_rhoxP_A1 = 9,
	MO_rhoxP_T1_x = 10,
	MO_rhoxP_T1_y = 11,
	MO_rhoxP_T1_z = 12,
	MO_rhoxP_T2_x = 13,
	MO_rhoxP_T2_y = 14,
	MO_rhoxP_T2_z = 15,
	MO_a1xP_A1 = 16,
	MO_a1xP_T2_x = 17,
	MO_a1xP_T2_y = 18,
	MO_a1xP_T2_z = 19,
	MO_a1xP_E_1 = 20,
	MO_a1xP_E_2 = 21,
	MO_b1xP_T1_x = 22,
	MO_b1xP_T1_y = 23,
	MO_b1xP_T1_z = 24,
	MO_b1xD_A2 = 25,
	MO_b1xD_T1_x = 26,
	MO_b1xD_T1_y = 27,
	MO_b1xD_T1_z = 28,
	MO_b1xD_T2_x = 29,
	MO_b1xD_T2_y = 30,
	MO_b1xD_T2_z = 31,
	MO_b1xD_E_1 = 32,
	MO_b1xD_E_2 = 33,
	MO_a0_primexD_x = 34,
	MO_a0_primexD_y = 35,
	MO_a0_primexD_z = 36,
	MO_rhoxB_T1_x = 37,
	MO_rhoxB_T1_y = 38,
	MO_rhoxB_T1_z = 39,
	MO_rhoxB_T2_x = 40,
	MO_rhoxB_T2_y = 41,
	MO_rhoxB_T2_z = 42,
	MO_a1xB_A1 = 43,
	MO_a1xB_T1_x = 44,
	MO_a1xB_T1_y = 45,
	MO_a1xB_T1_z = 46,
	MO_a1xB_T2_x = 47,
	MO_a1xB_T2_y = 48,
	MO_a1xB_T2_z = 49,
	MO_a1xD_A2 = 50,
	MO_a1xD_T1_x = 51,
	MO_a1xD_T1_y = 52,
	MO_a1xD_T1_z = 53,
	MO_a1xD_T2_x = 54,
	MO_a1xD_T2_y = 55,
	MO_a1xD_T2_z = 56,
	MO_a1xD_E_1 = 57,
	MO_a1xD_E_2 = 58,
	MO_rhoxD_A2 = 59,
	MO_rhoxD_T1_x = 60,
	MO_rhoxD_T1_y = 61,
	MO_rhoxD_T1_z = 62,
	MO_rhoxD_T2_x = 63,
	MO_rhoxD_T2_y = 64,
	MO_rhoxD_T2_z = 65,
	MO_pionxB_T1_x = 66,
	MO_pionxB_T1_y = 67,
	MO_pionxB_T1_z = 68,
	MO_pionxD_T2_x = 69,
	MO_pionxD_T2_y = 70,
	MO_pionxD_T2_z = 71,
	NUM_WMESON_OP_KIND = 72,
};
typedef enum WMesonOpKind WMesonOpKind;

enum WMesonState {
	MS_a0xP_x = 0,
	MS_a0xP_y = 1,
	MS_a0xP_z = 2,
	MS_pionxP_x = 3,
	MS_pionxP_y = 4,
	MS_pionxP_z = 5,
	MS_a0_primexP_x = 6,
	MS_a0_primexP_y = 7,
	MS_a0_primexP_z = 8,
	MS_rhoxP_A1_1 = 9,
	MS_rhoxP_T1_x = 10,
	MS_rhoxP_T1_y = 11,
	MS_rhoxP_T1_z = 12,
	MS_rhoxP_T2_x = 13,
	MS_rhoxP_T2_y = 14,
	MS_rhoxP_T2_z = 15,
	MS_a1xP_A1_1 = 16,
	MS_a1xP_T2_x = 17,
	MS_a1xP_T2_y = 18,
	MS_a1xP_T2_z = 19,
	MS_a1xP_E_1 = 20,
	MS_a1xP_E_2 = 21,
	MS_b1xP_T1_x = 22,
	MS_b1xP_T1_y = 23,
	MS_b1xP_T1_z = 24,
	MS_b1xD_A2_1 = 25,
	MS_b1xD_T1_x = 26,
	MS_b1xD_T1_y = 27,
	MS_b1xD_T1_z = 28,
	MS_b1xD_T2_x = 29,
	MS_b1xD_T2_y = 30,
	MS_b1xD_T2_z = 31,
	MS_b1xD_E_1 = 32,
	MS_b1xD_E_2 = 33,
	MS_a0_primexD_x = 34,
	MS_a0_primexD_y = 35,
	MS_a0_primexD_z = 36,
	MS_rhoxB_T1_x = 37,
	MS_rhoxB_T1_y = 38,
	MS_rhoxB_T1_z = 39,
	MS_rhoxB_T2_x = 40,
	MS_rhoxB_T2_y = 41,
	MS_rhoxB_T2_z = 42,
	MS_a1xB_A1_1 = 43,
	MS_a1xB_T1_x = 44,
	MS_a1xB_T1_y = 45,
	MS_a1xB_T1_z = 46,
	MS_a1xB_T2_x = 47,
	MS_a1xB_T2_y = 48,
	MS_a1xB_T2_z = 49,
	MS_a1xD_A2_1 = 50,
	MS_a1xD_T1_x = 51,
	MS_a1xD_T1_y = 52,
	MS_a1xD_T1_z = 53,
	MS_a1xD_T2_x = 54,
	MS_a1xD_T2_y = 55,
	MS_a1xD_T2_z = 56,
	MS_a1xD_E_1 = 57,
	MS_a1xD_E_2 = 58,
	MS_rhoxD_A2_1 = 59,
	MS_rhoxD_T1_x = 60,
	MS_rhoxD_T1_y = 61,
	MS_rhoxD_T1_z = 62,
	MS_rhoxD_T2_x = 63,
	MS_rhoxD_T2_y = 64,
	MS_rhoxD_T2_z = 65,
	MS_pionxB_T1_x = 66,
	MS_pionxB_T1_y = 67,
	MS_pionxB_T1_z = 68,
	MS_pionxD_T2_x = 69,
	MS_pionxD_T2_y = 70,
	MS_pionxD_T2_z = 71,
	NUM_WMESON_STATE = 72,
};
typedef enum WMesonState WMesonState;

enum WMesonOutputName {
	a0xP = 0,
	pionxP = 1,
	a0_primexP = 2,
	rhoxP_A1 = 3,
	rhoxP_T1 = 4,
	rhoxP_T2 = 5,
	a1xP_A1 = 6,
	a1xP_T2 = 7,
	a1xP_E = 8,
	b1xP_T1 = 9,
	b1xD_A2 = 10,
	b1xD_T1 = 11,
	b1xD_T2 = 12,
	b1xD_E = 13,
	a0_primexD = 14,
	rhoxB_T1 = 15,
	rhoxB_T2 = 16,
	a1xB_A1 = 17,
	a1xB_T1 = 18,
	a1xB_T2 = 19,
	a1xD_A2 = 20,
	a1xD_T1 = 21,
	a1xD_T2 = 22,
	a1xD_E = 23,
	rhoxD_A2 = 24,
	rhoxD_T1 = 25,
	rhoxD_T2 = 26,
	pionxB_T1 = 27,
	pionxD_T2 = 28,
	NUM_WMESON_OUTPUT = 29,
};
typedef enum WMesonOutputName WMesonOutputName;

enum WMesonCategory {
	NORMALMESON = 0,
	EXT_FIRSTDEV_MESON = 1,
	EXT_SECONDDEV_SYM_MESON = 2,
	EXT_SECONDDEV_ANTISYM_MESON = 3,
	EXT_SECONDDEV_DIAG_MESON = 4,
	MIXING = 5,
};
typedef enum WMesonCategory WMesonCategory;

enum WExtMesonBEOutputName {
	BE_pionxB = 0,
	BE_rhoxB_T1 = 0 + 1,
	NUM_WEXTMESON_BE_OUTPUT = 0 + 2,
};
typedef enum WExtMesonBEOutputName WExtMesonBEOutputName;

enum WExtMesonBEState {
	BE_MS_pionxB_x = 0,
	BE_MS_pionxB_y = 0 + 1,
	BE_MS_pionxB_z = 0 + 2,
	BE_MS_rhoxB_T1_x = 0 + 3,
	BE_MS_rhoxB_T1_y = 0 + 4,
	BE_MS_rhoxB_T1_z = 0 + 5,
	NUM_WEXTMESON_BE_STATES = 0 + 6,
};
typedef enum WExtMesonBEState WExtMesonBEState;

enum WExtMesonBEOp {
	BE_MO_pionxB_x = 0,
	BE_MO_pionxB_y = 0 + 1,
	BE_MO_pionxB_z = 0 + 2,
	BE_MO_rhoxB_T1_x = 0 + 3,
	BE_MO_rhoxB_T1_y = 0 + 4,
	BE_MO_rhoxB_T1_z = 0 + 5,
	NUM_WEXTMESON_BE_OPS = 0 + 6,
};
typedef enum WExtMesonBEOp WExtMesonBEOp;

enum WExtMesonBECategory {
	ELEC_HYBRID_BE = 0,
	MAG_HYBRID_BE = 0 + 1,
	MIXING_BE = 0 + 2,
};
typedef enum WExtMesonBECategory WExtMesonBECategory;

enum FieldTensorId {
	FB1 = 0,
	FB2 = 0 + 1,
	FB3 = 0 + 2,
	FE1 = 0 + 3,
	FE2 = 0 + 4,
	FE3 = 0 + 5,
	NUM_FLDS = 0 + 6,
	FUNIT = 0 + 7,
	SUM_MAGN = 0 + 8,
	SUM_ELEC = 0 + 9,
	SUM_MAGN_ELEC = 0 + 10,
	NUM_FLD_OPS = 0 + 11,
};
typedef enum FieldTensorId FieldTensorId;

/* the xdr functions */

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__STDC__) || defined(__cplusplus)
extern  bool_t vml_Float (VML *, char *instance, Float*);
extern  bool_t vml_IFloat (VML *, char *instance, IFloat*);
extern  bool_t vml_DirType (VML *, char *instance, DirType*);
extern  bool_t vml_FclassType (VML *, char *instance, FclassType*);
extern  bool_t vml_GclassType (VML *, char *instance, GclassType*);
extern  bool_t vml_StrOrdType (VML *, char *instance, StrOrdType*);
extern  bool_t vml_CnvFrmType (VML *, char *instance, CnvFrmType*);
extern  bool_t vml_FermionFieldDimension (VML *, char *instance, FermionFieldDimension*);
extern  bool_t vml_PreserveType (VML *, char *instance, PreserveType*);
extern  bool_t vml_StartConfType (VML *, char *instance, StartConfType*);
extern  bool_t vml_StartSeedType (VML *, char *instance, StartSeedType*);
extern  bool_t vml_ChkbType (VML *, char *instance, ChkbType*);
extern  bool_t vml_DagType (VML *, char *instance, DagType*);
extern  bool_t vml_BndCndType (VML *, char *instance, BndCndType*);
extern  bool_t vml_FixGaugeType (VML *, char *instance, FixGaugeType*);
extern  bool_t vml_SprojType (VML *, char *instance, SprojType*);
extern  bool_t vml_SigmaprojType (VML *, char *instance, SigmaprojType*);
extern  bool_t vml_RitzMatType (VML *, char *instance, RitzMatType*);
extern  bool_t vml_RatApproxType (VML *, char *instance, RatApproxType*);
extern  bool_t vml_MultiShiftSolveType (VML *, char *instance, MultiShiftSolveType*);
extern  bool_t vml_WbaryonFold (VML *, char *instance, WbaryonFold*);
extern  bool_t vml_SourceKind (VML *, char *instance, SourceKind*);
extern  bool_t vml_MomentumKind (VML *, char *instance, MomentumKind*);
extern  bool_t vml_DEVOperatorKind (VML *, char *instance, DEVOperatorKind*);
extern  bool_t vml_WMesonOpKind (VML *, char *instance, WMesonOpKind*);
extern  bool_t vml_WMesonState (VML *, char *instance, WMesonState*);
extern  bool_t vml_WMesonOutputName (VML *, char *instance, WMesonOutputName*);
extern  bool_t vml_WMesonCategory (VML *, char *instance, WMesonCategory*);
extern  bool_t vml_WExtMesonBEOutputName (VML *, char *instance, WExtMesonBEOutputName*);
extern  bool_t vml_WExtMesonBEState (VML *, char *instance, WExtMesonBEState*);
extern  bool_t vml_WExtMesonBEOp (VML *, char *instance, WExtMesonBEOp*);
extern  bool_t vml_WExtMesonBECategory (VML *, char *instance, WExtMesonBECategory*);
extern  bool_t vml_FieldTensorId (VML *, char *instance, FieldTensorId*);

#else /* K&R C */
extern  bool_t vml_Float (VML *, char *instance, Float*);
extern  bool_t vml_IFloat (VML *, char *instance, IFloat*);
extern  bool_t vml_DirType (VML *, char *instance, DirType*);
extern  bool_t vml_FclassType (VML *, char *instance, FclassType*);
extern  bool_t vml_GclassType (VML *, char *instance, GclassType*);
extern  bool_t vml_StrOrdType (VML *, char *instance, StrOrdType*);
extern  bool_t vml_CnvFrmType (VML *, char *instance, CnvFrmType*);
extern  bool_t vml_FermionFieldDimension (VML *, char *instance, FermionFieldDimension*);
extern  bool_t vml_PreserveType (VML *, char *instance, PreserveType*);
extern  bool_t vml_StartConfType (VML *, char *instance, StartConfType*);
extern  bool_t vml_StartSeedType (VML *, char *instance, StartSeedType*);
extern  bool_t vml_ChkbType (VML *, char *instance, ChkbType*);
extern  bool_t vml_DagType (VML *, char *instance, DagType*);
extern  bool_t vml_BndCndType (VML *, char *instance, BndCndType*);
extern  bool_t vml_FixGaugeType (VML *, char *instance, FixGaugeType*);
extern  bool_t vml_SprojType (VML *, char *instance, SprojType*);
extern  bool_t vml_SigmaprojType (VML *, char *instance, SigmaprojType*);
extern  bool_t vml_RitzMatType (VML *, char *instance, RitzMatType*);
extern  bool_t vml_RatApproxType (VML *, char *instance, RatApproxType*);
extern  bool_t vml_MultiShiftSolveType (VML *, char *instance, MultiShiftSolveType*);
extern  bool_t vml_WbaryonFold (VML *, char *instance, WbaryonFold*);
extern  bool_t vml_SourceKind (VML *, char *instance, SourceKind*);
extern  bool_t vml_MomentumKind (VML *, char *instance, MomentumKind*);
extern  bool_t vml_DEVOperatorKind (VML *, char *instance, DEVOperatorKind*);
extern  bool_t vml_WMesonOpKind (VML *, char *instance, WMesonOpKind*);
extern  bool_t vml_WMesonState (VML *, char *instance, WMesonState*);
extern  bool_t vml_WMesonOutputName (VML *, char *instance, WMesonOutputName*);
extern  bool_t vml_WMesonCategory (VML *, char *instance, WMesonCategory*);
extern  bool_t vml_WExtMesonBEOutputName (VML *, char *instance, WExtMesonBEOutputName*);
extern  bool_t vml_WExtMesonBEState (VML *, char *instance, WExtMesonBEState*);
extern  bool_t vml_WExtMesonBEOp (VML *, char *instance, WExtMesonBEOp*);
extern  bool_t vml_WExtMesonBECategory (VML *, char *instance, WExtMesonBECategory*);
extern  bool_t vml_FieldTensorId (VML *, char *instance, FieldTensorId*);

#endif /* K&R C */

#ifdef __cplusplus
}
#endif
CPS_END_NAMESPACE

#endif /* !_ENUM_H_RPCGEN */

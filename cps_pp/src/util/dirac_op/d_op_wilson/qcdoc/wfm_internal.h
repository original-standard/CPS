#ifndef _WFM_PAB_INTERNAL_H_
#define _WFM_PAB_INTERNAL_H_

/*--------------------------------------------------------------------------*/
/* Wilson assembler routines                                                */
/*--------------------------------------------------------------------------*/

/*
 * FP KERNELS
 */

extern "C" {

void dec_hsu3    (void *psi,void *gauge,void *len,void *tab);
void dec_hsu3_dag(void *psi,void *gauge,void *len,void *tab);

void rec_su3    (void *psi,void *gauge,void *chiin,void *len);
void rec_su3_dag(void *psi,void *gauge,void *chiin,void *len);


void vaxpy_norm (Float *scalep,
		 Float *InOutScale, 
		 Float *Add, 
		 int len,
		 Float *res);

void vaxpy (Float *scalep,Float *InOutScale, Float *Add, int len);
  /*Used to scatter a recv buf onto the appropriate face*/
void face_scatter(Float *TwoSpinor,
		  Float *RcvBuf, 
		  unsigned long*FaceTable,
		  unsigned long);
}


#endif

#include<config.h>
#include<qalloc.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Definition of glb_sum_five routine.

*/
//--------------------------------------------------------------------
//  CVS keywords
//  $Source: /space/cvs/cps/cps++/src/comms/qcdoc/scu/glb_sum2.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//--------------------------------------------------------------
//  glb_sum_internal2
//
//  An internal routine for glb_sum(), glb_sum_five()
//--------------------------------------------------------------


CPS_END_NAMESPACE
#include<comms/glb.h>
#include<comms/scu.h>
#include<util/gjp.h>
#include<util/checksum.h>
#include<util/data_shift.h>
#include<comms/double64.h>
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE

union DoubleBytes {
    Float dblval;
    char byte[8];
};

static Double64 *transmit_buf = NULL;
static Double64 *receive_buf = NULL;
static Double64 *gsum_buf = NULL;

//----------------------------------------------------------------------
/*!
  This routine need only be used by domain-wall fermion code where
  the 5th dimension is parallelised.
  
  \param float_p The number to be summed.
  \post The number pointed to by \a float_p is summed over all nodes
  and that sum is written back to \a float_p, which is identical on all nodes.

  \ingroup comms
*/
//---------------------------------------------------------------------- 

static int initted=0;
void glb_sum_internal2(Float * float_p,int ndir)
{
  static int NP[5] ={0,0,0,0,0};
  static int coor[5] ={0,0,0,0,0};
  static SCUDirArgIR *Send[5];
  static SCUDirArgIR *Recv[5];
  
  if (!initted){
      NP[0] = GJP.Xnodes();
      int max=NP[0];
      NP[1] = GJP.Ynodes();
      NP[2] = GJP.Znodes();
      NP[3] = GJP.Tnodes();
      NP[4] = GJP.Snodes();
      for (int i = 1;i<5;i++)
	if (max <NP[i]) max = NP[i];
      transmit_buf = (Double64 *)qalloc(QFAST|QNONCACHE,sizeof(Double64)*2);
	  if(!transmit_buf) ERR.Pointer("","glb_sum_internal2(p,i)","transmit_buf");
      receive_buf = transmit_buf+1;
//      gsum_buf = (Double64 *)qalloc(QFAST,sizeof(Double64)*max);
	  gsum_buf = (Double64 *)fmalloc("","glb_sum_internal2(p,i)","gsum_buf",sizeof(Double64)*max);
      for(int i = 0;i<5;i++)
      if (NP[i]>1){
      Send[i] = new SCUDirArgIR(transmit_buf, gjp_scu_dir[2*i+1], SCU_SEND, sizeof(Double64));
      Recv[i] = new SCUDirArgIR(receive_buf, gjp_scu_dir[2*i], SCU_REC, sizeof(Double64));
      }
  }
  initted = 1;

  // Sum over the "virtual" 5-dimensional mesh
  //------------------------------------------------------------
//  gsum_buf[0] = (Double64)*float_p;

  Double64 tmp_sum = (Double64)*float_p;

  // Save checksum of local floating points
  //---------------------------------------------------
  unsigned long sum_sum = local_checksum(float_p,1);
// CSM.AccumulateCsum(CSUM_GLB_LOC,sum_sum);


  for(int i = 0; i < ndir; ++i) 
  if (NP[i] >1) {
      int coor = (GJP.NodeCoor(i)-GDS.Origin(i)+NP[i])%NP[i];
//      fprintf(stdout,"GDS.Origin(%d)=%d coor[%d]=%d\n",i,GDS.Origin(i),i,coor);
      *transmit_buf = gsum_buf[coor]= tmp_sum;

      for (int itmp = 1; itmp < NP[i]; itmp++) {
	coor = (coor+1)%NP[i];
	Send[i]->StartTrans(); Recv[i]->StartTrans();
	Send[i]->TransComplete(); Recv[i]->TransComplete();

        gsum_buf[coor] = *receive_buf;
        *transmit_buf = *receive_buf;
      }
      tmp_sum = gsum_buf[0];
      for (int itmp = 1; itmp < NP[i]; itmp++) {
	    tmp_sum += gsum_buf[itmp];
      }
  }
  *float_p = (Float)tmp_sum;

  // accumulate final global sum checksum
  //------------------------------------
  sum_sum = local_checksum(float_p,1);
//  CSM.AccumulateCsum(CSUM_GLB_SUM,sum_sum);
  
}

static int initted_u=0;

static unsigned long long *transmit_buf_u = NULL;
static unsigned long long *receive_buf_u = NULL;
static unsigned long long *gsum_buf_u = NULL;

void glb_sum_internal2(unsigned int *uint_p, int ndir, int sum_flag) {
  static int NP[5] = {0,0,0,0,0};
  static int coor[5] = {0,0,0,0,0};
  static SCUDirArgIR *Send[5];
  static SCUDirArgIR *Recv[5];
  
  if (!initted_u) {
	NP[0] = GJP.Xnodes();
	int max = NP[0];
	NP[1] = GJP.Ynodes();
	NP[2] = GJP.Znodes();
	NP[3] = GJP.Tnodes();
	NP[4] = GJP.Snodes();
	for (int i = 1;i<5;i++)
	  if (max<NP[i]) max = NP[i];
	transmit_buf_u = (unsigned long long*)qalloc(QFAST|QNONCACHE,sizeof(unsigned long long)*2);
	if(!transmit_buf_u) ERR.Pointer("","glb_sum_internal2(p,i,i)","transmit_buf_u");
	receive_buf_u = transmit_buf_u+1;
	gsum_buf_u = (unsigned long long*)fmalloc("","glb_sum_internal2(p,i,i)","gsum_buf_u",sizeof(unsigned long long)*max);
//	  printf("transmit_buf_u=%p gsum_buf_u=%p\n",transmit_buf_u,gsum_buf_u);
	for(int i=0; i<5; i++)
      if (NP[i]>1) {
		Send[i] = new SCUDirArgIR(transmit_buf_u, gjp_scu_dir[2*i+1], SCU_SEND, sizeof(unsigned long long));
		Recv[i] = new SCUDirArgIR(receive_buf_u, gjp_scu_dir[2*i], SCU_REC, sizeof(unsigned long long));
      }
  }
  initted_u = 1;
  
  // Sum over the "virtual" 5-dimensional mesh
  //------------------------------------------------------------
  //  gsum_buf_u[0] = (unsigned long long)*float_p;
  
  unsigned int tmp_sum = *uint_p;
  
  for(int i=0; i<ndir; ++i) 
	if (NP[i] > 1) {
      int coor = GJP.NodeCoor(i);
//	  printf("coor[%d]=%d\n",i,coor);
      *transmit_buf_u = gsum_buf_u[coor] = (unsigned long long)tmp_sum;
	  
      for (int itmp = 1; itmp < NP[i]; itmp++) {
		coor = (coor+1)%NP[i];
//	    printf("coor[%d]=%d\n",i,coor);
		Send[i]->StartTrans(); Recv[i]->StartTrans();
		Send[i]->TransComplete(); Recv[i]->TransComplete();
		
        gsum_buf_u[coor] = *receive_buf_u;
        *transmit_buf_u = *receive_buf_u;
      }
      tmp_sum = (unsigned int)gsum_buf_u[0];
      for (int itmp = 1; itmp < NP[i]; itmp++) {
	if (sum_flag) tmp_sum += (unsigned int)gsum_buf_u[itmp];
	else  tmp_sum ^= (unsigned int)gsum_buf_u[itmp];
      }
	}
  *uint_p = tmp_sum;
}


CPS_END_NAMESPACE

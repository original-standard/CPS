#include<config.h>
#include<util/qcdio.h>
#include<util/data_shift.h>
#include<qalloc.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief Definition of glb_sum routine.
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /space/cvs/cps/cps++/src/comms/qcdoc/scu/glb_sum_internal.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//--------------------------------------------------------------
// glb_sum
//
// Sum over all nodes 
// {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes()}
//--------------------------------------------------------------
CPS_END_NAMESPACE
#include<comms/glb.h>
#include<comms/scu.h>
#include<util/gjp.h>
#include<comms/double64.h>
#include <comms/sysfunc_cps.h>
#include <comms/glb_sum_internal.h>
CPS_START_NAMESPACE



static Double64 *transmit_buf = NULL;
static Double64 *receive_buf = NULL;
static Double64 *gsum_buf = NULL;
static SCUDirArgIR *Send[5];
static SCUDirArgIR *Recv[5];
static SCUDirArgMulti *multi[5];
//const int MAX_BUF = 72;
static int output = 0;
static int initted = 0;
static int NP[5];
static int length[5];


//----------------------------------------------------------------------
/*!
  \param float_p The number to be summed.
  \post The number pointed to by \a float_p is summed over all nodes
  and that sum is written back to \a float_p, which is identical on all nodes.

  \ingroup comms
*/
//---------------------------------------------------------------------- 
void glb_sum_internal (Float * float_p, int dir,int len)
{
//  printf("float_p(%p)=%e dir=%d len=%d\n",float_p,*float_p,dir,len);
  if (dir<0 || dir>4)
  ERR.General("","glb_sum_internal","dir(%d) out of range",dir);
  if(len>MAX_BUF)
  ERR.General("","glb_sum_internal","len(%d)>MAX_BUF(%d)",len,MAX_BUF);
  if (!initted) {
    NP[0]=GJP.Xnodes();
	int max = NP[0];
    NP[1]=GJP.Ynodes();
    NP[2]=GJP.Znodes();
    NP[3]=GJP.Tnodes();
    NP[4]=GJP.Snodes();
    for (int i = 1;i<5;i++)
      if (max <NP[i]) max = NP[i]; 
    printf("max=%d\n",max);
    transmit_buf = (Double64 *)qalloc(QFAST|QNONCACHE,sizeof(Double64)*MAX_BUF*2);
    receive_buf = transmit_buf +MAX_BUF;
    gsum_buf = (Double64 *)qalloc(0,sizeof(Double64)*MAX_BUF*max);
    printf("gsum_buf=%p\n",gsum_buf);
    for(int i = 0;i<5;i++)
    if (NP[i] >1){
	  length[i] = 1;
      Send[i]= new SCUDirArgIR(transmit_buf,gjp_scu_dir[2*i+1],SCU_SEND,length[i]*sizeof(Double64));
       Recv[i]= new SCUDirArgIR(receive_buf,gjp_scu_dir[2*i],SCU_REC,length[i]*sizeof(Double64));
    }
    initted = 1;
  }
  if (NP[dir]==1) return;

  static Double64 tmp_sum[MAX_BUF];
  if (output) printf("glb_sum cpp before = %e ", (double)*float_p);
  if (length[dir] !=len && NP[dir]>1 ){
     printf("length[%d]=%d\n",dir,len);
		 length[dir] = len;
		 delete Send[dir]; 
         Send[dir]= new SCUDirArgIR(transmit_buf,gjp_scu_dir[2*dir+1],SCU_SEND,length[dir]*sizeof(Double64));
		 delete Recv[dir]; 
         Recv[dir]= new SCUDirArgIR(receive_buf,gjp_scu_dir[2*dir],SCU_REC,length[dir]*sizeof(Double64));
  }

  int coor = (GJP.NodeCoor(dir)-GDS.Origin(dir)+NP[dir])%NP[dir];
  //printf("coor(%d)=%d\n",dir,coor);
  for(int i = 0;i<len;i++){
  transmit_buf[i] = gsum_buf[coor*len+i] = (Double64)(float_p[i]);
 // printf("float_p[%d]=%e\n",i,float_p[i]);
  }

 //     *transmit_buf = *gsum_buf;

  if (NP[dir]>1){
    for (int itmp = 1; itmp < NP[dir]; itmp++) {
//      printf("dir=%d itmp=%d\n",dir,itmp);
      coor = (coor+1)%NP[dir];
      Send[dir]->StartTrans();
      Recv[dir]->StartTrans();
      Send[dir]->TransComplete();
      Recv[dir]->TransComplete();

      for(int i = 0;i<len;i++){
        gsum_buf[coor*len+i] = receive_buf[i];
        transmit_buf[i] = receive_buf[i];
      }
    }
  }
  for(int i = 0;i<len;i++){
    tmp_sum[i] = gsum_buf[i];
    if (NP[dir]>1)
    for( int itmp=1;itmp<NP[dir];itmp++){
      tmp_sum[i] += gsum_buf[itmp*len+i];
    }
    float_p[i] = (Float)tmp_sum[i];
//  printf("float_p[%d]=%e\n",i,float_p[i]);
  }
  if (output)   printf("after = %e\n", (double)*float_p);
}
CPS_END_NAMESPACE

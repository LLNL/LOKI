#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*
 Loki private version of the A++P++/PARTI supplied file.  To achieve good
 scaling of communication in A++P++ one should not throttle the number of
 messages in transit.  This may have been true when PARTI was written in the
 early 90's but it is no longer the case.  This throttling was controlled by
 setting GRPSIZE in this file to a small value.  In this version of the file
 it is set to MAX_NODES which is in turn defined by PARTI to be MAX_PROCESSORS,
 the maximum number of processors that A++/P++ was built to support, which is
 131072.  This is not an ideal fix for several reasons.  First there are a
 number of arrays in this file that are statically sized to GRPSIZE, now
 MAX_PROCESSORS.  So the functionality provided by this file now uses a lot
 more memory.  Secondly, if we really want to eliminate the limit on the number
 of messages in flight, this file can become tremendously simpler and most
 likely can be implemented with a minimal amount of required storage.  Finally,
 today 131072 is not an especially large number of processors but it will be
 sufficient for now.  Further hacking of A++/P++/PARTI headers to extend this
 limit is almost certainly not a good idea.
*/

/*
  This is not needed since the relavant macro is defined directly in port.h
  include <Parti_config.h>
 */

#include<stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "hash.h"
#include "bsparti.h"

/* This is redundent since port.h is included at the base of bsparti.h */
/* include"port.h" */

#define GRPSIZE		        MAX_NODES
#define NTYPE1			100
#define NODE_PID		0


#ifdef  SP1
int     PARTI_source;           /* source node of incoming message */
int     PARTI_type;             /* type of incoming message */
size_t     PARTI_nbytes;           /* number of bytes received */
#endif

/* everybody but PVM uses nonblocking sends/receives */

char *RmsgBuf[GRPSIZE],*SmsgBuf[GRPSIZE],	/* receive _and_ send buffers */
        *local_msgBuf;
int  RmsgBuf_size[GRPSIZE],SmsgBuf_size[GRPSIZE], 
	local_msgBuf_size;


void dDataMove(array, sched, dest)
  double            *array, *dest;
  SCHED    *sched;
{
/*
 * Nonblocking/Unbuffered version for the SP-1
 *	-Jim 6/94
 */

  int srcPrc, destPrc, lDest, lSrc, lPck, lWait, i,
      arraydims, sending, receiving, grp_size, GSsrc, GSdest, proc, 
      nproc, startPosn,Wprocs,allmesg,j1;

  int numPieces[MAX_DIM], k[MAX_DIM], step[MAX_DIM], currdim, j;

  int sGidVALID[GRPSIZE], rGidVALID[GRPSIZE], sizeR[GRPSIZE],
  sizeS[GRPSIZE], recv_done[GRPSIZE];

  double         *arrayPtr, *destPtr, *RbufPtr[GRPSIZE], *SbufPtr[GRPSIZE];

  double	 *tempPtr[MAX_DIM];

  SchedData *recvData, *sendData;

#ifdef PARTI_ENABLE_MP_INTERFACE_MPI
  MPI_Request sGid[GRPSIZE], rGid[GRPSIZE];
#else
  long sGid[GRPSIZE], rGid[GRPSIZE];
#endif

  if (sched == NULL)
    return;

  proc = PARTI_myproc();
  nproc = PARTI_numprocs();

  if( GRPSIZE > nproc )
    grp_size = nproc;
  else
    grp_size = GRPSIZE;

  Wprocs = (nproc / grp_size) * grp_size;

  /* Break 0..(nproc-1) loop into groups of size grp_size */
  for( srcPrc=0; srcPrc < nproc; srcPrc+=grp_size){
      if( srcPrc >= Wprocs   )
	GSsrc = nproc - Wprocs;
      else
	GSsrc = grp_size;

      GSdest = grp_size;

      if( (proc<(srcPrc+GSsrc))&&(proc>=srcPrc) )
	sending=1;
      else
	sending=0;

    for( destPrc=0; destPrc < nproc; destPrc+=grp_size){

      if( destPrc >= Wprocs   )
	GSdest = nproc - Wprocs;


      if( (proc<(destPrc+GSdest))&&(proc>=destPrc) )
	receiving=1;
      else
	receiving=0;


      /* check receive allocation */
      if( receiving ){
        for( lPck=0; lPck < GSsrc; lPck++ ){
          i = srcPrc+lPck;

          if (sched->rMsgSz[i] > 0) { 
            sizeR[lPck] = sched->rMsgSz[i] * DBL_SZ;
      
            if (i != proc) {          /* als - 5/19/92 */
              if (sizeR[lPck] > RmsgBuf_size[lPck]) {
            /* need a bigger message buffer (als - 4/93) */
                if (RmsgBuf[lPck] != NULL) free(RmsgBuf[lPck]);
                RmsgBuf[lPck] = (char *) malloc(sizeR[lPck]);
                if (RmsgBuf[lPck] == NULL) {
                  fatal_error("can't allocate space for message buffer in dataMove");
                }
                RmsgBuf_size[lPck] = sizeR[lPck];
              }
                RbufPtr[lPck]  = (double *) RmsgBuf[lPck];
/*		if ( proc == 0 ) printf("%d %p %lf\n", lPck, RbufPtr[lPck], *RbufPtr[lPck]); */

            } 
            else {
              /* copy from the local message buffer */
              if (sizeR[lPck] > local_msgBuf_size) {
                /* need a bigger local message buffer (als - 4/93) */
                if (local_msgBuf != NULL) free(local_msgBuf);
                local_msgBuf = (char *) malloc(sizeR[lPck]);
                if (local_msgBuf == NULL) {
                  fatal_error(
                    "can't allocate space for local message buffer in dataMove");
                }
                local_msgBuf_size = sizeR[lPck];
              }
              RbufPtr[lPck] = (double *) local_msgBuf;
/*	      if ( proc == 0 ) printf("%d %p %lf\n", lPck, RbufPtr[lPck], *RbufPtr[lPck]); */
           } 

         }
       }
     } 


      /* check sendbuffer allocation */
      if(sending){
        for( lPck=0; lPck < GSdest; lPck++ ){
	  i=lPck+destPrc;
          if (sched->sMsgSz[i] > 0) {
            sizeS[lPck] = sched->sMsgSz[i] * DBL_SZ;
            if (i != proc) {
            /* prepare a send */
              if (sizeS[lPck] > SmsgBuf_size[lPck]) {
                /* need a bigger message buffer (als - 4/93) */
                if (SmsgBuf[lPck] != NULL) free(SmsgBuf[lPck]);
                SmsgBuf[lPck] = (char *) malloc(sizeS[lPck]);
                if (SmsgBuf[lPck] == NULL) { fatal_error(
		  "can't allocate space for message buffer in dataMove"); }
                SmsgBuf_size[lPck] = sizeS[lPck];
              }
              SbufPtr[lPck] = (double *) SmsgBuf[lPck];
            }
            else {                    /* als - 5/19/92 */
            /* prepare the on-processor copy */
              if (sizeS[lPck] > local_msgBuf_size) {
              /* need a bigger local message buffer (als - 4/93) */
                if (local_msgBuf != NULL) free(local_msgBuf);
                local_msgBuf = (char *) malloc(sizeS[lPck]);
                if (local_msgBuf == NULL) {
                  fatal_error("can't allocate space for local message buffer in dataMove");
                }
                local_msgBuf_size = sizeS[lPck];
              }
              SbufPtr[lPck] = (double *) local_msgBuf;
            }

      /* pack sends */
            sendData = sched->sData[i];
        /* New Stuff goes here */
  
        /*  Gagan   07/03/93:  Rewritten to accomodate for
                               the new definition of schedule */

	    /* als 9/96 - added arbitrary dimensions for arrays */
  
            startPosn = sendData->startPosn ;
            arrayPtr  = array + startPosn ;
  
            arraydims = sendData->numDims ;

            if (arraydims > MAX_DIM) {
              fatal_error(
      	  "Datamove from array with more than MAX_DIM dimensions not allowed");
            }

	    for (j = 0; j < arraydims; j++) {
	      step[j] = sendData->str[j];
	      numPieces[j] = sendData->numelem[j] ;
	    }
	    /* set lower bounds for all loops - 0 here */
	    for (currdim = 0; currdim < arraydims; currdim++) {
	      k[currdim] = 0;
	      tempPtr[currdim] = arrayPtr;
	    }

	    for (currdim = 0; ; currdim = 0) {
	      /* reset indices and bump and reset buffer pointers for loops
		 that are done */
	      while (k[currdim] >= numPieces[currdim]) {
		if (currdim == arraydims-1) goto done1;

		k[currdim] = 0;
		currdim++;
		k[currdim]++;
		tempPtr[currdim] += step[currdim];
		/* reset all inner loop buffer pointers */
		for (j = currdim-1; j >= 0; j--) {
		  tempPtr[j] = tempPtr[j+1];
		}
	      }
	      /* inner loop code */
	      *SbufPtr[lPck] = *tempPtr[0] ;
	      SbufPtr[lPck]++ ;
	      tempPtr[0] += step[0] ;

	      k[0]++;

	    }
done1:      ;

            if(i != proc ){
	      /*
	      sGid[lPck] = PARTI_isend(NTYPE1+proc, SmsgBuf[lPck], 
                 sizeS[lPck], i, NULL );
	      */
#ifdef PARTI_ENABLE_MP_INTERFACE_MPI
	      PARTI_isend(NTYPE1+proc,SmsgBuf[lPck],sizeS[lPck],i,&sGid[lPck]);
#else
	      sGid[lPck] = PARTI_isend(NTYPE1+proc,SmsgBuf[lPck],sizeS[lPck],i,NULL);
#endif
	    }
        }
      }
    }
  

      /* receive */
      if(receiving){
        for( lSrc=0; lSrc < GSsrc; lSrc++ ){
          i=lSrc+srcPrc;
          if (sched->rMsgSz[i] > 0) {
  	      if(i!=proc){
		/*
	        rGid[lSrc] = PARTI_irecv(NTYPE1+i, RmsgBuf[lSrc], 
	          sizeR[lSrc], i, NULL );
		*/
#ifdef PARTI_ENABLE_MP_INTERFACE_MPI
	        PARTI_irecv(NTYPE1+i, RmsgBuf[lSrc], sizeR[lSrc], i, &rGid[lSrc]);
#else
	        rGid[lSrc] = PARTI_irecv(NTYPE1+i, RmsgBuf[lSrc], sizeR[lSrc], i, 
					 NULL );
#endif
	      }
            }
	}
      }

      if(receiving){      
	allmesg = 0;
	for(j1=0;j1<grp_size;j1++)
		recv_done[j1]=0;
	while(!allmesg){
	  allmesg=1;
	  for( lWait=0; lWait < GSsrc; lWait++ ){
	    if( !recv_done[lWait] ){
	      i=lWait+srcPrc;
	      if (sched->rMsgSz[i] > 0) {
		if(i!=proc){
		  /*
		  if(PARTI_msgdone(rGid[lWait])){
		  */
#ifdef PARTI_ENABLE_MP_INTERFACE_MPI
		  if(PARTI_msgdone(&rGid[lWait])){
#else
		  if(PARTI_msgdone(rGid[lWait])){
#endif
		    recv_done[lWait] = 1;
		  }
		  else{
		    allmesg=0;
		  }
		}
		else
		  recv_done[lWait] = 1;
	        if( recv_done[lWait]){


      /* unpack receive */
            recvData = sched->rData[i];
        /*
         * unbuffer received data
        */
  
        /*  Gagan   07/03/93:  Rewritten to accomodate for
                               the new definition of schedule */
  
	    /* als 9/96 - added arbitrary dimensions for arrays */
  
            startPosn = recvData->startPosn ;
            destPtr   = dest + startPosn ;
  
            arraydims = recvData->numDims ;

            if (arraydims > MAX_DIM) {
              fatal_error(
      	  "Datamove into array with more than MAX_DIM dimensions not allowed");
            }

	    for (j = 0; j < arraydims; j++) {
	      step[j] = recvData->str[j];
	      numPieces[j] = recvData->numelem[j] ;
	    }

	    /* set lower bounds for all loops - 0 here */
	    for (currdim = 0; currdim < arraydims; currdim++) {
	      k[currdim] = 0;
	      tempPtr[currdim] = destPtr;
	    }

	    for (currdim = 0; ; currdim = 0) {
	      /* reset indices and bump and reset buffer pointers for loops
		 that are done */
	      while (k[currdim] >= numPieces[currdim]) {
		if (currdim == arraydims-1) goto done2;

		k[currdim] = 0;
		currdim++;
		k[currdim]++;
		tempPtr[currdim] += step[currdim];
		/* reset all inner loop buffer pointers */
		for (j = currdim-1; j >= 0; j--) {
		  tempPtr[j] = tempPtr[j+1];
		}
	      }
	      /* inner loop code */
/*	      if ( proc == 0 ) printf("%d %p %lf\n", lWait, RbufPtr[lWait], *RbufPtr[lWait]); */
	      *tempPtr[0] = *RbufPtr[lWait]  ;
	      RbufPtr[lWait]++ ;
	      tempPtr[0] += step[0] ;
	      
	      k[0]++;
	    }
done2:      ;

	        }
	      }
	    }
	  }
	}
      }

      if(sending){
        for( lWait=0; lWait < GSdest; lWait++ ){
          i=lWait+destPrc;
          if (sched->sMsgSz[i] > 0) {
            if(i!=proc){
	      /*
              PARTI_msgwait(sGid[lWait]);
	      */
#ifdef PARTI_ENABLE_MP_INTERFACE_MPI
              PARTI_msgwait(&sGid[lWait]);
#else
              PARTI_msgwait(sGid[lWait]);
#endif
            }
          }
	}
      }
    }
  }
}


void iDataMove(array, sched, dest)
  int            *array, *dest;
  SCHED    *sched;
   {
  /*
   * Nonblocking/Unbuffered version for the SP-1
   *	-Jim 6/94
   */

     int srcPrc, destPrc, lDest, lSrc, lPck, lWait, i,
         arraydims, sending, receiving, grp_size, GSsrc, GSdest, proc, 
         nproc, startPosn,Wprocs,allmesg,j1;

     int numPieces[MAX_DIM], k[MAX_DIM], step[MAX_DIM], currdim, j;

     int sGidVALID[GRPSIZE], rGidVALID[GRPSIZE], sizeR[GRPSIZE],
         sizeS[GRPSIZE], recv_done[GRPSIZE];

     int         *arrayPtr, *destPtr, *RbufPtr[GRPSIZE], *SbufPtr[GRPSIZE];

     int         *tempPtr[MAX_DIM];

     SchedData *recvData, *sendData;

#ifdef PARTI_ENABLE_MP_INTERFACE_MPI
     MPI_Request sGid[GRPSIZE], rGid[GRPSIZE];
#else
     long sGid[GRPSIZE], rGid[GRPSIZE];
#endif

     if (PARTI_messagePassingInterpretationReport > 0)
        {
          printf ("Inside of iDataMove() (BEFORE SYNC) on processor %d \n",PARTI_myproc());
          MPI_Barrier(MPI_COMM_WORLD);
          printf ("Inside of iDataMove() (AFTER SYNC) on processor %d \n",PARTI_myproc());
/*
          printf ("Inside of iDataMove() Exiting (AFTER SYNC) on processor %d \n",PARTI_myproc());
          abort();
 */
        }

     if (sched == NULL)
        {
          if (PARTI_messagePassingInterpretationReport > 0)
               printf ("Leaving iDataMove because the sched == NULL on processor %d \n",PARTI_myproc());
          return;
        }
     
     proc = PARTI_myproc();
     nproc = PARTI_numprocs();

     if( GRPSIZE > nproc )
          grp_size = nproc;
       else
          grp_size = GRPSIZE;

     Wprocs = (nproc / grp_size) * grp_size;

  /* Break 0..(nproc-1) loop into groups of size grp_size */
     for( srcPrc=0; srcPrc < nproc; srcPrc+=grp_size)
        {
          if( srcPrc >= Wprocs   )
	       GSsrc = nproc - Wprocs;
            else
               GSsrc = grp_size;

          GSdest = grp_size;

          if( (proc<(srcPrc+GSsrc))&&(proc>=srcPrc) )
	       sending=1;
            else
	       sending=0;

          for( destPrc=0; destPrc < nproc; destPrc+=grp_size)
             {

               if( destPrc >= Wprocs   )
	            GSdest = nproc - Wprocs;

               if( (proc<(destPrc+GSdest))&&(proc>=destPrc) )
	            receiving=1;
                 else
                    receiving=0;

            /* check receive allocation */
               if( receiving )
                  {
                    if (PARTI_messagePassingInterpretationReport > 0)
                         printf ("Recieving messages in iDataMove() (setup buffers) on processor %d \n",PARTI_myproc());

                    for( lPck=0; lPck < GSsrc; lPck++ )
                       {
                         i = srcPrc+lPck;

                         if (sched->rMsgSz[i] > 0)
                            {
                              sizeR[lPck] = sched->rMsgSz[i] * INT_SZ;

                              if (i != proc)
                                 {          /* als - 5/19/92 */
                                   if (sizeR[lPck] > RmsgBuf_size[lPck])
                                      {
                                     /* need a bigger message buffer (als - 4/93) */
                                        if (RmsgBuf[lPck] != NULL)
                                             free(RmsgBuf[lPck]);
                                        RmsgBuf[lPck] = (char *) malloc(sizeR[lPck]);
                                        if (RmsgBuf[lPck] == NULL)
                                           {
                                             fatal_error("can't allocate space for message buffer in dataMove");
                                           }
                                        RmsgBuf_size[lPck] = sizeR[lPck];
                                      }
                                   RbufPtr[lPck]  = (int *) RmsgBuf[lPck];
                                 } 
                                else
                                 {
                                /* copy from the local message buffer */
                                   if (sizeR[lPck] > local_msgBuf_size)
                                      {
                                     /* need a bigger local message buffer (als - 4/93) */
                                        if (local_msgBuf != NULL) free(local_msgBuf);
                                             local_msgBuf = (char *) malloc(sizeR[lPck]);
                                        if (local_msgBuf == NULL)
                                           {
                                             fatal_error("can't allocate space for local message buffer in dataMove");
                                           }
                                        local_msgBuf_size = sizeR[lPck];
                                      }
                                   RbufPtr[lPck] = (int *) local_msgBuf;
                                 }
                            }
                       }
                  }
#if 1
                 else
                  {
                    if (PARTI_messagePassingInterpretationReport > 0)
                         printf ("NO MESSAGES TO RECIEVE in iDataMove() (setup buffers) on processor %d \n",PARTI_myproc());
                  }
#endif

            /* check sendbuffer allocation */
               if(sending)
                  {
                    if (PARTI_messagePassingInterpretationReport > 0)
                         printf ("Sending messages in iDataMove() (call PARTI_isend) on processor %d \n",PARTI_myproc());

                    for( lPck=0; lPck < GSdest; lPck++ )
                       {
                         i=lPck+destPrc;
                         if (sched->sMsgSz[i] > 0)
                            {
                              sizeS[lPck] = sched->sMsgSz[i] * INT_SZ;
                              if (i != proc)
                                 {
                                /* prepare a send */
                                   if (sizeS[lPck] > SmsgBuf_size[lPck])
                                      {
                                     /* need a bigger message buffer (als - 4/93) */
                                        if (SmsgBuf[lPck] != NULL)
                                             free(SmsgBuf[lPck]);
                                        SmsgBuf[lPck] = (char *) malloc(sizeS[lPck]);
                                        if (SmsgBuf[lPck] == NULL)
                                           {
                                             fatal_error("can't allocate space for message buffer in dataMove");
                                           }
                                        SmsgBuf_size[lPck] = sizeS[lPck];
                                      }
                                   SbufPtr[lPck] = (int *) SmsgBuf[lPck];
                                 }
                                else
                                 {                    /* als - 5/19/92 */
                                /* prepare the on-processor copy */
                                   if (sizeS[lPck] > local_msgBuf_size)
                                      {
                                     /* need a bigger local message buffer (als - 4/93) */
                                        if (local_msgBuf != NULL)
                                             free(local_msgBuf);
                                        local_msgBuf = (char *) malloc(sizeS[lPck]);
                                        if (local_msgBuf == NULL)
                                           {
                                             fatal_error("can't allocate space for local message buffer in dataMove");
                                           }
                                        local_msgBuf_size = sizeS[lPck];
                                      }
                                    SbufPtr[lPck] = (int *) local_msgBuf;
                                 }

                           /* pack sends */
                              sendData = sched->sData[i];
                           /* New Stuff goes here */
  
                          /*  Gagan   07/03/93:  Rewritten to accomodate for
                                                 the new definition of schedule */
  
                           /* als 9/96 - added arbitrary dimensions for arrays */
                              startPosn = sendData->startPosn ;
                              arrayPtr  = array + startPosn ;
  
                              arraydims = sendData->numDims ;

                              if (arraydims > MAX_DIM)
                                 {
                                   fatal_error("Datamove from array with more than MAX_DIM dimensions not allowed");
                                 }

                              for (j = 0; j < arraydims; j++)
                                 {
                                   step[j] = sendData->str[j];
                                   numPieces[j] = sendData->numelem[j];
                                 }
                           /* set lower bounds for all loops - 0 here */
                              for (currdim = 0; currdim < arraydims; currdim++)
                                 {
                                   k[currdim] = 0;
                                   tempPtr[currdim] = arrayPtr;
                                 }

                              for (currdim = 0; ; currdim = 0)
                                 {
                                /* reset indices and bump and reset buffer pointers for loops that are done */
                                   while (k[currdim] >= numPieces[currdim])
                                      {
                                        if (currdim == arraydims-1)
                                             goto done1;

                                        k[currdim] = 0;
                                        currdim++;
                                        k[currdim]++;
                                        tempPtr[currdim] += step[currdim];
                                     /* reset all inner loop buffer pointers */
                                        for (j = currdim-1; j >= 0; j--)
                                           {
                                             tempPtr[j] = tempPtr[j+1];
                                           }
                                      }
                                /* inner loop code */
                                   *SbufPtr[lPck] = *tempPtr[0] ;
                                   SbufPtr[lPck]++ ;
                                   tempPtr[0] += step[0] ;
                                   k[0]++;

                                 }
                       done1:      ;

                              if(i != proc )
                                 {
                                /* sGid[lPck] = PARTI_isend(NTYPE1+proc, SmsgBuf[lPck],sizeS[lPck], i, NULL ); */
#ifdef PARTI_ENABLE_MP_INTERFACE_MPI
                                   PARTI_isend(NTYPE1+proc,SmsgBuf[lPck],sizeS[lPck],i,&sGid[lPck]);
#else
                                   sGid[lPck] = PARTI_isend(NTYPE1+proc,SmsgBuf[lPck],sizeS[lPck],i,NULL);
#endif
                                 }
                            }
                       }
                  }
#if 1
                 else
                  {
                    if (PARTI_messagePassingInterpretationReport > 0)
                         printf ("NO MESSAGES TO SEND in iDataMove() (call PARTI_isend) on processor %d \n",PARTI_myproc());
                  }
#endif

            /* receive */
               if(receiving)
                  {
                    if (PARTI_messagePassingInterpretationReport > 0)
                         printf ("Recieving messages in iDataMove() (call PARTI_irecv) on processor %d \n",PARTI_myproc());

                    for( lSrc=0; lSrc < GSsrc; lSrc++ )
                       {
                         i=lSrc+srcPrc;
                         if (sched->rMsgSz[i] > 0)
                            {
                              if(i!=proc)
                                 {
                                /* rGid[lSrc] = PARTI_irecv(NTYPE1+i, RmsgBuf[lSrc], sizeR[lSrc], i, NULL ); */
#ifdef PARTI_ENABLE_MP_INTERFACE_MPI
                                   PARTI_irecv(NTYPE1+i, RmsgBuf[lSrc], sizeR[lSrc], i, &rGid[lSrc]);
#else
                                   rGid[lSrc] = PARTI_irecv(NTYPE1+i, RmsgBuf[lSrc], sizeR[lSrc], i, NULL );
#endif
                                 }
                            }
                       }
                  }
#if 1
                 else
                  {
                    if (PARTI_messagePassingInterpretationReport > 0)
                         printf ("NO MESSAGES TO RECIEVE in iDataMove() (call PARTI_irecv) on processor %d \n",PARTI_myproc());
                  }
#endif

               if(receiving)
                  {
                    if (PARTI_messagePassingInterpretationReport > 0)
                         printf ("Recieving messages in iDataMove() (check PARTI_msgdone) on processor %d \n",PARTI_myproc());

                    allmesg = 0;
                    for(j1=0;j1<grp_size;j1++)
                         recv_done[j1]=0;
                    while(!allmesg)
                       {
                         allmesg=1;
                         for( lWait=0; lWait < GSsrc; lWait++ )
                            {
                              if( !recv_done[lWait] )
                                 {
                                   i=lWait+srcPrc;
                                   if (sched->rMsgSz[i] > 0)
                                      {
                                        if(i!=proc)
                                           {
                                          /* if(PARTI_msgdone(rGid[lWait])){ */
#if 1
#ifdef PARTI_ENABLE_MP_INTERFACE_MPI
                                             if(PARTI_msgdone(&rGid[lWait]))
                                                {
#else
                                             if(PARTI_msgdone(rGid[lWait]))
                                                {
#endif
                                                  recv_done[lWait] = 1;
                                                }
                                               else
                                                {
                                                  allmesg=0;
                                                }
#else
                                             printf ("WARNING: PARTI_msgdone COMMENTED OUT on processor %d \n",PARTI_myproc());
#endif
                                           }
                                          else
                                             recv_done[lWait] = 1;
                                        if( recv_done[lWait])
                                           {
                                          /* unpack receive */

                                             recvData = sched->rData[i];
                                          /* unbuffer received data */

                                         /*  Gagan   07/03/93:  Rewritten to accomodate for
                                                                the new definition of schedule */
  
                                          /* als 9/96 - added arbitrary dimensions for arrays */

                                             startPosn = recvData->startPosn;
                                             destPtr   = dest + startPosn ;

                                             arraydims = recvData->numDims ;

                                             if (arraydims > MAX_DIM)
                                                {
                                                  fatal_error("Datamove into array with more than MAX_DIM dimensions not allowed");
                                                }

                                             for (j = 0; j < arraydims; j++)
                                                {
                                                  step[j] = recvData->str[j];
                                                  numPieces[j] = recvData->numelem[j] ;
                                                }

                                          /* set lower bounds for all loops - 0 here */
                                             for (currdim = 0; currdim < arraydims; currdim++)
                                                {
                                                  k[currdim] = 0;
                                                  tempPtr[currdim] = destPtr;
                                                }

                                             for (currdim = 0; ; currdim = 0)
                                                {
                                               /* reset indices and bump and reset buffer pointers for loops
                                                  that are done */
                                                  while (k[currdim] >= numPieces[currdim])
                                                     {
                                                       if (currdim == arraydims-1) goto done2;

                                                       k[currdim] = 0;
                                                       currdim++;
                                                       k[currdim]++;
                                                       tempPtr[currdim] += step[currdim];
                                                    /* reset all inner loop buffer pointers */
                                                       for (j = currdim-1; j >= 0; j--)
                                                          {
                                                            tempPtr[j] = tempPtr[j+1];
                                                          }
                                                     }
                                               /* inner loop code */
                                                  *tempPtr[0] = *RbufPtr[lWait];
                                                  RbufPtr[lWait]++ ;
                                                  tempPtr[0] += step[0] ;
                                                  k[0]++;
                                                }
                                        done2:  ;

                                           }
                                      }
                                 }
                            }
                       }
                  }
#if 1
                 else
                  {
                    if (PARTI_messagePassingInterpretationReport > 0)
                         printf ("NO MESSAGES TO SEND in iDataMove() (call PARTI_msgdone) on processor %d \n",PARTI_myproc());
                  }
#endif

               if(sending)
                  {
                    if (PARTI_messagePassingInterpretationReport > 0)
                         printf ("WAITING for messages (waiting on sends) in iDataMove() (call PARTI_msgwait) on processor %d \n",PARTI_myproc());
                    for( lWait=0; lWait < GSdest; lWait++ )
                       {
                         i=lWait+destPrc;
                         if (sched->sMsgSz[i] > 0)
                            {
                              if(i!=proc)
                                 {
                                /* PARTI_msgwait(sGid[lWait]); */
#ifdef PARTI_ENABLE_MP_INTERFACE_MPI
                                   PARTI_msgwait(&sGid[lWait]);
#else
                                   PARTI_msgwait(sGid[lWait]);
#endif
                                 }
                            }
                       }
                  }
#if 1
                 else
                  {
                    if (PARTI_messagePassingInterpretationReport > 0)
                         printf ("NO MESSAGES TO SEND in iDataMove() (call PARTI_msgwait) on processor %d \n",PARTI_myproc());
                  }
#endif
             }
        }

#if 0
     printf ("Before SYNC at base of iDataMove() on processor %d \n",PARTI_myproc());
     MPI_Barrier(MPI_COMM_WORLD);
     printf ("Exiting at base of datamove.C iDataMove() on processor %d \n",PARTI_myproc());
     abort();
#endif
   }


      
void fDataMove(array, sched, dest)
  float            *array, *dest;
  SCHED    *sched;
{
/*
 * Nonblocking/Unbuffered version for the SP-1
 *	-Jim 6/94
 */

  int srcPrc, destPrc, lDest, lSrc, lPck, lWait, i,
      arraydims, sending, receiving, grp_size, GSsrc, GSdest, proc, 
      nproc, startPosn,Wprocs,allmesg,j1;

  int numPieces[MAX_DIM], k[MAX_DIM], step[MAX_DIM], currdim, j;

  int sGidVALID[GRPSIZE], rGidVALID[GRPSIZE], sizeR[GRPSIZE],
  sizeS[GRPSIZE], recv_done[GRPSIZE];

  float         *arrayPtr, *destPtr, *RbufPtr[GRPSIZE], *SbufPtr[GRPSIZE];

  float         *tempPtr[MAX_DIM];


  SchedData *recvData, *sendData;

#ifdef PARTI_ENABLE_MP_INTERFACE_MPI
  MPI_Request sGid[GRPSIZE], rGid[GRPSIZE];
#else
  long sGid[GRPSIZE], rGid[GRPSIZE];
#endif

  if (sched == NULL)
    return;

  proc = PARTI_myproc();
  nproc = PARTI_numprocs();

  if( GRPSIZE > nproc )
    grp_size = nproc;
  else
    grp_size = GRPSIZE;

  Wprocs = (nproc / grp_size) * grp_size;

  /* Break 0..(nproc-1) loop into groups of size grp_size */
  for( srcPrc=0; srcPrc < nproc; srcPrc+=grp_size){
      if( srcPrc >= Wprocs   )
	GSsrc = nproc - Wprocs;
      else
	GSsrc = grp_size;

      GSdest = grp_size;

      if( (proc<(srcPrc+GSsrc))&&(proc>=srcPrc) )
	sending=1;
      else
	sending=0;

    for( destPrc=0; destPrc < nproc; destPrc+=grp_size){

      if( destPrc >= Wprocs   )
	GSdest = nproc - Wprocs;


      if( (proc<(destPrc+GSdest))&&(proc>=destPrc) )
	receiving=1;
      else
	receiving=0;


      /* check receive allocation */
      if( receiving ){
        for( lPck=0; lPck < GSsrc; lPck++ ){
          i = srcPrc+lPck;

          if (sched->rMsgSz[i] > 0) { 
            sizeR[lPck] = sched->rMsgSz[i] * FLT_SZ;
      
            if (i != proc) {          /* als - 5/19/92 */
              if (sizeR[lPck] > RmsgBuf_size[lPck]) {
            /* need a bigger message buffer (als - 4/93) */
                if (RmsgBuf[lPck] != NULL) free(RmsgBuf[lPck]);
                RmsgBuf[lPck] = (char *) malloc(sizeR[lPck]);
                if (RmsgBuf[lPck] == NULL) {
                  fatal_error("can't allocate space for message buffer in dataMove");
                }
                RmsgBuf_size[lPck] = sizeR[lPck];
              }
                RbufPtr[lPck]  = (float *) RmsgBuf[lPck];

            } 
            else {
              /* copy from the local message buffer */
              if (sizeR[lPck] > local_msgBuf_size) {
                /* need a bigger local message buffer (als - 4/93) */
                if (local_msgBuf != NULL) free(local_msgBuf);
                local_msgBuf = (char *) malloc(sizeR[lPck]);
                if (local_msgBuf == NULL) {
                  fatal_error(
                    "can't allocate space for local message buffer in dataMove");
                }
                local_msgBuf_size = sizeR[lPck];
              }
              RbufPtr[lPck] = (float *) local_msgBuf;
           } 

         }
       }
     } 


      /* check sendbuffer allocation */
      if(sending){
        for( lPck=0; lPck < GSdest; lPck++ ){
	  i=lPck+destPrc;
          if (sched->sMsgSz[i] > 0) {
            sizeS[lPck] = sched->sMsgSz[i] * FLT_SZ;
            if (i != proc) {
            /* prepare a send */
              if (sizeS[lPck] > SmsgBuf_size[lPck]) {
                /* need a bigger message buffer (als - 4/93) */
                if (SmsgBuf[lPck] != NULL) free(SmsgBuf[lPck]);
                SmsgBuf[lPck] = (char *) malloc(sizeS[lPck]);
                if (SmsgBuf[lPck] == NULL) { fatal_error(
		  "can't allocate space for message buffer in dataMove"); }
                SmsgBuf_size[lPck] = sizeS[lPck];
              }
              SbufPtr[lPck] = (float *) SmsgBuf[lPck];
            }
            else {                    /* als - 5/19/92 */
            /* prepare the on-processor copy */
              if (sizeS[lPck] > local_msgBuf_size) {
              /* need a bigger local message buffer (als - 4/93) */
                if (local_msgBuf != NULL) free(local_msgBuf);
                local_msgBuf = (char *) malloc(sizeS[lPck]);
                if (local_msgBuf == NULL) {
                  fatal_error("can't allocate space for local message buffer in dataMove");
                }
                local_msgBuf_size = sizeS[lPck];
              }
              SbufPtr[lPck] = (float *) local_msgBuf;
            }

      /* pack sends */
            sendData = sched->sData[i];
        /* New Stuff goes here */
  
        /*  Gagan   07/03/93:  Rewritten to accomodate for
                               the new definition of schedule */

	    /* als 9/96 - added arbitrary dimensions for arrays */

            startPosn = sendData->startPosn ;
            arrayPtr  = array + startPosn ;
  
            arraydims = sendData->numDims ;

            if (arraydims > MAX_DIM) {
              fatal_error(
      	  "Datamove from array with more than MAX_DIM dimensions not allowed");
            }

	    for (j = 0; j < arraydims; j++) {
	      step[j] = sendData->str[j];
	      numPieces[j] = sendData->numelem[j] ;
	    }
	    /* set lower bounds for all loops - 0 here */
	    for (currdim = 0; currdim < arraydims; currdim++) {
	      k[currdim] = 0;
	      tempPtr[currdim] = arrayPtr;
	    }

	    for (currdim = 0; ; currdim = 0) {
	      /* reset indices and bump and reset buffer pointers for loops
		 that are done */
	      while (k[currdim] >= numPieces[currdim]) {
		if (currdim == arraydims-1) goto done1;

		k[currdim] = 0;
		currdim++;
		k[currdim]++;
		tempPtr[currdim] += step[currdim];
		/* reset all inner loop buffer pointers */
		for (j = currdim-1; j >= 0; j--) {
		  tempPtr[j] = tempPtr[j+1];
		}
	      }
	      /* inner loop code */
	      *SbufPtr[lPck] = *tempPtr[0] ;
	      SbufPtr[lPck]++ ;
	      tempPtr[0] += step[0] ;

	      k[0]++;

	    }
done1:      ;

            if(i != proc ){
	      /*
	      sGid[lPck] = PARTI_isend(NTYPE1+proc, SmsgBuf[lPck], 
                 sizeS[lPck], i, NULL );
	      */
#ifdef PARTI_ENABLE_MP_INTERFACE_MPI
	      PARTI_isend(NTYPE1+proc,SmsgBuf[lPck],sizeS[lPck],i,&sGid[lPck]);
#else
	      sGid[lPck] = PARTI_isend(NTYPE1+proc,SmsgBuf[lPck],sizeS[lPck],i,NULL);
#endif
	    }
        }
      }
    }
  

      /* receive */
      if(receiving){
        for( lSrc=0; lSrc < GSsrc; lSrc++ ){
          i=lSrc+srcPrc;
          if (sched->rMsgSz[i] > 0) {
  	      if(i!=proc){
		/*
	        rGid[lSrc] = PARTI_irecv(NTYPE1+i, RmsgBuf[lSrc], 
	          sizeR[lSrc], i, NULL );
		*/
#ifdef PARTI_ENABLE_MP_INTERFACE_MPI
	        PARTI_irecv(NTYPE1+i, RmsgBuf[lSrc], sizeR[lSrc], i, &rGid[lSrc]);
#else
	        rGid[lSrc] = PARTI_irecv(NTYPE1+i, RmsgBuf[lSrc], sizeR[lSrc], i, 
					 NULL );
#endif
	      }
            }
          }
      }

      if(receiving){      
	allmesg = 0;
	for(j1=0;j1<grp_size;j1++)
		recv_done[j1]=0;
	while(!allmesg){
	  allmesg=1;
	  for( lWait=0; lWait < GSsrc; lWait++ ){
	    if( !recv_done[lWait] ){
	      i=lWait+srcPrc;
	      if (sched->rMsgSz[i] > 0) {
		if(i!=proc){
		  /*
		  if(PARTI_msgdone(rGid[lWait])){
		  */
#ifdef PARTI_ENABLE_MP_INTERFACE_MPI
		  if(PARTI_msgdone(&rGid[lWait])){
#else
		  if(PARTI_msgdone(rGid[lWait])){
#endif
		    recv_done[lWait] = 1;
		  }
		  else{
		    allmesg=0;
		  }
		}
		else
		  recv_done[lWait] = 1;
	        if( recv_done[lWait]){


      /* unpack receive */
		
            recvData = sched->rData[i];
        /*
         * unbuffer received data
        */
  
        /*  Gagan   07/03/93:  Rewritten to accomodate for
                               the new definition of schedule */

	    /* als 9/96 - added arbitrary dimensions for arrays */
  
            startPosn = recvData->startPosn ;
            destPtr   = dest + startPosn ;
  
            arraydims = recvData->numDims ;

            if (arraydims > MAX_DIM) {
              fatal_error(
      	  "Datamove into array with more than MAX_DIM dimensions not allowed");
            }

	    for (j = 0; j < arraydims; j++) {
	      step[j] = recvData->str[j];
	      numPieces[j] = recvData->numelem[j] ;
	    }

	    /* set lower bounds for all loops - 0 here */
	    for (currdim = 0; currdim < arraydims; currdim++) {
	      k[currdim] = 0;
	      tempPtr[currdim] = destPtr;
	    }

	    for (currdim = 0; ; currdim = 0) {
	      /* reset indices and bump and reset buffer pointers for loops
		 that are done */
	      while (k[currdim] >= numPieces[currdim]) {
		if (currdim == arraydims-1) goto done2;

		k[currdim] = 0;
		currdim++;
		k[currdim]++;
		tempPtr[currdim] += step[currdim];
		/* reset all inner loop buffer pointers */
		for (j = currdim-1; j >= 0; j--) {
		  tempPtr[j] = tempPtr[j+1];
		}
	      }
	      /* inner loop code */
	      *tempPtr[0] = *RbufPtr[lWait]  ;
	      RbufPtr[lWait]++ ;
	      tempPtr[0] += step[0] ;
	      
	      k[0]++;
	    }
done2:      ;

	        }
	      }
	    }
	  }
	}
      }

      if(sending){
        for( lWait=0; lWait < GSdest; lWait++ ){
          i=lWait+destPrc;
          if (sched->sMsgSz[i] > 0) {
            if(i!=proc){
	      /*
              PARTI_msgwait(sGid[lWait]);
	      */
#ifdef PARTI_ENABLE_MP_INTERFACE_MPI
              PARTI_msgwait(&sGid[lWait]);
#else
              PARTI_msgwait(sGid[lWait]);
#endif
            }
          }
	}
      }
    }
  }
}



void cDataMove(array, sched, dest)
  char            *array, *dest;
  SCHED    *sched;
{
/*
 * Nonblocking/Unbuffered version for the SP-1
 *	-Jim 6/94
 */

  int srcPrc, destPrc, lDest, lSrc, lPck, lWait, i,
      arraydims, sending, receiving, grp_size, GSsrc, GSdest, proc, 
      nproc, startPosn,Wprocs,allmesg,j1;

  int numPieces[MAX_DIM], k[MAX_DIM], step[MAX_DIM], currdim, j;

  int sGidVALID[GRPSIZE], rGidVALID[GRPSIZE], sizeR[GRPSIZE],
  sizeS[GRPSIZE], recv_done[GRPSIZE];

  char         *arrayPtr, *destPtr, *RbufPtr[GRPSIZE], *SbufPtr[GRPSIZE];

  char         *tempPtr[MAX_DIM];

  SchedData *recvData, *sendData;

#ifdef PARTI_ENABLE_MP_INTERFACE_MPI
  MPI_Request sGid[GRPSIZE], rGid[GRPSIZE];
#else
  long sGid[GRPSIZE], rGid[GRPSIZE];
#endif

  if (sched == NULL)
    return;

  proc = PARTI_myproc();
  nproc = PARTI_numprocs();

  if( GRPSIZE > nproc )
    grp_size = nproc;
  else
    grp_size = GRPSIZE;

  Wprocs = (nproc / grp_size) * grp_size;

  /* Break 0..(nproc-1) loop into groups of size grp_size */
  for( srcPrc=0; srcPrc < nproc; srcPrc+=grp_size){
      if( srcPrc >= Wprocs   )
	GSsrc = nproc - Wprocs;
      else
	GSsrc = grp_size;

      GSdest = grp_size;

      if( (proc<(srcPrc+GSsrc))&&(proc>=srcPrc) )
	sending=1;
      else
	sending=0;

    for( destPrc=0; destPrc < nproc; destPrc+=grp_size){

      if( destPrc >= Wprocs   )
	GSdest = nproc - Wprocs;


      if( (proc<(destPrc+GSdest))&&(proc>=destPrc) )
	receiving=1;
      else
	receiving=0;


      /* check receive allocation */
      if( receiving ){
        for( lPck=0; lPck < GSsrc; lPck++ ){
          i = srcPrc+lPck;

          if (sched->rMsgSz[i] > 0) { 
            sizeR[lPck] = sched->rMsgSz[i] * CHR_SZ;
      
            if (i != proc) {          /* als - 5/19/92 */
              if (sizeR[lPck] > RmsgBuf_size[lPck]) {
            /* need a bigger message buffer (als - 4/93) */
                if (RmsgBuf[lPck] != NULL) free(RmsgBuf[lPck]);
                RmsgBuf[lPck] = (char *) malloc(sizeR[lPck]);
                if (RmsgBuf[lPck] == NULL) {
                  fatal_error("can't allocate space for message buffer in dataMove");
                }
                RmsgBuf_size[lPck] = sizeR[lPck];
              }
                RbufPtr[lPck]  = (char *) RmsgBuf[lPck];

            } 
            else {
              /* copy from the local message buffer */
              if (sizeR[lPck] > local_msgBuf_size) {
                /* need a bigger local message buffer (als - 4/93) */
                if (local_msgBuf != NULL) free(local_msgBuf);
                local_msgBuf = (char *) malloc(sizeR[lPck]);
                if (local_msgBuf == NULL) {
                  fatal_error(
                    "can't allocate space for local message buffer in dataMove");
                }
                local_msgBuf_size = sizeR[lPck];
              }
              RbufPtr[lPck] = (char *) local_msgBuf;
           } 

         }
       }
     } 


      /* check sendbuffer allocation */
      if(sending){
        for( lPck=0; lPck < GSdest; lPck++ ){
	  i=lPck+destPrc;
          if (sched->sMsgSz[i] > 0) {
            sizeS[lPck] = sched->sMsgSz[i] * CHR_SZ;
            if (i != proc) {
            /* prepare a send */
              if (sizeS[lPck] > SmsgBuf_size[lPck]) {
                /* need a bigger message buffer (als - 4/93) */
                if (SmsgBuf[lPck] != NULL) free(SmsgBuf[lPck]);
                SmsgBuf[lPck] = (char *) malloc(sizeS[lPck]);
                if (SmsgBuf[lPck] == NULL) { fatal_error(
		  "can't allocate space for message buffer in dataMove"); }
                SmsgBuf_size[lPck] = sizeS[lPck];
              }
              SbufPtr[lPck] = (char *) SmsgBuf[lPck];
            }
            else {                    /* als - 5/19/92 */
            /* prepare the on-processor copy */
              if (sizeS[lPck] > local_msgBuf_size) {
              /* need a bigger local message buffer (als - 4/93) */
                if (local_msgBuf != NULL) free(local_msgBuf);
                local_msgBuf = (char *) malloc(sizeS[lPck]);
                if (local_msgBuf == NULL) {
                  fatal_error("can't allocate space for local message buffer in dataMove");
                }
                local_msgBuf_size = sizeS[lPck];
              }
              SbufPtr[lPck] = (char *) local_msgBuf;
            }

      /* pack sends */
            sendData = sched->sData[i];
        /* New Stuff goes here */
  
        /*  Gagan   07/03/93:  Rewritten to accomodate for
                               the new definition of schedule */
  
	    /* als 9/96 - added arbitrary dimensions for arrays */
  
            startPosn = sendData->startPosn ;
            arrayPtr  = array + startPosn ;
  
            arraydims = sendData->numDims ;

            if (arraydims > MAX_DIM) {
              fatal_error(
      	  "Datamove from array with more than MAX_DIM dimensions not allowed");
            }

	    for (j = 0; j < arraydims; j++) {
	      step[j] = sendData->str[j];
	      numPieces[j] = sendData->numelem[j] ;
	    }
	    /* set lower bounds for all loops - 0 here */
	    for (currdim = 0; currdim < arraydims; currdim++) {
	      k[currdim] = 0;
	      tempPtr[currdim] = arrayPtr;
	    }

	    for (currdim = 0; ; currdim = 0) {
	      /* reset indices and bump and reset buffer pointers for loops
		 that are done */
	      while (k[currdim] >= numPieces[currdim]) {
		if (currdim == arraydims-1) goto done1;

		k[currdim] = 0;
		currdim++;
		k[currdim]++;
		tempPtr[currdim] += step[currdim];
		/* reset all inner loop buffer pointers */
		for (j = currdim-1; j >= 0; j--) {
		  tempPtr[j] = tempPtr[j+1];
		}
	      }
	      /* inner loop code */
	      *SbufPtr[lPck] = *tempPtr[0] ;
	      SbufPtr[lPck]++ ;
	      tempPtr[0] += step[0] ;

	      k[0]++;

	    }
done1:      ;

            if(i != proc ){
	      /*
	      sGid[lPck] = PARTI_isend(NTYPE1+proc, SmsgBuf[lPck], 
                 sizeS[lPck], i, NULL );
	      */
#ifdef PARTI_ENABLE_MP_INTERFACE_MPI
	      PARTI_isend(NTYPE1+proc,SmsgBuf[lPck],sizeS[lPck],i,&sGid[lPck]);
#else
	      sGid[lPck] = PARTI_isend(NTYPE1+proc,SmsgBuf[lPck],sizeS[lPck],i,NULL);
#endif
	    }
        }
      }
    }
  

      /* receive */
      if(receiving){
        for( lSrc=0; lSrc < GSsrc; lSrc++ ){
          i=lSrc+srcPrc;
          if (sched->rMsgSz[i] > 0) {
  	      if(i!=proc){
		/*
	        rGid[lSrc] = PARTI_irecv(NTYPE1+i, RmsgBuf[lSrc], 
	          sizeR[lSrc], i, NULL );
		*/
#ifdef PARTI_ENABLE_MP_INTERFACE_MPI
	        PARTI_irecv(NTYPE1+i, RmsgBuf[lSrc], sizeR[lSrc], i, &rGid[lSrc]);
#else
	        rGid[lSrc] = PARTI_irecv(NTYPE1+i, RmsgBuf[lSrc], sizeR[lSrc], i, 
					 NULL );
#endif
	      }
            }
	}
      }

      if(receiving){      
	allmesg = 0;
	for(j1=0;j1<grp_size;j1++)
		recv_done[j1]=0;
	while(!allmesg){
	  allmesg=1;
	  for( lWait=0; lWait < GSsrc; lWait++ ){
	    if( !recv_done[lWait] ){
	      i=lWait+srcPrc;
	      if (sched->rMsgSz[i] > 0) {
		if(i!=proc){
		  /*
		  if(PARTI_msgdone(rGid[lWait])){
		  */
#ifdef PARTI_ENABLE_MP_INTERFACE_MPI
		  if(PARTI_msgdone(&rGid[lWait])){
#else
		  if(PARTI_msgdone(rGid[lWait])){
#endif
		    recv_done[lWait] = 1;
		  }
		  else{
		    allmesg=0;
		  }
		}
		else
		  recv_done[lWait] = 1;
	        if( recv_done[lWait]){


      /* unpack receive */
		
            recvData = sched->rData[i];
        /*
         * unbuffer received data
        */
  
        /*  Gagan   07/03/93:  Rewritten to accomodate for
                               the new definition of schedule */
  
	    /* als 9/96 - added arbitrary dimensions for arrays */
  
            startPosn = recvData->startPosn ;
            destPtr   = dest + startPosn ;
  
            arraydims = recvData->numDims ;

            if (arraydims > MAX_DIM) {
              fatal_error(
      	  "Datamove into array with more than MAX_DIM dimensions not allowed");
            }

	    for (j = 0; j < arraydims; j++) {
	      step[j] = recvData->str[j];
	      numPieces[j] = recvData->numelem[j] ;
	    }

	    /* set lower bounds for all loops - 0 here */
	    for (currdim = 0; currdim < arraydims; currdim++) {
	      k[currdim] = 0;
	      tempPtr[currdim] = destPtr;
	    }

	    for (currdim = 0; ; currdim = 0) {
	      /* reset indices and bump and reset buffer pointers for loops
		 that are done */
	      while (k[currdim] >= numPieces[currdim]) {
		if (currdim == arraydims-1) goto done2;

		k[currdim] = 0;
		currdim++;
		k[currdim]++;
		tempPtr[currdim] += step[currdim];
		/* reset all inner loop buffer pointers */
		for (j = currdim-1; j >= 0; j--) {
		  tempPtr[j] = tempPtr[j+1];
		}
	      }
	      /* inner loop code */
	      *tempPtr[0] = *RbufPtr[lWait]  ;
	      RbufPtr[lWait]++ ;
	      tempPtr[0] += step[0] ;
	      
	      k[0]++;
	    }
done2:      ;

	        }
	      }
	    }
	  }
	}
      }


      if(sending){
        for( lWait=0; lWait < GSdest; lWait++ ){
          i=lWait+destPrc;
          if (sched->sMsgSz[i] > 0) {
            if(i!=proc){
	      /*
              PARTI_msgwait(sGid[lWait]);
	      */
#ifdef PARTI_ENABLE_MP_INTERFACE_MPI
              PARTI_msgwait(&sGid[lWait]);
#else
              PARTI_msgwait(sGid[lWait]);
#endif
            }
          }
	}
      }
    }
  }
}

/****************************************************************
*
*  Frees all message buffers, resetting sizes to 0
*
*    added for Gopal (als - 3/96)
*
*****************************************************************/
void free_message_buffers()
{
  int i;

  for (i = 0; i < GRPSIZE; i++) {
    if (RmsgBuf[i] != NULL) {
      free(RmsgBuf[i]);
      RmsgBuf[i] = NULL;
      RmsgBuf_size[i] = 0;
    }
    if (SmsgBuf[i] != NULL) {
      free(SmsgBuf[i]);
      SmsgBuf[i] = NULL;
      SmsgBuf_size[i] = 0;
    }
  }

  if (local_msgBuf != NULL) {
    free(local_msgBuf);
    local_msgBuf = NULL;
    local_msgBuf_size = 0;
  }
}


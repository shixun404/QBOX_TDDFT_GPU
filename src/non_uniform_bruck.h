/*
 * non_uniform_bruck.h
 *
 *      Author: kokofan
 */

#ifndef SRC_NON_UNIFORM_BRUCK_H_
#define SRC_NON_UNIFORM_BRUCK_H_

#include "brucks.h"


void twophase_bruck_alltoallv(char *sendbuf, int *sendcounts, int *sdispls, MPI_Datatype sendtype, char *recvbuf, int *recvcounts, int *rdispls, MPI_Datatype recvtype, MPI_Comm comm);

void twophase_bruck_alltoallv_new(char *sendbuf, int *sendcounts, int *sdispls, MPI_Datatype sendtype, char *recvbuf, int *recvcounts, int *rdispls, MPI_Datatype recvtype, MPI_Comm comm);


#endif /* SRC_NON_UNIFORM_BRUCK_H_ */

from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

if size > 4:
    exit('Hai scelto troppi processi')

irecv = np.zeros(size, dtype=np.int)
isend = np.zeros(1,  dtype=np.int)
isend[0] = rank+1

comm.Gather(isend, irecv, root=0)

if rank == 0:
    print('irecv: ', irecv)
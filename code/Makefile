target: decide

CC=mpicc
CFLAGS=-g -fopenmp -O3 -std=c99

main.o: main.c projet.h
aux.o: aux.c projet.h

decide: main.o aux.o
	$(CC) $(CFLAGS) $^ -o $@

.PHONY: clean exec

clean:
	rm -f *.o decide

exec:
	export OMP_NUM_THREADS=2
	mpirun -n 2 -hostfile hostfile --bynode ./decide "4k//4K/4P w"
#	mpirun -n 2 -hostfile hostfile --bynode ./decide "///2kpK/7P w" >> ~/Bureau/MPI/results_project/MPI+OpenMP/2.txt
#	mpirun -n 2 -hostfile hostfile --bynode ./decide "4k//4K/4P w" >> ~/Bureau/MPI/results_project/MPI+OpenMP/3.txt
#	mpirun -n 2 -hostfile hostfile --bynode ./decide "4k//4K//4P w" >> ~/Bureau/MPI/results_project/MPI+OpenMP/4.txt
#	mpirun -n 2 -hostfile hostfile --bynode ./decide "/ppp//PPP//7k//7K w" >> ~/Bureau/MPI/results_project/MPI+OpenMP/5.txt
#	mpirun -n 2 -hostfile hostfile --bynode ./decide "4k//4K///4P w" >> ~/Bureau/MPI/results_project/MPI+OpenMP/6.txt
#	mpirun -n 2 -hostfile hostfile --bynode ./decide "7K//k1P/7p b" >> ~/Bureau/MPI/results_project/MPI+OpenMP/7.txt

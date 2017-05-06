target: decide

CFLAGS=-g -O3 -std=c99 -fopenmp

main.o: main.c projet.h
aux.o: aux.c projet.h


decide: main.o aux.o
	$(CC) $(CFLAGS) $^ -o $@

.PHONY: clean exec

clean:
	rm -f *.o decide

exec:
	export OMP_NUM_THREADS=4
	./decide "4k//4K/4P w"
#	./decide "///2kpK/7P w" >> ~/Bureau/MPI/results_project/OpenMP/2.txt
#	./decide "4k//4K/4P w" >> ~/Bureau/MPI/results_project/OpenMP/3.txt
#	./decide "4k//4K//4P w" >> ~/Bureau/MPI/results_project/OpenMP/4.txt
#	./decide "/ppp//PPP//7k//7K w" >> ~/Bureau/MPI/results_project/OpenMP/5.txt
#	./decide "4k//4K///4P w" >> ~/Bureau/MPI/results_project/OpenMP/6.txt
#	./decide "7K//k1P/7p b" >> ~/Bureau/MPI/results_project/OpenMP/7.txt

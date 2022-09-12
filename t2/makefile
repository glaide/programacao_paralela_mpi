CC = mpicc
#flags para warning, depuração e biblioteca de matematica
CFLAGS = -Wall --std=c99 -L. -lm -O3

#regras de ligação
lcs = lcs_mpi.o

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) 

#variavel geral
all: lcs_mpi
#regras de compilação
lcs_mpi: $(lcs)
	$(CC) -o $@ $^ $(CFLAGS) 

.PHONY:	clean
#remove os arquivos temporarios e executaveis
clean: 

	rm *.o lcs_mpi 
	
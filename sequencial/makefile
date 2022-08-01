
CC=gcc
#flags para warning, depuração e biblioteca de matematica
CFLAGS= -Wall -fopenmp -O3

#regras de ligação
lcs = lcs-sequencial.o

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) 

#variavel geral
all: lcs_seq
#regras de compilação
lcs_seq: $(lcs)
	$(CC) -o $@ $^ $(CFLAGS) 

.PHONY:	clean
#remove os arquivos temporarios e executaveis
clean: 

	rm *.o lcs_seq 
	rm -rf ./object

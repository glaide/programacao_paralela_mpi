
CC=gcc
#flags para warning, depuração e biblioteca de matematica
CFLAGS= -Wall -fopenmp -O3

#regras de ligação
lcs = find_lcs.o

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) 

#variavel geral
all: lcs_par
#regras de compilação
lcs_par: $(lcs)
	$(CC) -o $@ $^ $(CFLAGS) 

.PHONY:	clean
#remove os arquivos temporarios e executaveis
clean: 

	rm *.o lcs_par 
	rm -rf ./object

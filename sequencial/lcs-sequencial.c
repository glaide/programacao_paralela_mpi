#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef max
#define max(a, b) (((a) > (b)) ? (a) : (b))
#endif

typedef struct sequencia
{
	char *texto;
	long tam;
} t_sequencia;

t_sequencia ler_entrada(char *filename);

typedef unsigned short mtype;

/* Read sequence from a file to a char vector.
 Filename is passed as parameter */

t_sequencia ler_entrada(char *filename)
{
	FILE *arquivo_sequencia = NULL;
	int i = 0;
	arquivo_sequencia = fopen(filename, "rt");
	t_sequencia s;

	if (arquivo_sequencia == NULL)
	{
		printf("a reading file %s\n", filename);
		exit(1);
	}

	// verifica o tamanho da arquivo_sequencia
	fseek(arquivo_sequencia, 0L, SEEK_END);
	// salva o tamanho da arquivo_sequencia
	s.tam = ftell(arquivo_sequencia);
	// volta o arquivo para o inicio
	rewind(arquivo_sequencia);

	// aloca memoria para salvar valor da sequencia

	s.texto = (char *)malloc(s.tam + 1 * sizeof(char));
	if (s.texto == NULL)
	{
		printf("Erro ao alocar memória da sequencia");
		exit(1);
	}

	// leitura da sequencia
	// read sequence from file
	while (!feof(arquivo_sequencia))
	{
		s.texto[i] = fgetc(arquivo_sequencia);
		if ((s.texto[i] != '\n') && (s.texto[i] != EOF))
			i++;
	}

	strcpy("\0", &s.texto[i]);

	// fecha o arquivo
	fclose(arquivo_sequencia);
	return s;
}
mtype **allocateScoreMatrix(long sizeA, long sizeB)
{
	int i;
	// Allocate memory for LCS score matrix
	mtype **scoreMatrix = (mtype **)malloc((sizeB + 1) * sizeof(mtype *));
	for (i = 0; i < (sizeB + 1); i++)
		scoreMatrix[i] = (mtype *)malloc((sizeA + 1) * sizeof(mtype));
	return scoreMatrix;
}

void initScoreMatrix(mtype **scoreMatrix, long sizeA, long sizeB)
{
	int i, j;
	// Fill first line of LCS score matrix with zeroes
	for (j = 0; j < (sizeA + 1); j++)
		scoreMatrix[0][j] = 0;

	// Do the same for the first collumn
	for (i = 1; i < (sizeB + 1); i++)
		scoreMatrix[i][0] = 0;
}

int LCS(mtype **scoreMatrix, long sizeA, long sizeB, char *seqA, char *seqB)
{
	int i, j;
	for (i = 1; i < sizeB + 1; i++)
	{
		for (j = 1; j < sizeA + 1; j++)
		{
			if (seqA[j - 1] == seqB[i - 1])
			{
				/* if elements in both sequences match,
				 the corresponding score will be the score from
				 previous elements + 1*/
				scoreMatrix[i][j] = scoreMatrix[i - 1][j - 1] + 1;
			}
			else
			{
				/* else, pick the maximum value (score) from left and upper elements*/
				scoreMatrix[i][j] =
					max(scoreMatrix[i - 1][j], scoreMatrix[i][j - 1]);
			}
		}
	}
	return scoreMatrix[sizeB][sizeA];
}
void printMatrix(char *seqA, char *seqB, mtype **scoreMatrix, int sizeA,
				 int sizeB)
{
	int i, j;

	// print header
	printf("Score Matrix:\n");
	printf("========================================\n");

	// print LCS score matrix allong with sequences

	printf("    ");
	printf("%5c   ", ' ');

	for (j = 0; j < sizeA; j++)
		printf("%5c   ", seqA[j]);
	printf("\n");
	for (i = 0; i < sizeB + 1; i++)
	{
		if (i == 0)
			printf("    ");
		else
			printf("%c   ", seqB[i - 1]);
		for (j = 0; j < sizeA + 1; j++)
		{
			printf("%5d   ", scoreMatrix[i][j]);
		}
		printf("\n");
	}
	printf("========================================\n");
}

void freeScoreMatrix(mtype **scoreMatrix, long sizeB)
{
	int i;
	for (i = 0; i < (sizeB + 1); i++)
		free(scoreMatrix[i]);
	free(scoreMatrix);
}

int main(int argc, char **argv)
{

	t_sequencia sequencia1;
	t_sequencia sequencia2;

	for (int i = 0; i < argc; i++)
	{
		printf("%s\n\n", argv[i]);
	}

	printf("aaaaaaaaaaa\n\n");

	// read both sequences
	sequencia1 = ler_entrada(argv[1]);
	sequencia2 = ler_entrada(argv[2]);

	// 	// allocate LCS score matrix
	// 	mtype **scoreMatrix = allocateScoreMatrix(sequencia1.tam, sequencia2.tam);

	// 	// initialize LCS score matrix
	// 	initScoreMatrix(scoreMatrix, sequencia1.tam, sequencia2.tam);

	// 	// fill up the rest of the matrix and return final score (element locate at the last line and collumn)
	// 	mtype score = LCS(scoreMatrix, sequencia1.tam, sequencia2.tam, sequencia1.texto, sequencia2.texto);

	// 	// 	/* if you wish to see the entire score matrix,
	// 	// 	 for debug purposes, define DEBUGMATRIX. */
	// #ifdef DEBUGMATRIX
	// 	printMatrix(seqA, seqB, scoreMatrix, sequencia1.tam, sequencia2.tam);
	// #endif

	// 	// print score
	// 	printf("\nScore: %d\n", score);

	// 	// free score matrix
	// 	freeScoreMatrix(scoreMatrix, sequencia2.tam);

	// 	return EXIT_SUCCESS;
}

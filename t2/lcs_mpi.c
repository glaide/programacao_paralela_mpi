#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"
#include <time.h>
// macros
#define TAM_C_UNICOS 500
#define max(x, y) ((x) > (y) ? (x) : (y))

// function prototypes
int indice(char *str, char x, int len);
void print_matrix(int *x, int row, int col);
int lcs(int *atual, int *anterior, int *P, char *A, char *B, char *C, int len_a, int len_b, int len_c, int myrank, int chunk_size, int resto);
void inicia_matriz_p(int *P, char *b, int len_b, char *c, int len_c, int myrank, int chunk_size, int resto);

char *string_char_unicos(char *seqA, char *seqB)
{
    // Aloca um array com todos os possíveis chars que possamos ter de entrada (2⁸)
    int *charsetC = calloc(TAM_C_UNICOS, sizeof(int));
    int sizeC = 0;
    // leio a seqB e separo seu chars únicos escrevendo 1 em sua posição de mem. equivalente em charsetC
    for (int i = 0; i < strlen(seqA); i++)
    {
        if (charsetC[(int)seqA[i]] == 0)
        {
            sizeC += 1; // contabiliza o num. de chars diferentes existentes
        }
        charsetC[(int)seqA[i]] = 1;
    }
    // faço o mesmo para a seqB
    for (int i = 0; i < strlen(seqB); i++)
    {
        if (charsetC[(int)seqB[i]] == 0)
        {
            sizeC += 1;
        }
        charsetC[(int)seqB[i]] = 1;
    }

    // aloco espaço para a string C e coloco seus valores (será sempre < 256)
    char *seqC = calloc(sizeC, sizeof(char));
    int j = 0;
    for (int i = 0; i < TAM_C_UNICOS; i++)
    {
        if (charsetC[i] == 1)
        {
            seqC[j] = i;
            j++;
        }
    }

    return seqC;
}

char *read_seq(char *fname)
{
    // file pointer
    FILE *fseq = NULL;
    // sequence size
    long size = 0;
    // sequence pointer
    char *seq = NULL;
    // sequence index
    int i = 0;

    // open file
    fseq = fopen(fname, "rt");
    if (fseq == NULL)
    {
        printf("Error reading file %s\n", fname);
        exit(1);
    }

    // find out sequence size to allocate memory afterwards
    fseek(fseq, 0, SEEK_END);
    size = ftell(fseq);
    rewind(fseq);

    // allocate memory (sequence)
    seq = (char *)calloc(size + 1, sizeof(char));
    if (seq == NULL)
    {
        printf("Erro allocating memory for sequence %s.\n", fname);
        exit(1);
    }

    // read sequence from file
    while (!feof(fseq))
    {
        seq[i] = fgetc(fseq);
        if ((seq[i] != '\n') && (seq[i] != EOF))
            i++;
    }
    // insert string terminator
    seq[i] = '\0';

    // close file
    fclose(fseq);

    // return sequence pointer
    return seq;
}

int indice(char *str, char x, int len)
{
    for (int i = 0; i < len; i++)
    {
        if (str[i] == x)
        {
            return i;
        }
    }
    return -1; // not found the character x in str
}

void print_matrix(int *x, int row, int col)
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            printf("%d ", x[i * col + j]);
        }
        printf("\n");
    }
}

void inicia_matriz_p(int *P, char *b, int len_b, char *c, int len_c, int myrank, int chunk_size, int resto)
{
    char string_c_scatter[chunk_size];
    int scatter[chunk_size * (len_b + 1)];

    MPI_Scatter(c, chunk_size, MPI_CHAR, &string_c_scatter, chunk_size, MPI_CHAR, 0, MPI_COMM_WORLD);

    MPI_Scatter(P, chunk_size * (len_b + 1), MPI_INT, &scatter, chunk_size * (len_b + 1), MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(b, len_b, MPI_CHAR, 0, MPI_COMM_WORLD);

    int i, j;
    for (i = 0; i < chunk_size; i++)
    {
        for (j = 1; j < len_b + 1; j++)
        {
            if (b[j - 1] == string_c_scatter[i])
            {
                scatter[i * (len_b + 1) + j] = j;
            }
            else
            {
                scatter[i * (len_b + 1) + j] = scatter[i * (len_b + 1) + j - 1];
            }
        }
    }

    MPI_Gather(scatter, chunk_size * (len_b + 1), MPI_INT, P, chunk_size * (len_b + 1), MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    if (myrank == 0)
    {
        for (i = len_c - resto; i < len_c; i++)
        {
            for (j = 1; j < len_b + 1; j++)
            {
                if (b[j - 1] == c[i])
                {
                    P[i * (len_b + 1) + j] = j;
                }
                else
                {
                    P[i * (len_b + 1) + j] = P[i * (len_b + 1) + j - 1];
                }
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

int lcs(int *atual, int *anterior, int *P, char *A, char *B, char *C, int len_a, int len_b, int len_c, int myrank, int chunk_size, int resto)
{
    MPI_Bcast(P, (len_c * (len_b + 1)), MPI_INT, 0, MPI_COMM_WORLD);
    int *aux;
    int i;

    int start_id = (myrank * chunk_size);
    int end_id = (myrank * chunk_size) + chunk_size;

    int dp_i_receive[chunk_size];

    for (i = 1; i < len_a; i++)
    {
        int c = indice(C, A[i - 1], len_c);
        int p_c_j, j;

        MPI_Scatter(atual, chunk_size, MPI_INT, dp_i_receive, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);

        j = start_id;

        if (j == 0)
            j = 1;

        for (; j < end_id; j++)
        {
            p_c_j = P[c * (len_b + 1) + j];
            if (p_c_j)
            {
                dp_i_receive[j - start_id] = max(anterior[j], anterior[p_c_j - 1] + 1);
            }
            else
            {
                dp_i_receive[j - start_id] = anterior[j];
            }
        }

        MPI_Allgather(dp_i_receive, chunk_size, MPI_INT, atual, chunk_size, MPI_INT, MPI_COMM_WORLD);

        if (myrank == 0)
        {
            for (j = len_b + 1 - resto; j < len_b + 1; j++)
            {
                p_c_j = P[c * (len_b + 1) + j];
                if (p_c_j)
                {
                    atual[j] = max(anterior[j], anterior[p_c_j - 1] + 1);
                }
                else
                {

                    atual[j] = anterior[j];
                }
            }
        }

        aux = atual;
        atual = anterior;
        anterior = aux;

        MPI_Barrier(MPI_COMM_WORLD);
    }
    // o algoritmo faz uma troca entre as linhas a mais que devia na última iteração, então a corrigimos
    aux = atual;
    atual = anterior;
    anterior = aux;

    return atual[len_b];
}

int main(int argc, char *argv[])
{
    if (argc <= 1)
    {
        printf("Error: No input file specified! Please specify the input file, and run again!\n");
        return 0;
    }

    int my_rank, num_procs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);   // grab this process's rank
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs); // grab the total num of processes

    int len_a, len_b, len_c;
    double start_time, stop_time, start_time_yang, stop_time_yang;
    char *string_A = read_seq(argv[1]);
    char *string_B = read_seq(argv[2]);

    char *unique_chars_C = string_char_unicos(string_A, string_B);
    len_a = strlen(string_A);
    len_b = strlen(string_B);
    len_c = strlen(unique_chars_C);
    len_c++;

    int chunk_size_p = (len_c / num_procs);
    int resto_p = (len_c % num_procs);
    int chunk_size_dp = ((len_b + 1) / num_procs);
    int resto_dp = ((len_b + 1) % num_procs);

    if (my_rank == 0)
    {
        printf("chunk_p: %d chunk_dp: %d procs: %d\n", chunk_size_p, chunk_size_dp, num_procs);
    }

    int *atual = calloc((len_b + 1), sizeof(int));
    int *anterior = calloc((len_b + 1), sizeof(int));

    int *P_Matrix = calloc((len_c * (len_b + 1)), sizeof(int));

    start_time_yang = MPI_Wtime();
    inicia_matriz_p(P_Matrix, string_B, len_b, unique_chars_C, len_c, my_rank, chunk_size_p, resto_p);
    print_matrix(P_Matrix, len_c, len_b + 1);
    int res = lcs(atual, anterior, P_Matrix, string_A, string_B, unique_chars_C, len_a, len_b, len_c, my_rank, chunk_size_dp, resto_dp);

    stop_time_yang = MPI_Wtime();

    if (my_rank == 0)
    {
        printf("lcs is: %d\n", res);
        printf("time taken for lcs is: %lf\n", stop_time_yang - start_time_yang);
    }
    // deallocate pointers
    free(atual);
    free(anterior);
    free(P_Matrix);

    // Shutdown MPI (important - don't forget!)
    MPI_Finalize();
    return 0;
}

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"
#include <time.h>
// macros
#define TAM_MAX_CHAR 256 // tabela ascii
#define max(x, y) ((x) > (y) ? (x) : (y))
typedef struct sequencia
{
    char *texto;
    long tam;
} t_sequencia;

// function prototypes
int indice(char *str, char x, int len);
void print_matrix(int *x, int row, int col);
int lcs(int *atual, int *anterior, int *P, t_sequencia seqA, t_sequencia seqB, t_sequencia seqC, int myrank, int chunk_size, int resto);
void inicia_matriz_p(int *P, int myrank, t_sequencia seqB, t_sequencia seqC, int chunk_size, int resto);
t_sequencia calc_char_unicos(t_sequencia sequenciaA, t_sequencia sequenciaB);
t_sequencia ler_entrada(char *filename);

t_sequencia calc_char_unicos(t_sequencia seqA, t_sequencia seqB)
{
    int tam_char = 0;
    int index, index2 = 0;
    int *seq_char = calloc(TAM_MAX_CHAR, sizeof(int));

    t_sequencia s;

    for (int i = 0; i > seqA.tam; i++)
    {
        index = (int)seqA.texto[i];
        if (seq_char[index] == 0) // significa que ainda não passou por esse char
        {
            tam_char++;
        }
        seq_char[index] = 1;
    }

    for (int i = 0; i < seqB.tam; i++)
    {
        index2 = (int)seqB.texto[i];
        if (seq_char[index2] == 0) // significa que ainda não passou por esse char
        {
            tam_char++;
        }
        seq_char[index2] = 1;
    }

    s.texto = calloc(tam_char + 1, sizeof(char));
    s.tam = tam_char;

    int j = 0;
    for (int i = 0; i < TAM_MAX_CHAR + 1; i++)
    {
        if (seq_char[i] == 1)
        {
            s.texto[j] = i;
            j++;
        }
    }

    return s;
}

t_sequencia ler_entrada(char *filename)
{
    FILE *arquivo_sequencia = NULL;
    int i = 0;
    arquivo_sequencia = fopen(filename, "rt");
    t_sequencia s;

    if (arquivo_sequencia == NULL)
    {
        printf("Erro ao ler o arquivo: %s\n", filename);
        exit(1);
    }

    // verifica o tamanho da arquivo_sequencia
    fseek(arquivo_sequencia, 0L, SEEK_END);
    // salva o tamanho da arquivo_sequencia
    s.tam = ftell(arquivo_sequencia);
    // volta o arquivo para o inicio
    rewind(arquivo_sequencia);

    // aloca memoria para salvar valor da sequencia

    s.texto = (char *)calloc(s.tam + 1, sizeof(char));
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

    // fecha o arquivo
    fclose(arquivo_sequencia);
    return s;
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
    return 0; // not found the character x in str
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

void inicia_matriz_p(int *P, int myrank, t_sequencia seqB, t_sequencia seqC, int chunk_size, int resto)
{
    char string_c_scatter[chunk_size];
    int scatter[chunk_size * (seqB.tam + 1)];

    MPI_Scatter(seqC.texto, chunk_size, MPI_CHAR, &string_c_scatter, chunk_size, MPI_CHAR, 0, MPI_COMM_WORLD);

    MPI_Scatter(P, chunk_size * (seqB.tam + 1), MPI_INT, &scatter, chunk_size * (seqB.tam + 1), MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(seqB.texto, seqB.tam, MPI_CHAR, 0, MPI_COMM_WORLD);

    int i, j;
    for (i = 0; i < chunk_size; i++)
    {
        for (j = 1; j < seqB.tam + 1; j++)
        {
            if (seqB.texto[j - 1] == string_c_scatter[i])
            {
                scatter[i * (seqB.tam + 1) + j] = j;
            }
            else
            {
                scatter[i * (seqB.tam + 1) + j] = scatter[i * (seqB.tam + 1) + j - 1];
            }
        }
    }

    MPI_Gather(scatter, chunk_size * (seqB.tam + 1), MPI_INT, P, chunk_size * (seqB.tam + 1), MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    if (myrank == 0)
    {
        for (i = seqC.tam - resto; i < seqC.tam; i++)
        {
            for (j = 1; j < seqB.tam + 1; j++)
            {
                if (seqB.texto[j - 1] == seqC.texto[i])
                {
                    P[i * (seqB.tam + 1) + j] = j;
                }
                else
                {
                    P[i * (seqB.tam + 1) + j] = P[i * (seqB.tam + 1) + j - 1];
                }
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}
int lcs(int *atual, int *anterior, int *P, t_sequencia seqA, t_sequencia seqB, t_sequencia seqC, int myrank, int chunk_size, int resto)

{
    MPI_Bcast(P, (seqC.tam * (seqB.tam + 1)), MPI_INT, 0, MPI_COMM_WORLD);
    int *aux;
    int i;

    int start_id = (myrank * chunk_size);
    int end_id = (myrank * chunk_size) + chunk_size;

    int dp_i_receive[chunk_size];

    for (i = 1; i < seqA.tam; i++)
    {
        int c = indice(seqC.texto, seqA.texto[i - 1], seqC.tam);
        int p_c_j, j;

        MPI_Scatter(atual, chunk_size, MPI_INT, dp_i_receive, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);

        j = start_id;

        if (j == 0)
            j = 1;

        for (; j < end_id; j++)
        {
            p_c_j = P[c * (seqB.tam + 1) + j];
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
            for (j = seqB.tam + 1 - resto; j < seqB.tam + 1; j++)
            {
                p_c_j = P[c * (seqB.tam + 1) + j];
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

    return atual[seqB.tam];
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

    double start_time, stop_time;
    t_sequencia sequenciaA;
    t_sequencia sequenciaB;
    t_sequencia sequenciaC;

    sequenciaA = ler_entrada(argv[1]);
    sequenciaB = ler_entrada(argv[2]);
    sequenciaC = calc_char_unicos(sequenciaA, sequenciaB);

    int chunk_size_p = (sequenciaC.tam / num_procs);
    int resto_p = (sequenciaC.tam % num_procs);
    int chunk_size_dp = ((sequenciaB.tam) / num_procs);
    int resto_dp = ((sequenciaB.tam) % num_procs);

    if (my_rank == 0)
    {
        printf("chunk_p: %d chunk_dp: %d procs: %d\n", chunk_size_p, chunk_size_dp, num_procs);
    }

    int *atual = calloc((sequenciaB.tam + 1), sizeof(int));
    int *anterior = calloc((sequenciaB.tam + 1), sizeof(int));

    int *P_Matrix = calloc((sequenciaC.tam * (sequenciaB.tam + 1)), sizeof(int));

    printf("sequence: %i\n", sequenciaC.tam);
    printf("sequence: %s\n", sequenciaC.texto);
    start_time = MPI_Wtime();
    inicia_matriz_p(P_Matrix, my_rank, sequenciaB, sequenciaC, chunk_size_p, resto_p);
    print_matrix(P_Matrix, sequenciaC.tam, sequenciaB.tam + 1);
    int res = lcs(atual, anterior, P_Matrix, sequenciaA, sequenciaB, sequenciaC, my_rank, chunk_size_dp, resto_dp);

    stop_time = MPI_Wtime();

    if (my_rank == 0)
    {
        printf("lcs is: %d\n", res);
        printf("time taken for lcs is: %lf\n", stop_time - start_time);
    }
    // deallocate pointers
    free(atual);
    free(anterior);
    free(P_Matrix);

    // Shutdown MPI (important - don't forget!)
    MPI_Finalize();
    return 0;
}

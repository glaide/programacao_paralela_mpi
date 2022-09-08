#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"
#include <time.h>
// macros
#define ALPHABET_LENGTH 4
#define max(x, y) ((x) > (y) ? (x) : (y))

// global variables
char *string_A;
char *string_B;
char *unique_chars_C; // unique alphabets
int c_len;
int *P_Matrix;
int *atual; // to store the DP values
int *anterior;
// function prototypes
int get_index_of_character(char *str, char x, int len);
void print_matrix(int **x, int row, int col);
void calc_matriz_dist(int *P, char *b, int len_b, char *c, int len_c, int myrank, int chunk_size, int resto_p, int num_procs);
int calc_lcs(int *DP, int *prev_row, int *P, char *A, char *B, char *C, int m, int n, int u, int myrank, int chunk_size);
int lcs(int **DP, char *A, char *B, int m, int n);

int get_index_of_character(char *str, char x, int len)
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

void print_matrix(int **x, int row, int col)
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            printf("%d ", x[i][j]);
        }
        printf("\n");
    }
}

void calc_matriz_dist(int *P, char *b, int len_b, char *c, int len_c, int myrank, int chunk_size, int resto_p, int num_procs)
{
    char receive_array_for_scatter_c[chunk_size];
    int receive_array_for_scatter_p[chunk_size * (len_b + 1)];
    // Scatter the char array chunks by sending each process a particular chunk
    MPI_Scatter(c, chunk_size, MPI_CHAR, &receive_array_for_scatter_c, chunk_size, MPI_CHAR, 0, MPI_COMM_WORLD);
    // Scatter the char array chunks by sending each process a particular chunk
    MPI_Scatter(P, chunk_size * (len_b + 1), MPI_INT, &receive_array_for_scatter_p, chunk_size * (len_b + 1), MPI_INT, 0, MPI_COMM_WORLD);
    // Broadcast the whole b  array to everybody
    MPI_Bcast(b, len_b, MPI_CHAR, 0, MPI_COMM_WORLD);

    for (int i = 0; i < chunk_size; i++)
    {
        for (int j = 1; j < len_b + 1; j++)
        {
            if (b[j - 1] == receive_array_for_scatter_c[i])
            {
                receive_array_for_scatter_p[(i * (len_b + 1)) + j] = j;
            }
            else
            {
                receive_array_for_scatter_p[(i * (len_b + 1)) + j] = receive_array_for_scatter_p[(i * (len_b + 1)) + j - 1];
            }
        }
    }
    // now gather all the calculated values of P matrix in process 0
    MPI_Gather(receive_array_for_scatter_p, chunk_size * (len_b + 1), MPI_INT, P, chunk_size * (len_b + 1), MPI_INT, 0, MPI_COMM_WORLD);
}

int calc_lcs(int *DP, int *prev_row, int *P, char *A, char *B, char *C, int m, int n, int u, int myrank, int chunk_size)
{
    MPI_Bcast(P, (u * (n + 1)), MPI_INT, 0, MPI_COMM_WORLD);
    for (int i = 1; i < m + 1; i++)
    {
        int c_i = get_index_of_character(C, A[i - 1], u);
        int dp_i_receive[chunk_size];
        MPI_Scatter(DP, chunk_size, MPI_INT, &dp_i_receive, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);
        int start_id = (myrank * chunk_size);
        int end_id = (myrank * chunk_size) + chunk_size;

        int t, s;

        for (int j = start_id; j < end_id; j++)
        {
            if (j == start_id && myrank == 0)
                j = j + 1;
            t = (0 - P[(c_i * (n + 1)) + j]) < 0;
            s = (0 - (prev_row[j] - (t * prev_row[P[(c_i * (n + 1)) + j] - 1])));
            dp_i_receive[j - start_id] = ((t ^ 1) || (s ^ 0)) * (prev_row[j]) + (!((t ^ 1) || (s ^ 0))) * (prev_row[P[(c_i * (n + 1)) + j] - 1] + 1);
        }
        // now gather all the calculated values of P matrix in process 0
        MPI_Allgather(dp_i_receive, chunk_size, MPI_INT, DP, chunk_size, MPI_INT, MPI_COMM_WORLD);

        for (int j = 1; j < n + 1; j++)
        {
            prev_row[j] = DP[j];
        }
    }
    return DP[n];
}

int lcs(int **DP, char *A, char *B, int m, int n)
{

    for (int i = 1; i < (m + 1); i++)
    {
        for (int j = 1; j < (n + 1); j++)
        {

            if (A[i - 1] == B[j - 1])
            {
                DP[i][j] = DP[i - 1][j - 1] + 1;
            }
            else
            {
                DP[i][j] = max(DP[i - 1][j], DP[i][j - 1]);
            }
        }
    }

    return DP[m][n];
}

int main(int argc, char *argv[])
{
    if (argc <= 1)
    {
        printf("Error: No input file specified! Please specify the input file, and run again!\n");
        return 0;
    }

    int my_rank;
    int num_procs;
    int chunk_size_p, chunk_size_dp, resto_p; // chunk_size for P matrix and DP matrix
    int res;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);   // grab this process's rank
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs); // grab the total num of processes

    FILE *fp;
    int len_a, len_b;
    double start_time, stop_time, start_time_yang, stop_time_yang;

    if (my_rank == 0)
        printf("\nYour input file: %s \n", argv[1]);
    fp = fopen(argv[1], "r");
    fscanf(fp, "%d %d %d", &len_a, &len_b, &c_len);
    string_A = (char *)calloc((len_a + 1), sizeof(char *));
    string_B = (char *)calloc((len_b + 1), sizeof(char *));
    unique_chars_C = (char *)calloc((c_len + 1), sizeof(char *));

    fscanf(fp, "%s %s %s", string_A, string_B, unique_chars_C);

    chunk_size_p = (c_len / num_procs);
    resto_p = (c_len % num_procs);
    chunk_size_dp = ((len_b + 1) / num_procs);

    if (my_rank == 0)
    {
        printf("chunk_p: %d chunk_dp: %d procs: %d\n", chunk_size_p, chunk_size_dp, num_procs);
    }

    atual = (int *)calloc((len_b + 1), sizeof(int));
    anterior = (int *)calloc((len_b + 1), sizeof(int));

    P_Matrix = (int *)calloc((c_len * (len_b + 1)), sizeof(int));

    start_time_yang = MPI_Wtime();

    calc_matriz_dist(P_Matrix, string_B, len_b, unique_chars_C, c_len, my_rank, chunk_size_p, resto_p, num_procs);

    res = calc_lcs(atual, anterior, P_Matrix, string_A, string_B, unique_chars_C, len_a, len_b, c_len, my_rank, chunk_size_dp);

    stop_time_yang = MPI_Wtime();

    if (my_rank == 0)
    {
        printf("calc_lcs is: %d\n", res);
        printf("time taken for calc_lcs is: %lf\n", stop_time_yang - start_time_yang);
    }

    free(P_Matrix);
    free(atual);

    MPI_Finalize();
    return 0;
}
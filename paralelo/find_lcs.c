
/****
    Author: Rayhan Shikder,
    email: shikderr@myumanitoba.ca
    MSc Student,
    Department of Computer Science,
    University of Manitoba, Winnipeg, MB, Canada
****/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "omp.h"

// macros
#define TAM_MAX_CHAR 256

#define max(x, y) ((x) > (y) ? (x) : (y))

// global variables
char *string_A;
char *string_B;
char *unique_chars_C; // unique alphabets
int c_len;
int **P_Matrix;
int *DP_Results;  // to store the DP values of current row
int *dp_prev_row; // to store the DP values of previous row
typedef struct sequencia
{
    char *texto;
    int tam;
} t_sequencia;

t_sequencia ler_entrada(char *filename);

t_sequencia tam_char_seq(char *seqA, char *seqB, int *seq_char)
{
    int tam_char = 0;
    int index = 0;
    t_sequencia s;
    for (int i = 0; i < strlen(seqA); i++)
    {
        if (seq_char[(int)seqA[i]] == 0) // significa que ainda não passou por esse char
        {
            tam_char++;
        }
        index = (int)seqA[i];
        seq_char[index] = 1;
    }

    for (int i = 0; i < strlen(seqB); i++)
    {
        if (seq_char[(int)seqB[i]] == 0) // significa que ainda não passou por esse char
        {
            tam_char++;
        }
        index = (int)seqB[i];
        seq_char[index] = 1;
    }

    s.texto = calloc(tam_char + 1, sizeof(char));
    int j = 0;
    for (int i = 0; i < TAM_MAX_CHAR; i++)
    {
        if (seq_char[i] == 1)
        {
            s.texto[j] = i;
            j++;
        }
    }
    s.tam = tam_char;
    s.texto[s.tam + 1] = '\0';

    return s;
}

// sequencia b para calculo do caracteres unicos
t_sequencia calc_char_unicos(char *sequenciaA, char *sequenciaB)
{
    int *seq_char = calloc(TAM_MAX_CHAR, sizeof(int));

    return tam_char_seq(sequenciaA, sequenciaB, seq_char);
}

t_sequencia ler_entrada(char *filename)
{
    FILE *arquivo_sequencia = NULL;
    int i = 0;
    arquivo_sequencia = fopen(filename, "rt");
    t_sequencia s;

    if (arquivo_sequencia == NULL)
    {
        printf("Error reading file %s\n", filename);
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
    s.texto[i] = '\0';

    // fecha o arquivo
    fclose(arquivo_sequencia);
    return s;
}
// function prototypes
int get_index_of_character(char *str, char x, int len);
void print_matrix(int **x, int row, int col);
void calc_P_matrix_v2(int **P, char *b, int len_b, char *c, int len_c);
int lcs_yang_v2(int *DP, int *prev_row, int **P, char *A, char *B, char *C, int m, int n, int u);
int lcs(int *DP, int *prev_row, char *A, char *B, int m, int n);

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

void calc_P_matrix_v2(int **P, char *b, int len_b, char *c, int len_c)
{
#pragma omp parallel for
    for (int i = 0; i < len_c; i++)
    {
        for (int j = 1; j < len_b + 1; j++)
        {
            if (b[j - 1] == c[i])
            {
                P[i][j] = j;
            }
            else
            {
                P[i][j] = P[i][j - 1];
            }
        }
    }
}

int lcs_yang_v2(int *DP, int *prev_row, int **P, char *A, char *B, char *C, int m, int n, int u)
{
    for (int i = 1; i < m + 1; i++)
    {
        int c_i = get_index_of_character(C, A[i - 1], u);
        int t, s;

#pragma omp parallel for private(t, s) schedule(static)
        for (int j = 0; j < n + 1; j++)
        {
            t = (0 - P[c_i][j]) < 0;
            s = (0 - (prev_row[j] - (t * prev_row[P[c_i][j] - 1])));
            DP[j] = ((t ^ 1) || (s ^ 0)) * (prev_row[j]) + (!((t ^ 1) || (s ^ 0))) * (prev_row[P[c_i][j] - 1] + 1);
        }

#pragma omp parallel for schedule(static)
        for (int j = 0; j < n + 1; j++)
        {
            prev_row[j] = DP[j];
        }
    }
    printf("dpn%i \n", DP[n]);
    return DP[n];
}

int lcs(int *DP, int *prev_row, char *A, char *B, int m, int n)
{
    for (int i = 1; i < (m + 1); i++)
    {
        for (int j = 1; j < (n + 1); j++)
        {
            if (A[i - 1] == B[j - 1])
            {
                DP[j] = prev_row[j - 1] + 1;
            }
            else
            {
                DP[j] = max(prev_row[j], DP[j - 1]);
            }
        }

        for (int j = 0; j < n + 1; j++)
        {
            prev_row[j] = DP[j];
        }
    }

    return DP[n];
}

int main(int argc, char *argv[])
{

    t_sequencia sequenciaA;
    t_sequencia sequenciaB;
    t_sequencia sequenciaC;

    double start_time,
        stop_time, parcent_match;
    int len_a, len_b;

    sequenciaA = ler_entrada(argv[1]);
    string_A = sequenciaA.texto;
    len_a = sequenciaA.tam;
    sequenciaB = ler_entrada(argv[2]);
    string_B = sequenciaB.texto;
    len_b = sequenciaB.tam;

    sequenciaC = calc_char_unicos(sequenciaA.texto, sequenciaB.texto);
    c_len = sequenciaC.tam;

    printf("Length of sequence 1: %d bp\n", len_a);
    printf("Length of sequence 2: %d bp\n", len_b);

    unique_chars_C = sequenciaC.texto;
    printf("%s \n", string_A);
    printf("%s \n", string_B);
    printf("%s \n", unique_chars_C);
    // printf("\n##################################\n");
    printf("\n######## Results ########\n");
    // printf("##################################\n");
    // looking at the number of available threads
#pragma omp parallel
    {
#pragma omp single
        {
            printf("Number of threads used: %d\n", omp_get_num_threads());
        }
    }
    // allocate memory for DP Results
    DP_Results = (int *)malloc((len_b + 1) * sizeof(int));
    dp_prev_row = (int *)malloc((len_b + 1) * sizeof(int));

    // allocate memory for P_Matrix array
    P_Matrix = (int **)malloc(c_len * sizeof(int *));
    for (int k = 0; k < c_len; k++)
    {
        P_Matrix[k] = (int *)calloc((len_b + 1), sizeof(int));
    }

    start_time = omp_get_wtime();
    calc_P_matrix_v2(P_Matrix, string_B, len_b, unique_chars_C, c_len);

    int res = lcs_yang_v2(DP_Results, dp_prev_row, P_Matrix, string_A, string_B, unique_chars_C, len_a, len_b, c_len);
    stop_time = omp_get_wtime();
    printf("Length of the LCS is: %d\n", res);
    parcent_match = ((double)res / (double)len_a) * 100.0;
    printf("%.2f%% of the first sequence matches with second one\n", parcent_match);
    printf("Total time taken: %lf seconds\n", stop_time - start_time);

    free(P_Matrix);
    free(DP_Results);
    return 0;
}

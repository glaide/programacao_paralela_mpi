
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "omp.h"

// tamanho tabela ASCII para possíveis entradas
#define TAM_MAX_CHAR 256
#define NUM_THREADS 8
#define max(x, y) ((x) > (y) ? (x) : (y))

typedef struct sequencia
{
    char *texto;
    long tam;
} t_sequencia;

t_sequencia tam_char_seq(char *seqA, char *seqB, int *seq_char);
t_sequencia ler_entrada(char *filename);
int indice_char_unico(char *str, char x, int len);
void calc_matriz_P(int **P, char *b, long len_b, char *c, long len_c);
int calc_dif_distancia(int *atual, int *anterior, int **P, t_sequencia A, t_sequencia B, t_sequencia C);
int lcs(int *DP, int *prev_row, t_sequencia a, t_sequencia b);

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

int indice_char_unico(char *str, char x, int len)
{
    for (int i = 0; i < len; i++)
    {
        if (str[i] == x)
        {

            return i;
        }
    }
    return -1;
}

void calc_matriz_P(int **P, char *b, long len_b, char *c, long len_c)
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

int calc_dif_distancia(int *atual, int *anterior, int **P, t_sequencia A, t_sequencia B, t_sequencia C)
{
    for (int i = 1; i < A.tam + 1; i++)
    {
        int c_i = indice_char_unico(C.texto, A.texto[i - 1], C.tam);
        int t, s;

#pragma omp parallel for private(t, s) schedule(static)
        for (int j = 0; j < B.tam + 1; j++)
        {
            t = (0 - P[c_i][j]) < 0;
            s = (0 - (anterior[j] - (t * anterior[P[c_i][j] - 1])));
            atual[j] = ((t ^ 1) || (s ^ 0)) * (anterior[j]) + (!((t ^ 1) || (s ^ 0))) * (anterior[P[c_i][j] - 1] + 1);
        }

#pragma omp parallel for schedule(static)
        for (int j = 0; j < B.tam + 1; j++)
        {
            anterior[j] = atual[j];
        }
    }

    return atual[B.tam];
}

int lcs(int *DP, int *prev_row, t_sequencia a, t_sequencia b)
{
    for (int i = 1; i < (a.tam + 1); i++)
    {
        for (int j = 1; j < (b.tam + 1); j++)
        {
            if (a.texto[i - 1] == b.texto[j - 1])
            {
                DP[j] = prev_row[j - 1] + 1;
            }
            else
            {
                DP[j] = max(prev_row[j], DP[j - 1]);
            }
        }

        for (int j = 0; j < b.tam + 1; j++)
        {
            prev_row[j] = DP[j];
        }
    }

    return DP[b.tam];
}

int main(int argc, char *argv[])
{

    t_sequencia sequenciaA;
    t_sequencia sequenciaB;
    t_sequencia sequenciaC;
    int *linha_atual;
    int *linha_ant;

    int **P_Matrix;
    double start_time,
        stop_time;

    sequenciaA = ler_entrada(argv[1]);

    sequenciaB = ler_entrada(argv[2]);

    sequenciaC = calc_char_unicos(sequenciaA.texto, sequenciaB.texto);

    printf("Tamanho da sequencia A: %ld bp\n", sequenciaA.tam);
    printf("Tamanho da sequencia B: %ld bp\n", sequenciaB.tam);

    // printf("\n##################################\n");
    printf("\n######## Resultados ########\n");
    omp_set_num_threads(NUM_THREADS);

#pragma omp parallel
    {
#pragma omp single
        {
            printf("Número de threads: %d\n", omp_get_num_threads());
        }
    }

    linha_atual = (int *)malloc((sequenciaB.tam + 1) * sizeof(int));
    linha_ant = (int *)malloc((sequenciaB.tam + 1) * sizeof(int));

    // allocate memory for P_Matrix array
    P_Matrix = (int **)malloc(sequenciaC.tam * sizeof(int *));
    for (int k = 0; k < sequenciaC.tam; k++)
    {
        P_Matrix[k] = (int *)calloc((sequenciaB.tam + 1), sizeof(int));
    }

    start_time = omp_get_wtime();
    calc_matriz_P(P_Matrix, sequenciaB.texto, sequenciaB.tam, sequenciaC.texto, sequenciaC.tam);

    int res = calc_dif_distancia(linha_atual, linha_ant, P_Matrix, sequenciaA, sequenciaB, sequenciaC);

    stop_time = omp_get_wtime();
    printf("Tamanho do LCS é: %d\n", res);

    printf("tempo total: %lf seconds\n", stop_time - start_time);

    free(P_Matrix);
    free(linha_atual);
    return 0;
}

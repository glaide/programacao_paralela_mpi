#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main()
{
    char *texto = "abcdefghijklmnopqrstuvwxyz";
    for (int i = 0; i < strlen(texto); i++)
    {
        printf("%i ", (int)texto[i]);
    }
    printf("\n\n tamanho: %li\n", strlen(texto));
}
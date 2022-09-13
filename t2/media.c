#include <stdio.h>
#include <stdlib.h>
#define tam_max 5

double calcu_media(double valores[])
{
    double soma = 0.00000;
    for (int i = 0; i < tam_max; i++)
    {
        soma += valores[i];
    }
    return (soma / 5);
}
int main()
{
    double sequencial[tam_max];

    for (int i = 0; i < tam_max; i++)
    {
        scanf("%lf ", &sequencial[i]);
    }

    double media = calcu_media(sequencial);

    printf("Média: %lf \n", media);
}
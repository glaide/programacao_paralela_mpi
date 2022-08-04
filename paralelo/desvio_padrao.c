#include <stdio.h>
#include <stdlib.h>
#define tam_max 20

double calcu_media(double valores[])
{
    double soma = 0.00000;
    for (int i = 0; i < tam_max; i++)
    {
        soma += valores[i];
    }
    return (soma / 20);
}
int main()
{
    double sequencial[tam_max],
        th_2[tam_max], th_4[tam_max], media_1, media_2, media_3;

    for (int i = 0; i < tam_max; i++)
    {
        scanf("%lf ", &sequencial[i]);
    }

    media_1 = calcu_media(sequencial);

    printf("Média: %lf \n", media_1);
}
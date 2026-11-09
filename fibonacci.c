#include <stdio.h>
#include <stdlib.h>

int main()
{
    int n;
    unsigned long long fib[3];

    fib[0] = 0;
    fib[1] = 1;

    printf("Digite qual numero de Fibonacci deseja calcular: ");
    scanf("%d", &n);

    for (int i = 2; i <= n; i++)
        fib[i % 3] = fib[(i + 1) % 3] + fib[(i + 2) % 3];

    printf("%llu\n", fib[n % 3]);

    return 0;
}

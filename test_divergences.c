#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Teste da função maior_valor (versão com bug)
double maior_valor_buggy(double **array, int linha, int coluna, int i, int j){
    double max = fabs(array[i][j]);
    for(int i = 3; i < linha; i++){
        for(int j = 2; j < coluna; j++){
            if(array[i][j] > max)  // BUG: compara valor direto com absoluto
                max = fabs(array[i][j]);
        }
    }
    return max;
}

// Versão corrigida
double maior_valor_correto(double **array, int linha, int coluna, int i, int j){
    double max = fabs(array[i][j]);
    for(int i = 3; i < linha; i++){
        for(int j = 2; j < coluna; j++){
            if(fabs(array[i][j]) > max)  // CORRETO: compara absolutos
                max = fabs(array[i][j]);
        }
    }
    return max;
}

int main(){
    int linha = 6, coluna = 6;
    
    // Aloca matriz teste
    double **array = malloc(linha * sizeof(double*));
    for(int i = 0; i < linha; i++){
        array[i] = malloc(coluna * sizeof(double));
    }
    
    // Preenche com valores teste (incluindo negativos grandes)
    for(int i = 0; i < linha; i++){
        for(int j = 0; j < coluna; j++){
            array[i][j] = 1.0;
        }
    }
    
    // Coloca um valor negativo grande que deveria ser o máximo
    array[4][3] = -15.0;  // |15| > qualquer outro valor
    array[3][2] = 5.0;    // Valor inicial
    
    double resultado_buggy = maior_valor_buggy(array, linha, coluna, 3, 2);
    double resultado_correto = maior_valor_correto(array, linha, coluna, 3, 2);
    
    printf("Resultado com bug: %.2f\n", resultado_buggy);
    printf("Resultado correto: %.2f\n", resultado_correto);
    printf("Divergência: %s\n", (resultado_buggy != resultado_correto) ? "SIM" : "NÃO");
    
    // Libera memória
    for(int i = 0; i < linha; i++){
        free(array[i]);
    }
    free(array);
    
    return 0;
}
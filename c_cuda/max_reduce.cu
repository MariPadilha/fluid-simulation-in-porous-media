#include <stdio.h>
#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <math.h>

// Functor para calcular o máximo dos valores absolutos
struct MaxAbsOp {
    __device__ __forceinline__
    double operator()(const double &a, const double &b) const {
        return fmax(fabs(a), fabs(b));
    }
};

double max_reduce(double* dev_matriz_linearizada, int imax, int jmax){
    int tamanho = imax * jmax;
    double *dev_max;
    cudaMalloc(&dev_max, sizeof(double));
    void *d_temp_storage = NULL;
    size_t temp_storage_bytes = 0;

    // Inicializa com o valor absoluto do primeiro elemento
    double primeiro_valor;
    cudaMemcpy(&primeiro_valor, dev_matriz_linearizada, sizeof(double), cudaMemcpyDeviceToHost);
    double valor_inicial = fabs(primeiro_valor);
    cudaMemcpy(dev_max, &valor_inicial, sizeof(double), cudaMemcpyHostToDevice);

    // Usa operador customizado para máximo dos valores absolutos
    MaxAbsOp max_abs_op;
    cub::DeviceReduce::Reduce(d_temp_storage, temp_storage_bytes, dev_matriz_linearizada, dev_max, tamanho, max_abs_op, valor_inicial);
    cudaMalloc(&d_temp_storage, temp_storage_bytes);

    cub::DeviceReduce::Reduce(d_temp_storage, temp_storage_bytes, dev_matriz_linearizada, dev_max, tamanho, max_abs_op, valor_inicial);

    double h_max;
    cudaMemcpy(&h_max, dev_max, sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(dev_max);
    cudaFree(d_temp_storage);

    return h_max;
}

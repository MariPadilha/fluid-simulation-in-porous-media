#include <stdio.h>
#include <cuda_runtime.h>

// Kernel simples
__global__ void helloFromGPU() {
    printf("Hello from GPU! (thread %d, block %d)\n", threadIdx.x, blockIdx.x);
}

int main() {
    // 1. Verifica se há dispositivo CUDA disponível
    int deviceCount = 0;
    cudaError_t err = cudaGetDeviceCount(&deviceCount);
    if (err != cudaSuccess) {
        printf("Erro ao acessar CUDA: %s\n", cudaGetErrorString(err));
        return 1;
    }

    if (deviceCount == 0) {
        printf("Nenhuma GPU CUDA detectada.\n");
        return 1;
    }

    printf("GPUs CUDA detectadas: %d\n", deviceCount);

    // 2. Mostra informações da GPU
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    printf("Usando GPU: %s\n", prop.name);
    printf("Compute capability: %d.%d\n", prop.major, prop.minor);
    printf("Memória global total: %.2f GB\n", prop.totalGlobalMem / (1024.0 * 1024 * 1024));

    // 3. Executa o kernel simples
    helloFromGPU<<<2, 4>>>();
    cudaDeviceSynchronize();

    // 4. Verifica erros de execução
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("Erro durante execução do kernel: %s\n", cudaGetErrorString(err));
        return 1;
    }

    printf("Execução concluída com sucesso!\n");
    return 0;
}

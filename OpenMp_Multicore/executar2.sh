#!/bin/bash
#==============================================================================
# Script de Benchmark com perf e OpenMP
# Coleta tempos individuais, média e métricas de hardware por thread count
#==============================================================================

PROGRAM=./cylinder_solver.out
OUTPUT_DIR=perf_results
OUTPUT_CSV="$OUTPUT_DIR/resultados_perf.csv"
OUTPUT_RESUMO="$OUTPUT_DIR/resumo_perf.txt"
NUM_REPETICOES=1
THREADS_LIST=(32)

mkdir -p "$OUTPUT_DIR"
echo "threads,run1_time,run2_time,run3_time,mean_time,cycles,instructions,cache_misses,branch_misses" > "$OUTPUT_CSV"
echo "📊 Iniciando benchmarks com perf" | tee "$OUTPUT_RESUMO"

for threads in "${THREADS_LIST[@]}"; do
    echo -e "\n───────────────────────────────────────────────" | tee -a "$OUTPUT_RESUMO"
    echo "Executando com $threads thread(s) - $NUM_REPETICOES repetições" | tee -a "$OUTPUT_RESUMO"
    echo "───────────────────────────────────────────────" | tee -a "$OUTPUT_RESUMO"

    export OMP_NUM_THREADS=$threads
    tempos_reais=()

    #----------------------------------------------
    # Execuções individuais
    #----------------------------------------------
    for i in $(seq 1 $NUM_REPETICOES); do
        echo "  → Execução $i de $NUM_REPETICOES com $threads threads" | tee -a "$OUTPUT_RESUMO"

        PERF_RUN_OUTPUT="${OUTPUT_DIR}/perf_${threads}threads_run${i}.txt"
        SAIDA_RUN_OUTPUT="${OUTPUT_DIR}/saida_${threads}threads_run${i}.txt"

        perf stat -d \
            -e cycles,instructions,cache-misses,branch-misses \
            $PROGRAM > "$SAIDA_RUN_OUTPUT" 2> "$PERF_RUN_OUTPUT"

        tempo_real_i=$(grep "seconds time elapsed" "$PERF_RUN_OUTPUT" | awk '{print $1}' | tr ',' '.')
        tempos_reais+=($tempo_real_i)
        echo "    Tempo real (execução $i): ${tempo_real_i}s" | tee -a "$OUTPUT_RESUMO"
    done

    #----------------------------------------------
    # Calcular média
    #----------------------------------------------
    soma=0
    for t in "${tempos_reais[@]}"; do
        soma=$(echo "$soma + $t" | bc -l)
    done
    tempo_medio=$(echo "scale=6; $soma / ${#tempos_reais[@]}" | bc -l)
    echo "    → Tempo médio: ${tempo_medio}s" | tee -a "$OUTPUT_RESUMO"

    #----------------------------------------------
    # Executar perf agregado (para pegar contadores médios)
    #----------------------------------------------
    PERF_OUTPUT="${OUTPUT_DIR}/perf_${threads}threads.txt"
    perf stat -d -r $NUM_REPETICOES \
        -e cycles,instructions,cache-misses,branch-misses \
        $PROGRAM > /dev/null 2> "$PERF_OUTPUT"

    # Extrair métricas do perf médio
    cycles=$(grep "cycles" "$PERF_OUTPUT" | awk '{print $1}' | tr -d ',')
    instr=$(grep "instructions" "$PERF_OUTPUT" | awk '{print $1}' | tr -d ',')
    cachemiss=$(grep "cache-misses" "$PERF_OUTPUT" | awk '{print $1}' | tr -d ',')
    branchmiss=$(grep "branch-misses" "$PERF_OUTPUT" | awk '{print $1}' | tr -d ',')

    #----------------------------------------------
    # Salvar no CSV (tempos individuais + média + métricas)
    #----------------------------------------------
    echo "$threads,${tempos_reais[0]},${tempos_reais[1]},${tempos_reais[2]},$tempo_medio,$cycles,$instr,$cachemiss,$branchmiss" >> "$OUTPUT_CSV"
done

echo -e "\n✅ Resultados salvos em:"
echo "   - $OUTPUT_RESUMO"
echo "   - $OUTPUT_CSV"

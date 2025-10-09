#!/bin/bash

#==============================================================================
# Script de Benchmark - Análise de Performance com OpenMP
#==============================================================================
# Executa testes de performance variando número de threads
# - 10 execuções de aquecimento (descartadas)
# - 5 execuções oficiais (calcula média)
# - Salva resultados em CSV com todos os tempos (real, user, sys)
#==============================================================================

#------------------------------------------------------------------------------
# CONFIGURAÇÕES
#------------------------------------------------------------------------------
OUTPUT_CSV="resultados_benchmark.csv"
OUTPUT_DETALHADO="resultados_detalhados.csv"
OUTPUT_RESUMO="resumo_benchmark.txt"
PROGRAM="./cylinder_solver.out"
THREADS_ARRAY=(1 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32)
NUM_WARMUP=10
NUM_EXECUCOES=5

#------------------------------------------------------------------------------
# FUNÇÃO: Converter tempo do formato "2m25.339s" para segundos
#------------------------------------------------------------------------------
converter_tempo() {
    local tempo=$1
    
    # Remove espaços em branco
    tempo=$(echo "$tempo" | tr -d ' ')
    
    # Se contém 'm' (minutos)
    if [[ $tempo == *m* ]]; then
        local minutos=$(echo "$tempo" | cut -d'm' -f1)
        local segundos=$(echo "$tempo" | cut -d'm' -f2 | sed 's/s$//')
        echo "$minutos * 60 + $segundos" | bc -l
    else
        # Apenas segundos (remove 's' do final)
        echo "$tempo" | sed 's/s$//'
    fi
}

#------------------------------------------------------------------------------
# FUNÇÃO: Calcular média de um array
#------------------------------------------------------------------------------
calcular_media() {
    local -n arr=$1
    local soma=0
    local count=${#arr[@]}
    
    for valor in "${arr[@]}"; do
        soma=$(echo "$soma + $valor" | bc -l)
    done
    
    echo "scale=6; $soma / $count" | bc -l
}

#------------------------------------------------------------------------------
# FUNÇÃO: Calcular desvio padrão
#------------------------------------------------------------------------------
calcular_desvio_padrao() {
    local -n arr=$1
    local media=$2
    local soma_quadrados=0
    local count=${#arr[@]}
    
    for valor in "${arr[@]}"; do
        local diff=$(echo "$valor - $media" | bc -l)
        soma_quadrados=$(echo "$soma_quadrados + ($diff * $diff)" | bc -l)
    done
    
    echo "scale=6; sqrt($soma_quadrados / $count)" | bc -l
}

#------------------------------------------------------------------------------
# FUNÇÃO: Executar programa e capturar tempos
#------------------------------------------------------------------------------
executar_e_medir() {
    local temp_file=$(mktemp)
    
    # Executa e captura tempos
    { time $PROGRAM > /dev/null 2>&1; } 2> "$temp_file"
    
    # Extrai tempos
    local real_bruto=$(grep "real" "$temp_file" | awk '{print $2}')
    local user_bruto=$(grep "user" "$temp_file" | awk '{print $2}')
    local sys_bruto=$(grep "sys" "$temp_file" | awk '{print $2}')
    
    # Converte para segundos
    local real=$(converter_tempo "$real_bruto")
    local user=$(converter_tempo "$user_bruto")
    local sys=$(converter_tempo "$sys_bruto")
    
    rm -f "$temp_file"
    
    # Retorna valores separados por vírgula
    echo "$real,$user,$sys"
}

#------------------------------------------------------------------------------
# INICIALIZAÇÃO
#------------------------------------------------------------------------------
DATA_HORA=$(date '+%Y-%m-%d %H:%M:%S')
echo "=============================================================================="
echo "  BENCHMARK - Análise de Performance com OpenMP"
echo "=============================================================================="
echo "Data/Hora: $DATA_HORA"
echo "Programa: $PROGRAM"
echo "Aquecimento: $NUM_WARMUP execuções"
echo "Medições: $NUM_EXECUCOES execuções por configuração"
echo "Threads testadas: ${THREADS_ARRAY[@]}"
echo "=============================================================================="
echo ""

# Verificar se programa existe
if [ ! -f "$PROGRAM" ]; then
    echo "ERRO: Programa '$PROGRAM' não encontrado!"
    exit 1
fi

#------------------------------------------------------------------------------
# CRIAR ARQUIVOS CSV
#------------------------------------------------------------------------------

# CSV com médias (formato simples)
{
    echo "# Benchmark gerado em: $DATA_HORA"
    echo "# Programa: $PROGRAM"
    echo "# Aquecimento: $NUM_WARMUP execuções | Medições: $NUM_EXECUCOES execuções"
    echo "#"
    echo "num_threads,real_avg,user_avg,sys_avg,real_stddev,speedup_real,eficiencia"
} > "$OUTPUT_CSV"

# CSV detalhado com todas as execuções
{
    echo "# Benchmark detalhado - gerado em: $DATA_HORA"
    echo "# Todas as execuções individuais"
    echo "#"
    echo -n "num_threads,real_avg,user_avg,sys_avg,real_stddev,user_stddev,sys_stddev"
    for i in $(seq 1 $NUM_EXECUCOES); do
        echo -n ",real_exec${i},user_exec${i},sys_exec${i}"
    done
    echo ""
} > "$OUTPUT_DETALHADO"

# Arquivo de resumo em texto
{
    echo "=============================================================================="
    echo "  RELATÓRIO DE BENCHMARK - Análise de Performance OpenMP"
    echo "=============================================================================="
    echo "Data/Hora: $DATA_HORA"
    echo "Programa: $PROGRAM"
    echo "Configuração:"
    echo "  - Execuções de aquecimento: $NUM_WARMUP"
    echo "  - Execuções medidas: $NUM_EXECUCOES"
    echo "  - Threads testadas: ${THREADS_ARRAY[@]}"
    echo "=============================================================================="
    echo ""
} > "$OUTPUT_RESUMO"

#------------------------------------------------------------------------------
# VARIÁVEL PARA ARMAZENAR TEMPO BASE (1 thread)
#------------------------------------------------------------------------------
TEMPO_BASE=0

#------------------------------------------------------------------------------
# LOOP PRINCIPAL - Testa cada configuração de threads
#------------------------------------------------------------------------------
for num_threads in "${THREADS_ARRAY[@]}"; do
    export OMP_NUM_THREADS=$num_threads
    
    echo "" | tee -a "$OUTPUT_RESUMO"
    echo "┌─────────────────────────────────────────────────────────────────────────┐" | tee -a "$OUTPUT_RESUMO"
    echo "│ Threads: $num_threads" | tee -a "$OUTPUT_RESUMO"
    echo "└─────────────────────────────────────────────────────────────────────────┘" | tee -a "$OUTPUT_RESUMO"
    
    #--------------------------------------------------------------------------
    # FASE 1: Aquecimento
    #--------------------------------------------------------------------------
    echo "  [1/2] Aquecimento ($NUM_WARMUP execuções)..."
    for warmup in $(seq 1 $NUM_WARMUP); do
        $PROGRAM > /dev/null 2>&1
        printf "    Progresso: [%2d/%2d]\r" $warmup $NUM_WARMUP
    done
    echo "    Progresso: [$NUM_WARMUP/$NUM_WARMUP] ✓"
    
    #--------------------------------------------------------------------------
    # FASE 2: Medições Oficiais
    #--------------------------------------------------------------------------
    echo "  [2/2] Medições oficiais ($NUM_EXECUCOES execuções)..." | tee -a "$OUTPUT_RESUMO"
    
    declare -a tempos_real
    declare -a tempos_user
    declare -a tempos_sys
    
    for exec_num in $(seq 1 $NUM_EXECUCOES); do
        resultado=$(executar_e_medir)
        
        real=$(echo "$resultado" | cut -d',' -f1)
        user=$(echo "$resultado" | cut -d',' -f2)
        sys=$(echo "$resultado" | cut -d',' -f3)
        
        tempos_real+=($real)
        tempos_user+=($user)
        tempos_sys+=($sys)
        
        printf "    Execução %d: real=%8.3fs | user=%8.3fs | sys=%6.3fs\n" \
               $exec_num $real $user $sys | tee -a "$OUTPUT_RESUMO"
    done
    
    #--------------------------------------------------------------------------
    # FASE 3: Calcular Estatísticas
    #--------------------------------------------------------------------------
    media_real=$(calcular_media tempos_real)
    media_user=$(calcular_media tempos_user)
    media_sys=$(calcular_media tempos_sys)
    
    stddev_real=$(calcular_desvio_padrao tempos_real $media_real)
    stddev_user=$(calcular_desvio_padrao tempos_user $media_user)
    stddev_sys=$(calcular_desvio_padrao tempos_sys $media_sys)
    
    # Guardar tempo base (1 thread)
    if [ $num_threads -eq 1 ]; then
        TEMPO_BASE=$media_real
    fi
    
    # Calcular speedup e eficiência
    if (( $(echo "$media_real > 0" | bc -l) )); then
        speedup=$(echo "scale=6; $TEMPO_BASE / $media_real" | bc -l)
        eficiencia=$(echo "scale=6; ($speedup / $num_threads) * 100" | bc -l)
    else
        speedup=0
        eficiencia=0
    fi
    
    echo "" | tee -a "$OUTPUT_RESUMO"
    echo "  ┌─────────────────────────────────────────────────────────────────────┐" | tee -a "$OUTPUT_RESUMO"
    printf "  │ MÉDIAS:       real=%8.3fs | user=%8.3fs | sys=%6.3fs     │\n" \
           $media_real $media_user $media_sys | tee -a "$OUTPUT_RESUMO"
    printf "  │ DESVIO PADRÃO: real=%8.3fs | user=%8.3fs | sys=%6.3fs     │\n" \
           $stddev_real $stddev_user $stddev_sys | tee -a "$OUTPUT_RESUMO"
    printf "  │ SPEEDUP: %.2fx | EFICIÊNCIA: %.2f%%                                │\n" \
           $speedup $eficiencia | tee -a "$OUTPUT_RESUMO"
    echo "  └─────────────────────────────────────────────────────────────────────┘" | tee -a "$OUTPUT_RESUMO"
    
    #--------------------------------------------------------------------------
    # FASE 4: Salvar nos CSVs
    #--------------------------------------------------------------------------
    
    # CSV simples (médias)
    printf "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.2f\n" \
           $num_threads $media_real $media_user $media_sys $stddev_real $speedup $eficiencia \
           >> "$OUTPUT_CSV"
    
    # CSV detalhado (todas as execuções)
    {
        printf "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f" \
               $num_threads $media_real $media_user $media_sys $stddev_real $stddev_user $stddev_sys
        for i in $(seq 0 $((NUM_EXECUCOES-1))); do
            printf ",%.6f,%.6f,%.6f" "${tempos_real[$i]}" "${tempos_user[$i]}" "${tempos_sys[$i]}"
        done
        echo ""
    } >> "$OUTPUT_DETALHADO"
    
    # Limpar arrays para próxima iteração
    unset tempos_real tempos_user tempos_sys
done

#------------------------------------------------------------------------------
# FINALIZAÇÃO
#------------------------------------------------------------------------------
echo "" | tee -a "$OUTPUT_RESUMO"
echo "==============================================================================" | tee -a "$OUTPUT_RESUMO"
echo "  BENCHMARK CONCLUÍDO!" | tee -a "$OUTPUT_RESUMO"
echo "==============================================================================" | tee -a "$OUTPUT_RESUMO"
echo "Arquivos gerados:" | tee -a "$OUTPUT_RESUMO"
echo "  1. $OUTPUT_CSV (médias e speedup)" | tee -a "$OUTPUT_RESUMO"
echo "  2. $OUTPUT_DETALHADO (todas as execuções)" | tee -a "$OUTPUT_RESUMO"
echo "  3. $OUTPUT_RESUMO (relatório completo)" | tee -a "$OUTPUT_RESUMO"
echo "" | tee -a "$OUTPUT_RESUMO"
echo "Resumo de Speedup:" | tee -a "$OUTPUT_RESUMO"
awk -F',' 'NR>4 {printf "  %2d threads: tempo=%8.3fs | speedup=%5.2fx | eficiência=%6.2f%%\n", $1, $2, $6, $7}' "$OUTPUT_CSV" | tee -a "$OUTPUT_RESUMO"
echo "==============================================================================" | tee -a "$OUTPUT_RESUMO"
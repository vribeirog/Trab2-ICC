#!/bin/bash

METRICAS="L3 L2CACHE FLOPS_DP FLOPS_AVX" # seus grupos LIKWID
OUTDIR="results"          # pasta de destino das saídas
CPU=3

mkdir -p "$OUTDIR"

NSIZES="32 64 128 256 512 1000" # 2000 4000 8000 9000 10000 20000"

V1="v1"
V2="v2"

for n in $NSIZES
do
    for m in ${METRICAS}
    do

        LIKWID_OUT="${OUTDIR}/${V1}/${m}_${n}.txt"
        echo "metrica=$m para N=$n" > /dev/tty

    # Executa o programa com likwid-perfctr + inputs via heredoc
        likwid-perfctr -C ${CPU} -g ${m} -o ${LIKWID_OUT} -m ./v1/cgSolver \
<<EOF
$n
EOF
        done
done

# Geração dos .csv
echo "Gerando CSVs"

LIKWID_AVX_CSV="${OUTDIR}/${V1}/FLOPS_AVX.csv"
echo "N,gradiente_FLOPS_AVX,residuo_FLOPS_AVX" > ${LIKWID_AVX_CSV}

for m in ${METRICAS}
do
    #Criando CSV
    LIKWID_CSV="${OUTDIR}/${V1}/${m}.csv"

    # string a ser procurada no arquivo .txt
    METRIC_TO_GREP=""
    if   [ "$m" = "L3" ]; then
        METRIC_TO_GREP="L3 bandwidth"
    elif [ "$m" = "L2CACHE" ]; then
        METRIC_TO_GREP="L2 miss ratio"
    elif [ "$m" = "FLOPS_DP" ]; then
        METRIC_TO_GREP="DP MFLOP/s"
    fi

echo -e "\nProcessando $m. \nMétrica: $METRIC_TO_GREP. \nSaída: $LIKWID_CSV\n" >/dev/tty

    # Escreve o cabeçalho no CSV final
    echo "N,gradiente_${m},residuo_${m}" > ${LIKWID_CSV}

    # Realiza a procura em cada um dos testes realizados
    for n in $NSIZES
    do
        LIKWID_OUT="${OUTDIR}/${V1}/${m}_${n}.txt"

        if [ ! -f "$LIKWID_OUT" ]; then
            echo "Arquivo não encontrado: $LIKWID_OUT" >/dev/tty
            continue
        fi

        # Extrai a métrica para GRADIENTE_CONJUGADO (somando todas as iterações)
        METRICS=$(grep "$METRIC_TO_GREP" "$LIKWID_OUT")


        if [ "$m" = "FLOPS_DP" ]; then
            METRIC_GRADIENTE=$(grep "$METRIC_TO_GREP" "$LIKWID_OUT" | grep -v "AVX" | sed -n '1p' | grep -oE '[0-9]+\.[0-9]+')
            METRIC_RESIDUO=$(grep "$METRIC_TO_GREP" "$LIKWID_OUT"   | grep -v "AVX" | sed -n '2p' | grep -oE '[0-9]+\.[0-9]+')

            METRIC_GRADIENTE_AVX=$(grep "AVX DP MFLOP/s" "$LIKWID_OUT" | sed -n '1p' | grep -oE '[0-9]+\.[0-9]+')
            METRIC_RESIDUO_AVX=$(grep "AVX DP MFLOP/s" "$LIKWID_OUT"   | sed -n '2p' | grep -oE '[0-9]+\.[0-9]+')

            echo "${n},${METRIC_GRADIENTE_AVX},${METRIC_RESIDUO_AVX}" >> ${LIKWID_AVX_CSV}
        else
            METRIC_GRADIENTE=$(grep "$METRIC_TO_GREP" "$LIKWID_OUT" | sed -n '1p' | grep -oE '[0-9]+\.[0-9]+')
	    METRIC_RESIDUO=$(grep "$METRIC_TO_GREP" "$LIKWID_OUT"   | sed -n '2p' | grep -oE '[0-9]+\.[0-9]+')
          
        fi

        if [ -z "$METRIC_GRADIENTE" ]; then METRIC_GRADIENTE="0"; fi
        if [ -z "$METRIC_RESIDUO"   ]; then METRIC_RESIDUO="0"; fi
        if [ -z "$METRIC_GRADIENTE_AVX" ]; then METRIC_GRADIENTE_AVX="0"; fi
        if [ -z "$METRIC_RESIDUO_AVX"   ]; then METRIC_RESIDUO_AVX="0"; fi

        echo "${n},${METRIC_GRADIENTE},${METRIC_RESIDUO}" >> ${LIKWID_CSV}
    done
done

TIME_CSV="${OUTDIR}/${V1}/Exec_time.csv"

echo "N,Tempo_Gradiente_ms,Tempo_Residuo_ms" > ${TIME_CSV}

for n in $NSIZES
do
    SAIDA=$(./v1/cgSolver <<EOF
$n
EOF
)
    
    TIME_GRADIENTE=$(echo "$SAIDA"  | tail -n 2 | sed -n '1p')
    TIME_RESIDUO=$(echo "$SAIDA" | tail -n 2 | sed -n '2p')
    
    echo "${n},${TIME_GRADIENTE},${TIME_RESIDUO}" >> ${TIME_CSV}
done


mkdir -p "${OUTDIR}/${V1}/likwid_results/"
mv "${OUTDIR}/${V1}/"*.txt "${OUTDIR}/${V1}/likwid_results/"

# ================================

for n in $NSIZES
do
    for m in ${METRICAS}
    do

        LIKWID_OUT="${OUTDIR}/${V2}/${m}_${n}.txt"
        echo "metrica=$m para N=$n" > /dev/tty

    # Executa o programa com likwid-perfctr + inputs via heredoc
        likwid-perfctr -C ${CPU} -g ${m} -o ${LIKWID_OUT} -m ./v2/cgSolver \
<<EOF
$n
EOF
        done
done

# Geração dos .csv
echo "Gerando CSVs"

LIKWID_AVX_CSV="${OUTDIR}/${V2}/FLOPS_AVX.csv"
echo "N,gradiente_FLOPS_AVX,residuo_FLOPS_AVX" > ${LIKWID_AVX_CSV}

for m in ${METRICAS}
do
    #Criando CSV
    LIKWID_CSV="${OUTDIR}/${V2}/${m}.csv"

    # string a ser procurada no arquivo .txt
    METRIC_TO_GREP=""
    if   [ "$m" = "L3" ]; then
        METRIC_TO_GREP="L3 bandwidth"
    elif [ "$m" = "L2CACHE" ]; then
        METRIC_TO_GREP="L2 miss ratio"
    elif [ "$m" = "FLOPS_DP" ]; then
        METRIC_TO_GREP="DP MFLOP/s"
    fi

echo -e "\nProcessando $m. \nMétrica: $METRIC_TO_GREP. \nSaída: $LIKWID_CSV\n" >/dev/tty

    # Escreve o cabeçalho no CSV final
    echo "N,gradiente_${m},residuo_${m}" > ${LIKWID_CSV}

    # Realiza a procura em cada um dos testes realizados
    for n in $NSIZES
    do
        LIKWID_OUT="${OUTDIR}/${V2}/${m}_${n}.txt"

        if [ ! -f "$LIKWID_OUT" ]; then
            echo "Arquivo não encontrado: $LIKWID_OUT" >/dev/tty
            continue
        fi

        # Extrai a métrica para GRADIENTE_CONJUGADO (somando todas as iterações)
        METRICS=$(grep "$METRIC_TO_GREP" "$LIKWID_OUT")


        if [ "$m" = "FLOPS_DP" ]; then
            METRIC_GRADIENTE=$(grep "$METRIC_TO_GREP" "$LIKWID_OUT" | grep -v "AVX" | sed -n '1p' |  grep -oE '[0-9]+\.[0-9]+')
            METRIC_RESIDUO=$(grep "$METRIC_TO_GREP" "$LIKWID_OUT"   | grep -v "AVX" | sed -n '2p' |  grep -oE '[0-9]+\.[0-9]+')

            METRIC_GRADIENTE_AVX=$(grep "AVX DP MFLOP/s" "$LIKWID_OUT" | sed -n '1p' |  grep -oE '[0-9]+\.[0-9]+')
            METRIC_RESIDUO_AVX=$(grep "AVX DP MFLOP/s" "$LIKWID_OUT"   | sed -n '2p' |  grep -oE '[0-9]+\.[0-9]+')

            echo "${n},${METRIC_GRADIENTE_AVX},${METRIC_RESIDUO_AVX}" >> ${LIKWID_AVX_CSV}
        else
	    METRIC_GRADIENTE=$(grep "$METRIC_TO_GREP" "$LIKWID_OUT" | sed -n '1p' | grep -oE '[0-9]+\.[0-9]+')
            METRIC_RESIDUO=$(grep "$METRIC_TO_GREP" "$LIKWID_OUT"   | sed -n '2p' | grep -oE '[0-9]+\.[0-9]+')
          
        fi

        if [ -z "$METRIC_GRADIENTE" ]; then METRIC_GRADIENTE="0"; fi
        if [ -z "$METRIC_RESIDUO"   ]; then METRIC_RESIDUO="0"; fi
        if [ -z "$METRIC_GRADIENTE_AVX" ]; then METRIC_GRADIENTE_AVX="0"; fi
        if [ -z "$METRIC_RESIDUO_AVX"   ]; then METRIC_RESIDUO_AVX="0"; fi

        echo "${n},${METRIC_GRADIENTE},${METRIC_RESIDUO}" >> ${LIKWID_CSV}
    done
done

TIME_CSV="${OUTDIR}/${V2}/Exec_time.csv"

echo "N,Tempo_Gradiente_ms,Tempo_Residuo_ms" > ${TIME_CSV}

for n in $NSIZES
do
    SAIDA=$(./v2/cgSolver <<EOF
$n
EOF
)
    
    TIME_GRADIENTE=$(echo "$SAIDA"  | tail -n 2 | sed -n '1p')
    TIME_RESIDUO=$(echo "$SAIDA"  | tail -n 2 | sed -n '2p')
    
    echo "${n},${TIME_GRADIENTE},${TIME_RESIDUO}" >> ${TIME_CSV}
done


mkdir -p "${OUTDIR}/${V2}/likwid_results/"
mv "${OUTDIR}/${V2}/"*.txt "${OUTDIR}/${V2}/likwid_results/"

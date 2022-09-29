#!/bin/bash

gfa=$1
fa=$2
WD=$3
T=$4
BIN=$5

run () {
    i=$1
    tmp_fa="$WD/$i.read.fa"
    out_fa="$WD/$i.consensus.msa.fa"
    out_gfa="$WD/$i.gfa"
    tlog="$WD/$i.time"
    log="$WD/$i.log"
    log_gfa="$WD/$i.gfa.log"
    head -n $((2*i)) $fa | tail -n 2 > $tmp_fa
    /usr/bin/time -vo $tlog $BIN -r 2 -m 0 -s -i $gfa $tmp_fa > $out_fa 2> $log
    # $BIN -r 3 -m 0 -s -i $gfa $tmp_fa > $out_gfa 2> $log_gfa
}

mkdir -p $WD
N=$(grep -c "^>" $fa)

for i in $(seq 0 $((N/T - 1)))
do
    for t in $(seq 0 $((T-1)))
    do
	run $((i*T + $t + 1)) &
    done
    wait
done

REM=$((N % T))
if [[ ! $REM -eq "0" ]]
then
    i=$((i+1))
    for t in $(seq 0 $((REM-1)))
    do
	run $((i*T + $t + 1)) &
    done
    wait
fi

grep "abpoa_main" $WD/*.log | cut -f2- -d':'

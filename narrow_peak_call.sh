#!/bin/bash
cd /mnt/176/Changhai_ATAC/results/bwa/mergedLibrary/macs/narrow_peaks

atac_bam_files=`ls /mnt/176/Changhai_ATAC/results/bwa/mergedLibrary/*.bam`

for i in $atac_bam_files
do
    NAME=`echo -e  "$i\n" | cut -d . -f 1 | cut -d / -f 5`
    macs2 callpeak -t $i -g hs -n $NAME -B --nolambda --nomodel --call-summits \
    --outdir '/mnt/176/Changhai_ATAC/results/bwa/mergedLibrary/macs/narrow_peaks' --keep-dup all \
    --shift -75 --extsize 150 -p 0.01

done

echo 'done'

# -f BAMPE
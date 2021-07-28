#!/bin/bash

#path of sample contig
db="/ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/user/wangdaxi/pipeline_host_identitication/cytob.fasta"
sample=$1;
output=$2;

rm -f $output.psl && /ldfssz1/ST_INFECTION/P17Z10200N0536_Echinococcus/USER/wangdaxi/Programs/miniconda3/bin/blat -noHead $db $sample $output.psl
cat $output.psl|perl -ne '@a=split /\s+/;next if ($a[0]<500);$r=$a[0]/($a[0]+$a[1]);print "$a[9]\t$a[13]\t$r\t$a[10]\t$a[14]\n";'|sort -k3nr >$output.sorted.list
cat $output.sorted.list|head -n1|cut -f 2 >$output.sorted.shortlist

cat $db|grep -f $output.sorted.shortlist|perl -ne 'if(/\>\S+\s+(\S+\s+\S+)/){print "$1\n";}' >$output.matched.taxa

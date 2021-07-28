month=`date|awk '{print $2}'`
date=`date|awk '{print $3}'`
year=`date|awk '{print $6}'`
date=${year}${month}${date}

basedir=$1
samplesbatchsdate=$2


statfile=$basedir/$samplesbatchsdate/all.sum.assemble

#组装后的contigs总长度，总contig数，contig长度的中位值,insertsize平均值
AssemblyContigs=$basedir/$samplesbatchsdate/02.assemble/0203.megahit/final.contigs.fa
seqkitAllContigs=$basedir/$samplesbatchsdate/03.QC_asssemble_stats/megahit_seqkit_stats.txt
samtoolsAllContigs=$basedir/$samplesbatchsdate/03.QC_asssemble_stats/allcontigs_samtools_stats.txt
#1.组装后的contigs总长度
AllContigslenth=`awk '$1~!/^f/{print $5}' $seqkitAllContigs|sed '{s/,//g}'`
#2.总的contigs数
AllContigsNum=`grep -c ">" $AssemblyContigs`
#3.contig长度中位值
seqkitfx2tabAllContigs=$basedir/$samplesbatchsdate/03.QC_asssemble_stats/megahit_seqkit_fx2tab.txt
AllContigzhong=$basedir/$samplesbatchsdate/03.QC_asssemble_stats/megahit_median.contig.length
script=/ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/scripts/get_median.size.R
Rscript $script $seqkitfx2tabAllContigs $AllContigzhong
AllContigslengthzhong=`awk '{print $2}' $AllContigzhong`

#4.聚类后的条数
#4.insertsize
#AllContigsinsertsizq=`grep 'insert size average' $samtoolsAllContigs|awk '{print $5}'`

#注释得到的RNA病毒数目，比对上RdRp的数目，两者取交集的数目
#RNA_contigs=$basedir/$samplesbatchsdate/07.megahit/0701.diamond/${samplesbatchsdate}.CatRNAVirusContigs.id
#RdRp_contigs=$basedir/$samplesbatchsdate/07.megahit/0701.diamond/${samplesbatchsdate}.CatContigsProtein.blastp2rdrp.cid
#RNA_RdRp_contigs=$basedir/$samplesbatchsdate/07.megahit/0701.diamond/${samplesbatchsdate}.CatRNARdRpContigs.id
#RNA_contigs_sum=`wc -l $RNA_contigs`
#RdRp_contigs_sum=`wc -l $RdRp_contigs`
#RNA_RdRp_contigs_sum=`wc -l $RNA_RdRp_contigs`

#比对上多少个RdRp科
#RdRp_family_sum=`ls $basedir/$samplesbatchsdate/07.megahit/0702.rdrp|grep seqID|wc -l`

#checkv评估高完整性的contigs数目

#virfinder数目，virfinder和CAT取并集的数目，去完假阳后

echo -e "$AllContigslenth\t$AllContigsNum\t$AllContigslengthzhong" >> $statfile

echo end $samplesbatchsdate at:`date`

zdir=$1
sample=$2


#$sortmerna_fq1=$zdir/$sample/01.QC/0104.sortmeRNA/sortmerna_rmrRNA_fwd.fq.gz
#$s$ortmerna_fq2=$zdir/$sample/01.QC/0104.sortmeRNA/sortmerna_rmrRNA_rev.fq.gz

cox1_dir=$zdir/$sample/04.assemble_megahit_nr/0408.cox1_findhost
p1=$zdir/$sample/04.assemble_megahit_nr/0408.cox1_findhost/${sample}.centrifuge.KSout
if [ -f $p1 ];
then
    bash /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/scripts/cox1_searchbylineages.sh -t Chiroptera -s ${sample}.centrifuge.KSout -l S -i $cox1_dir -o $cox1_dir
    bash /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/scripts/cox1_searchbylineages.sh -t Rodentia -s ${sample}.centrifuge.KSout -l S -i $cox1_dir -o $cox1_dir
    bash /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/scripts/cox1_summary.sh -s ${sample}.centrifuge.KSout -l P -i $cox1_dir -o $cox1_dir
    echo $sample
#得到老鼠
fi
#mkdir -p $zdir/$sample/04.assemble_megahit_nr
#mkdir -p $zdir/$sample/04.assemble_megahit_nr/0408.cox1_findhost
#sh /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/scripts/cox1_centrifuge.sh $zdir $sample $cox1_dir $sortmerna_fq1 $sortmerna_fq2
#得到每个级别最高的
#sh /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/scripts/cox1_top_abundant_taxa.sh $cox1_dir {params.luname} {params.cox1}
#得到蝙蝠


#在脊椎动物
#bash /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/scripts/cox1_searchbylineages.sh -t Vertebrata -s ${sample}.centrifuge.KSout -l S -i $cox1_dir -o $cox1_dir
#得到丰度比较高的，在种级别
#bash /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/scripts/cox1_summary.sh -s ${sample}.centrifuge.KSout -l S -i $cox1_dir -o $cox1_dir
#在门级别


#


poutdir=$1 #basedir,也就是在哪个大目录下分析的这些样本
sample=$2 #样本分析目录的basename


blastp2rdrp_dir=$poutdir/$sample/04.assemble_megahit_nr/0402.CATpreProtein_blastp_rdrp
blastp2RdRp_out=$blastp2rdrp_dir/${sample}.CAT_prodigal_proteins_blastp2RdRp
#blastp2RdRp_out_besthit=$blastp2rdrp_dir/${sample}.CAT_prodigal_proteins.blastp2RdRp.filt
#rm -f $blastp2RdRp_out_besthit
echo 1
echo $blastp2RdRp_out
head $blastp2RdRp_out
#cat $blastp2RdRp_out|awk '{FS="\t" ; OFS="\t" ; if($11/$4 >= 0.5&&$12>50){print $0}}' > $blastp2RdRp_out_besthit
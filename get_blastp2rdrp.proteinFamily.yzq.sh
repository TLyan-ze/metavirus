poutdir=$1 #basedir,也就是在哪个大目录下分析的这些样本
sample=$2 #样本分析目录的basename
seqkit=/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/bin/seqkit
Rscript=/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/envs.multi-user/R/4.0.3/bin/Rscript
export R_LIBS="/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/envs.multi-user/R/4.0.3/lib/R/library:$R_LIBS"
get_rdrpFamily_proteinSets=/ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yanzeqin/metavirus/Lastest/get_rdrpFamily-proteinSets.R

CAT_dir=$poutdir/$sample/04.assemble_megahit_nr/0401.CAT_note_RNA
protein_fa=$CAT_dir/${sample}.predicted_proteins.faa
blastp2rdrp_dir=$poutdir/$sample/04.assemble_megahit_nr/0402.CATpreProtein_blastp_rdrp
blastp2RdRp_out=$blastp2rdrp_dir/${sample}.CAT_prodigal_proteins_blastp2RdRp
#筛选比对上rdrp的一致性大于50，覆盖该rdrp超过一半
#blastp2RdRp_out_besthit=$blastp2rdrp_dir/${sample}.CAT_prodigal_proteins.blastp2RdRp.filt
#rm -f $blastp2RdRp_out_besthit
#cat $blastp2RdRp_out|awk '{FS="\t" ; OFS="\t" ; if($11/$4 >= 0.5&&$12>50){print $0}}' > $blastp2RdRp_out_besthit

##获取contigId-proteinId-Family信息,为contig唯一个orf的bitscore
proteinId2Family=$blastp2rdrp_dir/${sample}.proteinId2Family_unfilt.list
rm -f $proteinId2Family
$Rscript $get_rdrpFamily_proteinSets $blastp2RdRp_out $proteinId2Family

##获取每个Family的蛋白序列
family_list=$blastp2rdrp_dir/${sample}.family.list
sed '1d' $proteinId2Family|awk '{print $3}'|sort|uniq > $family_list
#比对到的参考序列科目录
familyProtein_Fa_dir=$blastp2rdrp_dir/proteinFamily_fasta
rm -rf $familyProtein_Fa_dir
mkdir -p $familyProtein_Fa_dir

cat $family_list|while read family
do
	proteinFamily_id=$familyProtein_Fa_dir/${sample}.${family}.proteinId.list
	proteinFamily_fa=$familyProtein_Fa_dir/${sample}.${family}.proteins.fasta
	if [ $family != 'unassigned' ];then 
		grep -w $family $proteinId2Family|awk '{print $2}'|sort|uniq > $proteinFamily_id
		$seqkit grep -f $proteinFamily_id $protein_fa > $proteinFamily_fa
	fi
done

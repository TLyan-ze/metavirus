zsdir=$1
batch1datadir=$2
run_shell=$3
dthing=$4

awk '{print $2}' $batch1datadir|while read line
do
	outdir=$zsdir/shell
	sample=$zsdir/$line/02.assemble/megahit/${line}.contigs.fa
	output=$zsdir/$line/02.assemble/megahit/${line}
	if [ -f $outdir ];then
    	echo -e "$outdir\texit"
	else
    	mkdir -p $outdir
	fi
	runshell=$outdir/run_${dthing}_${line}.sh
	echo "sh $run_shell $sample $output" > $runshell
	cd $outdir
	qsub -clear -cwd -q st.q -P P20Z10200N0206 -binding linear:2 -l vf=1G,num_proc=1 $runshell
	#echo $runshell
done

zsdir=$1
batch1datadir=$2
run_shell=$3
dthing=$4

awk '{print $1}' $batch1datadir|while read samplebatchdate
do
	outdir=$zsdir/$samplebatchdate/shell
	if [ -f $outdir ];then
    	echo -e "$outdir\texit"
	else
    	mkdir -p $outdir
	fi

	runshell=$outdir/run_${dthing}_${samplebatchdate}.sh
	echo "sh $run_shell $zsdir $samplebatchdate" > $runshell
	cd $outdir
	qsub -clear -cwd -q st_short.q -P P20Z10200N0206 -binding linear:2 -l vf=1G,num_proc=1 $runshell
	#echo $runshell
done

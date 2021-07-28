snakemake=$1
cluster=$2
config=$3
outdir=$4

month=`date|awk '{print $2}'`
date=`date|awk '{print $3}'`
year=`date|awk '{print $6}'`
data=${year}${month}${date}
hour=`date|awk '{print $4}'|awk -F ':' '{print $1}'`


if [ ! -f $outdir/$data ];then
	mkdir $outdir/$data
fi
#集群运行目录
mkdir $outdir/$data/logs

out_nohup=$outdir/$data/snakemake_${hour}_nohup.txt

touch $out_nohup
chmod 775 $out_nohup

/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/bin/snakemake -s $snakemake -n --configfile $config --rerun-incomplete

echo start at:`date`

cp -f $snakemake $outdir/$data
cp -f $cluster $outdir/$data
cp -f $config $outdir/$data

 
cd $outdir/$data 
/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/bin/snakemake -s $snakemake --unlock --configfile $config --rerun-incomplete
nohup /hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/bin/snakemake \
--cluster "qsub -S /bin/bash -cwd -q {cluster.queue} -P {cluster.project} -l vf={cluster.mem},p={cluster.cores} -binding linear:{cluster.cores} -o {cluster.output} -e {cluster.error}" \
-s $snakemake \
--jobs 1600 \
--cluster-config $cluster \
--configfile $config \
--rerun-incomplete \
> $out_nohup 2>&1 &

echo end at:`date`

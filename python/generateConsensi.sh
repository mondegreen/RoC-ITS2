task(){
  BASE=`basename $1 -subreads.fq`
  DIR=`dirname $1`

  if test -f ${DIR}/${BASE}.done; then
    echo "skipping ${DIR}/${BASE}"
  else
    CONSENT-correct --in ${DIR}/${BASE}-subreads.fq \
      --out ${DIR}/${BASE}-corrected.fa --type ONT

    prank -F -o=${DIR}/${BASE}-corrected \
      -d=${DIR}/${BASE}-corrected.fa > /dev/null

    python3 smush-msa.py --msa ${DIR}/${BASE}-corrected.best.fas \
      > ${DIR}/${BASE}-corrected-smushed.best.fas

    python3 improve-msa.py --msa ${DIR}/${BASE}-corrected-smushed.best.fas \
      --max-cores 2 --word-size 12 > ${DIR}/${BASE}-smushed-improved.fa

    echo "working on ${DIR}/${BASE}-subreads"
    python3 make-consensus.py --msa ${DIR}/${BASE}-smushed-improved.fa \
      --read ${BASE} > ${DIR}/${BASE}-consensus.fa

    touch ${DIR}/${BASE}.done
  fi
}

# directory to extract subreads from and should include read# directory:
# for example barcode06/reads5
DIRR=$1
# number of processors to use
N=$2
# name of output file
out=$3

printenv PATH

for FILE in $DIRR/*/*-subreads.fq
  do
    ((i=i%N));
    ((i++==0)) && wait
    task $FILE&
  done

if test -f $out; then
    rm $out
fi
for FILE in $DIRR/*/*-consensus.fa
  do
    cat $FILE >> $out
  done

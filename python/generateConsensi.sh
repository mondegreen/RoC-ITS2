task(){
  BASE=`basename $1 -subreads.fq`
  DIR=`dirname $1`

  CONSENT-correct --in ${DIR}/${BASE}-subreads.fq \
    --out ${DIR}/${BASE}-corrected.fa --type ONT

  prank -F -o=${DIR}/${BASE}-corrected \
    -d=${DIR}/${BASE}-corrected.fa > /dev/null

  python smush-msa.py --msa ${DIR}/${BASE}-corrected.best.fas \
    > ${DIR}/${BASE}-corrected-smushed.best.fas

  python improve-msa.py --msa ${DIR}/${BASE}-corrected-smushed.best.fas \
    --max-cores 2 --word-size 12 > ${DIR}/${BASE}-smushed-improved.fa

  echo "working on ${DIR}/${BASE}-subreads"
  python make-consensus.py --msa ${DIR}/${BASE}-smushed-improved.fa \
    --read ${BASE} > ${DIR}/${BASE}-consensus.fa
}

N=25

for FILE in reads5/*/*-subreads.fq
  do
    ((i=i%N));
    ((i++==0)) && wait
    task $FILE&
  done
  cat reads5/*/*-consensus.fa > consensi-5plus-subreads.fa
done

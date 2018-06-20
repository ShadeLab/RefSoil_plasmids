#!/bin/bash -login

## This script searches protein sequence data using HMMs of antibiotic resistance genes
## This script does not involve assembly
## This overwrites previous search results

if [ $# -ne 2 ]; then
        echo "Requires two inputs : /path/hmmsearch_setenv.sh dataset"
        echo "  hmmsearch_setenv.sh is a file containing the parameter settings, requires absolute path."
        echo '  dataset should contain one or more protein sequence datasets to process with quotes around'
        echo 'Example command: /path/hmmsearch_setenv.sh "sequence1.faa sequence2.faa "'
        exit 1
fi

source $1
dataset=$2

for dataset in ${dataset}
	do 
		echo "### Search aa sequences ${dataset}"
		for gene in ${ARG}
			do
			${HMMSEARCH} --cpu ${THREADS} -E ${evalue} -o ${WORKDIR}/${gene}.${dataset}.stdout.txt --domtblout ${WORKDIR}/${gene}.${dataset}.${evalue}.tbl.txt -A ${WORKDIR}/${gene}.${dataset}.alignment.${evalue}.txt ${REF_DIR}/${gene}.hmm ${SEQDIR}/${dataset} || { echo "hmmsearch failed for ${gene} in ${dataset}" ; continue; }
			if [[ -s ${WORKDIR}/${gene}.${dataset}.alignment.${evalue}.txt ]]; then 
                		echo "${gene} in ${dataset} was detected"
        		else
                		echo "${gene} ${dataset}" >> ${WORKDIR}/no_hits.txt
                		rm ${WORKDIR}/${gene}.${dataset}.*
			fi
		done
	done

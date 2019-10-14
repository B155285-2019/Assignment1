#!/bin/bash

### My output folder destination

###FILES used during the script

path=$(pwd)
MY_FILE="$path"     
genome="/localdisk/data/BPSM/Assignment1/Tbb_genome/Tb927_genome.fasta.gz" 
bed_file="/localdisk/data/BPSM/Assignment1/Tbbgenes.bed"

### Arranging input fasta files according to the their TYPE
###Attaching type to the fasta files


##Assign types of the sample 
echo "================Types of Samples===================="
tn=$(awk '{IFS=$'\t';print $2}' /localdisk/data/BPSM/Assignment1/fastq/fqfiles | sort | uniq )
echo "$tn"
lines=$tn

fqfls=FASTA_FILES
chmod 700 ${MY_FILE}/
mkdir ${MY_FILE}/${fqfls}

cd /localdisk/data/BPSM/Assignment1/fastq
for file in ./*
do
	while read num type r1 r2
	do
	f="$(basename $file)"
	if [ $r1 == "$f" ]
		then 
		cp $file ${MY_FILE}/${fqfls}/${type}_$r1
	elif [ "$r2" == "$f" ]
		then
              	cp $file ${MY_FILE}/${fqfls}/${type}_$r2
	fi
	done < /localdisk/data/BPSM/Assignment1/fastq/fqfiles
done

cd ${MY_FILE}


###running fastqc
echo "==================FASTQ assessment and quality check=============="

fq_out=FASTQC_OUT
mkdir ${MY_FILE}/${fq_out}
fastqc -t 40 ${MY_FILE}/${fqfls}/* -o ${MY_FILE}/${fq_out}



###bowtie indexing and running
echo "====================Indexing the Reference Genome========================"

gn_name=$(basename "$genome")
cp $genome ${MY_FILE}/$gn_name
gunzip ${MY_FILE}/$gn_name
gn_n=${gn_name/.gz/}
gn_idx=${gn_name/_genome.fasta.gz/}
bowtie2-build --threads 40 -f ${MY_FILE}/$gn_n ${MY_FILE}/$gn_idx

echo "=====================Alignment using Bowtie2==================================="

declare -a fname_arr=()

for file in ${MY_FILE}/${fqfls}/*;
do
	name=$(basename "$file")
	fname=${name/_[0-9].fq.gz/} 
	if [[ ! "${fname_arr[@]}" =~ "${fname}" ]];
	then
		fname_arr+=("${fname}")
	fi	
done


bamout=BAM_OUT
mkdir ${MY_FILE}/${bamout}

for fname in "${fname_arr[@]}";
do
	bowtie2 -p 60 -q -x ${MY_FILE}/$gn_idx -1 ${MY_FILE}/${fqfls}/${fname}_1.fq.gz -2 ${MY_FILE}/${fqfls}/${fname}_2.fq.gz -S ${MY_FILE}/${bamout}/${fname}.sam
	samtools view -b ${MY_FILE}/${bamout}/${fname}.sam > ${MY_FILE}/${bamout}/${fname}.bam
	samtools sort -o ${MY_FILE}/${bamout}/${fname}.st.bam -@ 40 ${MY_FILE}/${bamout}/${fname}.bam
	samtools index ${MY_FILE}/${bamout}/${fname}.st.bam
done

echo "=============================Generating gene counts data=============================="

bed_out=bed_OUT
mkdir ${MY_FILE}/${bed_out}

for line in $lines;
do 	
	bedtools multicov -bams ${MY_FILE}/${bamout}/${line}_*.st.bam -bed $bed_file > ${MY_FILE}/${bed_out}/${line}.bed
done


out=RESULT
mkdir ${MY_FILE}/${out}

declare -a type_arr=()
for line in $lines;
do
	type_arr+=("${line}")
done

for line in $lines;
do
	while read chr start end gene term strand score1 score2 score3
	do
		sum=$(( ${score1} + ${score2} + ${score3} ))
		ave=$(( ${sum}/3 ))
		echo -e "${gene}" >> ${MY_FILE}/${out}/${line}_gene_list.txt
		echo -e "${ave}" >> ${MY_FILE}/${out}/${line}_gene_counts.txt
	done < ${MY_FILE}/${bed_out}/${line}.bed
done


paste ${MY_FILE}/${out}/${type_arr[1]}_gene_list.txt ${MY_FILE}/${out}/${type_arr[0]}_gene_counts.txt ${MY_FILE}/${out}/${type_arr[1]}_gene_counts.txt > ${MY_FILE}/${out}/combined_genes_ave_count.txt

head=(Genes ${type_arr[0]}_mean ${type_arr[1]}_mean) 
(IFS=$'\t'; echo "${head[*]}"; cat ${MY_FILE}/${out}/combined_genes_ave_count.txt) > ${MY_FILE}/${out}/gene_counts.txt


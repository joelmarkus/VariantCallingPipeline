#!/bin/bash

#############################################################################
#NOTE: READ BEFORE RUNNING!
#Current jar location of GATK in the script is ~/bin/GenomeAnalysisTK.jar. Change according to the location of your jar file.
#Added no. of threads as an argument. Use -t [INT].
#############################################################################

usage="Running $(basename "$0")... 


[-a str] [-b str] [-r str] [-e int] [-o str] [-f str] [-z int] [-v str] [-i int] [-t int]


-- Wrapper tool to generate variant calling .vcf and .bed files from raw reads.


WARNING: 
	1) Please specify the entire path of your files.
	2) The jar location of GATK should be ~/bin/GenomeAnalysisTK.jar.


Options:
   INPUT LOCATIONS:
	-a paired end read 1 location. [required]
	-b paired end read 2 location. [required]
	-r reference genome location. [required]
   
   
   ADVANCED OPTIONS:
	-o output file name.
	-f mills file location.
	
	-e realignment option (1 = re-alignment, 0 = no re-alignment) [0]
	-z file compression (1 = compressed output vcf, 0 = no compression) [0]
	-i BAM indexing (1 = index BAM, 0 = no indexing) [0]
	
	-v verbose (1 = verbose, 0 = silent) [0]
	
	-t no. of threads [1] 
	-k echo inputs (1 = yes, 0 = no) [0]
	"
    
 
#Input arguments
e=0;z=0;v=0;i=0;k=0;t=1;
while getopts “a:b:r:o:f::e::z::v::i::k::t::h” option 
do 
		case $option in 
			h) echo "$usage";exit;;
			a) reads1=$OPTARG;;
			b) reads2=$OPTARG;;
			r) ref=$OPTARG;;
			o) output=$OPTARG;;
			f) millsFile=$OPTARG;;
			e) realign=$OPTARG;;
			z) gunzip=$OPTARG;;
			v) v=$OPTARG;;
			i) index=$OPTARG;;
			k) answer=$OPTARG;;
			t) threads=$OPTARG;;
		esac 
done 

e=$realign;z=$gunzip;i=$index;k=$answer;t=$threads #for the -h argument printing
#h input: to print the arguments it took in

#to print the arguments the user has provided
if [[ $answer -eq 1 ]]
then
	echo -e "Arguments taken:
	-a $reads1 paired end read 1 input successful.
	-b $reads2 paired end read 2 input successful.
	-r $ref reference genome input successful.
	-e $e submitted (1 = re-alignment, 0 = no re-alignment).
	-o $output output file.
	-f $millsFile mills file location.
	-z $z submitted (1 = compressed output vcf, 0 = no compression).
	-v $v submitted (1 = verbose, 0 = no verbose).
	-i $i submitted (1 = index BAM, 0 = no indexing).
	-t no. of threads: $t. \n"
fi 

#file checks
if [[ ! -f $ref ]]; then
    echo -e "ERROR: Reference genome submitted doesn't exist. Exiting... \n"; exit 1;
elif [[ ! -f $reads1 ]]; then
    echo -e "ERROR: Read file (Read 1) doesn't exist. Exiting... \n"; exit 1;
elif [[ ! -f $reads2 ]]; then
    echo -e "ERROR: Read file (Read 2) doesn't exist. Exiting... \n"; exit 1;
fi

#read vcf file
if [[ -f $output ]]; then
	read -p "vcf file exits, type \"overwrite\" to overwirte the file or \"exit\" to exit the program" vcf_overwrite
	if [[ $vcf_overwrite -eq "exit" ]]
		then
		echo "Exiting..."; exit 1;
	elif [[ $vcf_overwrite -eq "overwrite" ]]
		then
		echo "VCF file will be overwritten"
		rm ${output}.vcf.gz
	else
		echo "Please supply a correct string. Exiting..."; exit 1;
	fi
fi	

#Create a path string for everything
PathStr=$(echo $reads1 | sed -r -e 's/\.fq(\.gz|)$//g')

#create index
if [[ $v -eq 1 ]]; then echo "Indexing the reference genome..."; fi
bwa index $ref
if [[ $v -eq 1 ]]; then echo "Indexing done."; fi

#creating SAM file
if [[ $v -eq 1 ]]; then echo "Mapping reads..."; fi
bwa mem -t $threads -R '@RG\tID:foo\tSM:bar\tLB:library1' $ref $reads1 $reads2 > ${PathStr}.sam
if [[ $v -eq 1 ]]; then echo "Mapping done...SAM present in ${PathStr}.sam"; fi

#conversion from sam to bam
if [[ $v -eq 1 ]]; then echo "Creating BAM from the SAM file..."; fi
samtools fixmate -O bam --threads $threads ${PathStr}.sam ${PathStr}.bam
if [[ $v -eq 1 ]]; then echo "BAM created and present in ${PathStr}.bam"; fi

#BAM sort 
if [[ $v -eq 1 ]]; then echo "Sorting initiated..."; fi
mkdir ${PathStr}_temp
samtools sort -O bam --threads $threads -o ${PathStr}.sorted.bam -T ${PathStr}_temp/ $PathStr.bam 
if [[ $v -eq 1 ]]; then echo "Sorted file present in ${PathStr}.sorted.bam"; fi

#create fai and dict
if [[ $v -eq 1 ]]; then echo "Creating .fai and .dict files"; fi
dict=$(echo $ref | sed -E 's/.fa(|sta)/.dict/g')
samtools faidx $ref -o $ref.fai
samtools dict $ref -o $dict
if [[ $v -eq 1 ]]; then echo ".fai and .dict files created"; fi

#realignment/improvement
if [[ $realign -eq 1 ]]; then
if [[ $v -eq 1 ]]; then 
echo "Running re-alignment, this may take a while..."; fi
samtools index ${PathStr}.sorted.bam
java -Xmx2g -jar ~/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref -I ${PathStr}.sorted.bam -o $PathStr.intervals --known $millsFile -log ${PathStr}.RTC.log
java -Xmx4g -jar ~/bin/GenomeAnalysisTK.jar -T IndelRealigner -R $ref -I ${PathStr}.sorted.bam -targetIntervals ${PathStr}.intervals -o ${PathStr}.sorted.bam -log ${PathStr}.IR.log
rm ${PathStr}.sorted.bai
elif [[ $realign -eq 0 ]]; then true ; else echo "ERROR: Please supply either 0 or 1 for -e. Exiting... "; exit 1; fi 
if [[ $v -eq 1 ]]; then echo "Re-alignment done"; fi

#index bam
if [[ $index -eq 1 ]];	then
if [[ $v -eq 1 ]]; then echo "Indexing BAM..."; fi
samtools index ${PathStr}.sorted.bam
elif [[ $index -eq 0 ]]; then true ; else echo "ERROR: Please supply either 0 or 1 for -i. Exiting... "; exit 1; fi 
if [[ $v -eq 1 ]]; then echo "BAM indexed"; fi

#variant calling
if [[ $v -eq 1 ]]; then echo "Calling variants..."; fi
bcftools mpileup -Ou -f $ref ${PathStr}.sorted.bam --threads $threads| bcftools call -vmO z -o ${output}.vcf.gz --threads $threads
if [[ $v -eq 1 ]]; then echo "VCF file generated, present as ${output}.vcf.gz"; fi

#index using tabix
if [[ $v -eq 1 ]]; then echo "Indexing VCF..."; fi
tabix -p vcf ${output}.vcf.gz
if [[ $v -eq 1 ]]; then echo "Converting vcf to bed..."; fi

#converting vcf to bed
gunzip -k ${output}.vcf.gz
cat ${output}.vcf | awk '!/#/ {print $1,$2,length($5)-length($4)+$2,length($5)-length($4)}' | sed -e 's/chr//g' | tee ${output}.bed | awk '{ if($4==0) { print $0 > "snps.txt" } else { print $0 > "indels.txt" } }'
if [[ $v -eq 1 ]]; then echo "bed file created with INDEL and SNP files."; fi
if [[ $gunzip -eq 1 ]]; then	rm ${output}.vcf; elif [[ $gunzip -eq 0 ]]; then rm ${output}.vcf.gz; else echo "ERROR: Please supply either 0 or 1 for -z. Exiting... "; exit 1; fi

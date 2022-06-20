#!/bin/bash

#Default Values 
PLOTS="./user_input_files/plot_data.csv"
PROBES="./user_input_files/probes.bed"
META_DATA="./user_input_files/data_metadata.csv"
NANO_BLOT_RSCRIPT="./scripts/nano_blot_generation.R"
PRINT_HELP=FALSE
SUBSET_BAMS=TRUE
MAKE_PLOT=TRUE
CDNA=FALSE



while getopts ":HFPCR:M:B:T:" opt; do
  case ${opt} in
    H ) 
    PRINT_HELP=TRUE
      ;;
    F ) 
    SUBSET_BAMS=FALSE
      ;;
    P ) 
    MAKE_PLOT=FALSE
      ;;
    C ) 
    CDNA=TRUE
    echo "THIS AINT DO SHIT RIGHT NOW"
      ;;
    R ) 
    	NANO_BLOT_RSCRIPT=${OPTARG}
    	;;
    M ) 
    	META_DATA=${OPTARG}
    	;;
    B ) 
    	PLOTS=${OPTARG}
    	;;
    T ) 
    	PROBES=${OPTARG}
    	;;
    \? ) echo "Usage: "
      ;;
  esac
done

echo "R Script: $NANO_BLOT_RSCRIPT";
echo "Meta Data File: $META_DATA";
echo "Probes Bed File: $PROBES";
echo "Plots File: $PLOTS";

#Check if user asked for the help text and echo the help text.
if [[ "$PRINT_HELP" == TRUE ]]
then
	echo "
Nanoblot (Version 1.0)

O       o O       o O       o O       o
| O   o | | O   o | | O   o | | O   o |
| | O | | | | O | | | | O | | | | O | |
| o   O | | o   O | | o   O | | o   O |
o       O o       O o       O o       O
             xxxxxxx
         x xxxxxxxxxxxxx x
      x     xxxxxxxxxxx     x
            xxxxxxxxx
  x          xxxxxxx          x
              xxxxx
 x             xxx             x
                x
xxxxxxxxxxxxxxx   xxxxxxxxxxxxxxx
 xxxxxxxxxxxxx     xxxxxxxxxxxxx
  xxxxxxxxxxx       xxxxxxxxxxx
   xxxxxxxxx         xxxxxxxxx
     xxxxxx           xxxxxx
       xxx             xxx
            x         x
                 x
O       o O       o O       o O       o
| O   o | | O   o | | O   o | | O   o |
| | O | | | | O | | | | O | | | | O | |
| o   O | | o   O | | o   O | | o   O |
o       O o       O o       O o       O

For an explanation of the required input files see the README.md

==========================================
-H  |  Print help menu
-T  |  Probes bed file 
-B  |  Blots metadata file
-M  |  Location of metadata file
-R  |  Use custem R script
-N  |  Normalize Data (Work in progress)
-C  |  Treat reads as cDNA (disregard strand) (WIP)
-F  |  Skip subsetting BAM files for plot generation
-P  |  Skip nanoblots generation
==========================================
"

exit
fi

declare -i END_PLOT=$(wc -l < $PLOTS)

for (( c=2; c<=$END_PLOT; c++ ))
do
  P_LINE=$(head -n $c $PLOTS | tail -n -1)
  echo "=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~="
  PLOT_NUM=$((c-1))
  echo "Generating plot" $PLOT_NUM
  echo "======="
  
  IFS=$'\t'; read -a fields <<<"$P_LINE"
  
  echo 'Getting Probe:'${fields[2]}
  echo 'Duplication Factor:'${fields[3]}
  
  DUP_FACTOR=${fields[3]}
  TARGET=${fields[2]}
  BAMS=${fields[1]}
  
  declare -i END_PROBE=$(wc -l < $PROBES)
  END_PROBE=$((END_PROBE+1))
  
  for (( d=1; d<=$END_PROBE; d++ ))
  do
     if [[ $d == $END_PROBE ]]
     then
       echo "Probe not found. Check bed file and blots metadata file."
       break
     fi
     
     PR_LINE=$(head -n $d $PROBES | tail -n -1)
     
     IFS=$'\t'; read -a feels <<<"$PR_LINE"
     TARGET_PROBE=${feels[3]}
     
     if [[ "$TARGET_PROBE" == "$TARGET" ]]
     then
       echo ${fields[2]}": "${feels[0]}:${feels[1]}-${feels[2]}

       declare -i END_META=$(wc -l < $META_DATA)
       IFS=','; read -a samples <<<"$BAMS"
      
       if [[ "$SUBSET_BAMS"  == TRUE ]]
       then
	       for (( e=2; e<=$END_META; e++ ))
	       do
	         DATA_LINE=$(head -n $e $META_DATA | tail -n -1)
	         IFS=$'\t'; read -a EELS <<<"$DATA_LINE"
	         DATA_LOCATION=${EELS[1]}
	         TEMP_NAME=${EELS[0]}_$TARGET_PROBE.bam
	         echo "Subsetting: "$DATA_LOCATION"
Naming Subset: " $TEMP_NAME
					 
					 #Get number of reads and write to a file.
					 if [[ "$NORM_DATA"  == TRUE ]]
					 then
						READ_NUM=$(samtools view -c -F 260 $DATA_LOCATION) #This line should output the number of reads in a bam file
						echo "Number of reads in orgnial file: "$READ_NUM
					 fi
					 
					 #Check if data should be treated as strand specific
					 if [[ "$CDNA"  == TRUE ]]
					 then
			         samtools view -b $DATA_LOCATION ${feels[0]}:${feels[1]}-${feels[2]} > "./temp/$TEMP_NAME"
					 else
							 echo "Strand: "${feels[5]}
							 echo "Yell at sam to figure out how to filter by strand"
							 #There gunna have to be an if statement here checking the strand of the probe
							 if [[ "${feels[5]}"  == "+" ]]
							 then
							   samtools view -F 20 -b $DATA_LOCATION ${feels[0]}:${feels[1]}-${feels[2]} > "./temp/$TEMP_NAME"
							 elif [[ "${feels[5]}"  == "-" ]]
							 then
							   samtools view -f 16 -b $DATA_LOCATION ${feels[0]}:${feels[1]}-${feels[2]} > "./temp/$TEMP_NAME"
							 else
							   echo "Warning: No strand found for probe treating region disregarding strand information"
							   samtools view -b $DATA_LOCATION ${feels[0]}:${feels[1]}-${feels[2]} > "./temp/$TEMP_NAME"
							 fi
					 fi
	       done
       else 
         echo "Skipping filtering BAM files. If filtering is desired remove -F flag."
       fi
       break
     fi
  done
  
  #Normalize Reads if desired
  
  if [[ "$MAKE_PLOT" == TRUE ]]
  then
  	echo ""
  	echo "======="
  	echo "Running R script"
  	echo "======="
  	Rscript $NANO_BLOT_RSCRIPT $BAMS $TARGET $DUP_FACTOR
  	echo "======="
  	echo "======="
  else
  	echo "Skipping plot generation. If plot generation is desired remove -P flag."
  fi
echo "=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~="
echo ""
done 


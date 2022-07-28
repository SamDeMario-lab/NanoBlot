#!/bin/bash

#Default Values 
PLOTS="./user_input_files/plot_data.csv" #Be careful of what these variables equal, since the value of
# these variables are still strings of file locations 
PROBES="./user_input_files/probes.bed"
META_DATA="./user_input_files/data_metadata.csv"
NANO_BLOT_RSCRIPT="./scripts/nano_blot_generation.R"
PRINT_HELP=FALSE
SUBSET_BAMS=TRUE
MAKE_PLOT=TRUE
CDNA=FALSE
NORM=TRUE
CLEAN_ALL=FALSE

#The option-string tells getopts which options to expect and which of them must have an argument. 
#Every character is simply named as is, and when you want it to expect an argument, just place a colon
# after the option flag 
# The first colon basically means getopts switches to "silent error reporting mode" --> allows you to 
# handle errors yourself without being disturbed by annoying messages 
while getopts ":HFPCNWR:M:B:T:" opt; do
	case $opt in
		H ) 
		PRINT_HELP=TRUE
		echo "Help printed"
			;;
		F ) 
		SUBSET_BAMS=FALSE
		echo "Subsetting BAM files skipped"
			;;
		P ) 
		MAKE_PLOT=FALSE
		echo "Nanoblot generation skipped"
			;;
		C ) 
		CDNA=TRUE
		echo "Treating reads as cDNA"
			;;
		R ) 
		NANO_BLOT_RSCRIPT=$OPTARG
		echo "Using custom R script $NANO_BLOT_RSCRIPT"
			;;
		M ) 
		META_DATA=$OPTARG
		echo "Using location of metadata file $METADATA"
			;;
		B ) 
		PLOTS=$OPTARG
		echo "Using blots metadata file $PLOTS"
			;;
		T ) 
		PROBES=$OPTARG
		echo "Using probes bed file $PROBES"
			;;
		N ) 
		NORM=FALSE
		echo "Skipping data normalization"
			;;
		W )
		CLEAN_ALL=TRUE
		echo "Clear all files from temp after plot generation"
			;;
		\? ) echo "Invalid option: -$OPTARG"
			;;
	esac
done

echo -e "Starting Nanoblot" '\u2622' 
echo "======="
sleep 1 #Pauses for one second to actually let the user know that the program is being started
echo "R Script: $NANO_BLOT_RSCRIPT";
echo "Meta Data File: $META_DATA";

#using the debugger

declare -i END_META=$(wc -l < $META_DATA) #stores the number of lines of META_DATA into a variable
# that is declared called END_META, declare -i basically makes it so that the variable can only be
# changed to another integer, declares the variable type 

echo "Probes Bed File: $PROBES";
echo "Plots File: $PLOTS";

#Check if user asked for the help text and echo the help text.
if [ $PRINT_HELP = TRUE ]
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
-N  |  Skip data normalization (Work in progress)
-C  |  Treat reads as cDNA (disregard strand) 
-F  |  Skip subsetting BAM files for plot generation
-P  |  Skip nanoblots generation
-W  |  Clear all files from ./temp/ after plot generation
==========================================
"
exit
fi

declare -i END_PLOT=$(wc -l < $PLOTS) # Again, declares the line numbers of PLOTS and sets it to 
# a variable called END_PLOT 

if [ $NORM = TRUE ] #edited from the brute force method 
then
	echo -e "=======\nNormalization"
	for (( c=2; c<=$END_PLOT; c++ ))
	do
		P_LINE=$(head -n $c $PLOTS | tail -n -1) #This line gets the individual row for each row of the plot_data
		IFS=$'\t'; # This sets the internal field separator to tab, which is used to determine
		# word splitting --> basically overrides the natural word boundary of a whitespace
		read -a fields <<<"$P_LINE" # I think this code creates an array called fields which separates each 
		# column of the individual row of $PLOTS 
		BAMS=${fields[1]} # This then gets the 1st index, 2nd column of each row? which is the loading order?
		NORM_FOLDER="./temp/"$(echo "$BAMS" | sed -e 's/,/_/g')"_NORM" #This basically creates a new folder called norm
		# and then changes the comma in the loading order to an underscore, the s stands for substitute function, whereas
		# the g stands for a global change, since s by itself only applies to the first match
		if [ ! -d "$NORM_FOLDER" ] # checks if the norm folder exists and is a directory 
		then
			echo $NORM_FOLDER " not found" #prints that the normalization folder is not found
			mkdir -p $NORM_FOLDER; #The manual says that this creates intermediate directories as required. If this option is not specified, 
			# the full path prefix of each operant must already exist --> this allows it so that you can create sub directories of a directory
			# It will create the parent directory first, and if it does exist, it will move to create further sub-directories 
			NORM_METADATA_FILE=$NORM_FOLDER"/data_metadata.csv" #Creates the metadata file type in the norm folder 
			NORM_METADATA_HEADER="Sample_name (This must be unique for each sample)  Type (FAST5 or BAM)  Location (For BAM inputs the path to the bam file should be given.)"
			# Above code is the header 
			echo "Starting Normalization"
			echo "Samples: "$BAMS
			IFS=',' #Changes the input field separator to commas now 
			echo ${NORM_METADATA_HEADER} >$NORM_METADATA_FILE # This creates the first line to be inputted into the empty csv file, which 
			# basically only contains the header first 
			COMPARE=999999999999999 #This is probobly not the best way to do this, agreed, this is hard coding
			for f in $BAMS; #Remember, the IFS changed to commas, so now this will create a for loop run for each of the loading orders
			do
				for (( o=2; o<=$END_META; o++ )) #Nested for loop, feel like there might be some easier way to do this as well 
				do
					DATA_LINE_o=$(head -n $o $META_DATA | tail -n -1) #Takes single rows again but this time out of the 
					# metadata folder --> basically now it knows how to pull the BAM files 
					IFS=$'\t'; #changes the IFS again to tabs
					read -a EELS <<<"$DATA_LINE_o" #Creates a new array called EELS --> this basically is each row in the metadata file
					DATA_LOCATION=${EELS[1]} #Takes the 1st index, which is the value for the data location
					SAMP_NAME=${EELS[0]} # Takes the 0th index, which is the value for the sample name
					
					#This if statement can just be calculated in one awk command
					#awk -v var="$f" '$1==var {print $0}'
					
					if [[ $f == $SAMP_NAME ]] #If it finds the sample name, so like WT or RRP6 in the SAMP_NAME of the data metadata
					then
						echo "Found "$f 
						READ_COUNT=$(samtools view -c -F 260 $DATA_LOCATION) # samtools view function will view and convert SAM/BAM/CRAM files
						# the -c option counts the number of alignments instead of printing out the alignments 
						# the -F flag does not output any alignments with the bits set in the Flag argument 
						# a -F flag of 260 means that
						# 256 - The given alignment is a secondary alignment 4 - The read was not aligned
						# In essence, this makes it so that all unaligned and non-first choice mappings are removed
						# In the future, we need a better way to store the sample files for reference 
						if [[ $READ_COUNT -lt $COMPARE ]] # -lt is less than
						then
							COMPARE=$READ_COUNT # An easier way is to just declare the variable, I get that you want this variable to 
							# exist outside the scope of the for loop 
							# Will take the lowest read count if there are multiple read counts from multiple bam file locations for a single
							# sample name 
						fi
					fi
				done
			done
			echo $COMPARE # Prints out the lowest read count across all samples in one plot
			IFS=',' #Changes the IFS again
			
			#After it finds all the metadata it needs to normalize the loading_order, it then normalizes the metadata sorted
			# bam files to that normalized number 
			for f in $BAMS;
			do
				for (( e=2; e<=$END_META; e++ ))
				do
					DATA_LINE=$(head -n $e $META_DATA | tail -n -1)
					IFS=$'\t'; read -a EELS <<<"$DATA_LINE"
					DATA_LOCATION=${EELS[1]}
					SAMP_NAME=${EELS[0]}
					if [[ $f == $SAMP_NAME ]]
					then
						echo "Subsetting "$SAMP_NAME" to "$COMPARE" reads."
						OUTPUT_NORM_SAM=$NORM_FOLDER"/"$f"_NORM.sam"
						OUTPUT_NORM_BAM=$NORM_FOLDER"/"$f"_NORM.bam"
						cat <(samtools view -H $DATA_LOCATION) <(samtools view ${DATA_LOCATION} | shuf -n $COMPARE) > $OUTPUT_NORM_SAM #THIS TAKES A LOT OF RAM
						# above line, -H means to output the header only, I think this is just doing header and body separately? 
						# so the way this is normalizing is it just randomly takes the lowest number of reads across all sequenced files
						# There has to be a better way to normalize using a package
						samtools view -S -b $OUTPUT_NORM_SAM > $OUTPUT_NORM_BAM #-S, previously required if input was in SAM format
						# now, it automatically detects -b outputs in bam format 
						samtools sort $OUTPUT_NORM_BAM -o $OUTPUT_NORM_BAM #sorts, which sorts all the reads by genomic position
						samtools index $OUTPUT_NORM_BAM #indexes the bam file, mainly for igv use, allows navigating it without
						# loading all the file into memory
						echo $SAMP_NAME"		"$OUTPUT_NORM_BAM >> $NORM_METADATA_FILE
					fi
				done
			done
		else #else loop runs if the norm directory already exists 
			echo "" 
			echo $NORM_FOLDER" already exists."
			echo "If data normalization should be re-run delete "$NORM_FOLDER
			echo ""
		fi #finishes the if to check if the directory exists already for normalization
	done # finishes the first loop
fi #finishes the first if to check if normalization is already done 

# Probe data + configuring plots 
for (( j=2; j<=$END_PLOT; j++ )) # For loop that goes through each row of plot_data
do
	P_LINE=$(head -n $j $PLOTS | tail -n -1) #Gets each individual row 
	echo "=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~="
	PLOT_NUM=$((j-1)) #Sets an integer that is the plot_num, starting at plot_num 1
	echo "Generating plot" $PLOT_NUM
	echo "======="
	
	IFS=$'\t'; read -a fields <<<"$P_LINE"
	
	echo 'Probe:'${fields[2]}
	if [ -z "${fields[4]}" ] #checks to see if there is an antiprobe, -z checks for empty string 
	then
		echo "No Negative Probe Used" #Print statement
	else
		echo 'Negative Probe:'${fields[4]} #Print statement
	fi
	echo 'Duplication Factor:'${fields[3]} #Print statement
	echo "======="
	
	# Setting variables 
	DUP_FACTOR=${fields[3]}
	TARGET=${fields[2]}
	TARG_NEGS=${fields[4]}
	BAMS=${fields[1]}
	
	# If normalization was true, then set the meta_data variable to the normalized meta_data 
	if [ $NORM = TRUE ]
	then
		NORM_FOLDER="./temp/"$(echo "$BAMS" | sed -e 's/,/_/g')"_NORM"
		NORM_METADATA_FILE=$NORM_FOLDER"/data_metadata.csv"
		META_DATA=$NORM_METADATA_FILE
	fi

#Apply the first filter
	declare -i END_PROBE=$(wc -l < $PROBES)
	END_PROBE=$((END_PROBE+1)) # I don't actually know if this +1 to the END_PROBE integer is necessary?

	for (( d=1; d<=$END_PROBE; d++ ))
	do
		if [[ $d == $END_PROBE ]] # Checks to see if the end_probe file is an empty text file 
		then
			echo "Probe not found. Check bed file and blots metadata file."
		fi

		IFS=$','
		PR_LINE=$(head -n $d $PROBES | tail -n -1)

		IFS=$'\t'; read -a feels <<<"$PR_LINE"
		TARGET_PROBE=${feels[3]}

		if [[ "$TARGET_PROBE" == "$TARGET" ]]
		then
			echo 'Probe' #Note to self: I guess there is no check for antisense probe direction yet
			echo ${fields[2]}": "${feels[0]}:${feels[1]}-${feels[2]} #These array subsets do not belong to the 
			# same variables 

			IFS=$','; echo $PR_LINE > "./temp/temp_bed.bed" #Overwrites the current matched probe_line and 
			# puts it into a temp bed file

			if [[ "$SUBSET_BAMS"  == TRUE ]]
			then
				for (( t=2; t<=$END_META; t++ )) 
				do
					DATA_LINE_T=$(head -n $t $META_DATA | tail -n -1)
					IFS=$'\t'; read -a EELS <<<"$DATA_LINE_T"
					DATA_LOCATION=${EELS[1]}
					TEMP_NAME=${EELS[0]}_$TARGET_PROBE.bam
					echo "Subsetting: "$DATA_LOCATION"
Naming Subset: " $TEMP_NAME
					if [[ "$CDNA"  == TRUE ]]
					then
						bedtools intersect -a $DATA_LOCATION -b "./temp/temp_bed.bed" -wa -split -nonamecheck > "./temp/$TEMP_NAME"
						samtools index "./temp/$TEMP_NAME"
					else
						bedtools intersect -a $DATA_LOCATION -b "./temp/temp_bed.bed" -wa -split -s -nonamecheck > "./temp/$TEMP_NAME"
						samtools index "./temp/$TEMP_NAME"
					fi
				done
			else 
				echo "Skipping filtering BAM files. If filtering is desired remove -F flag."
			fi
				break
		fi
	done
	echo "======="
	if [ -z "${fields[4]}" ]
	then
		echo "No negative probe used"
	else
		echo 'Negative Probe'

		declare -i END_PROBE=$(wc -l < $PROBES)
		END_PROBE=$((END_PROBE+1))

		for (( h=1; h<=$END_PROBE; h++ ))
		do
			if [[ $h == $END_PROBE ]]
			then
				echo "Probe not found. Check bed file and blots metadata file."
				break
			fi

			IFS=$','
			APR_LINE=$(head -n $h $PROBES | tail -n -1)
			echo $APR_LINE > "./temp/temp_anti_bed.bed"

			IFS=$'\t'; read -a deels <<<"$APR_LINE"
			TARG_ANTI=${deels[3]}

			if [[ $TARG_NEGS == $TARG_ANTI ]]
			then
				echo $TARG_ANTI": "${deels[0]}:${deels[1]}-${deels[2]}
				if [[ "$SUBSET_BAMS"  == TRUE ]]
				then
					for (( p=2; p<=$END_META; p++ ))
					do
						DL_ANTI=$(head -n $p $META_DATA | tail -n -1)
						IFS=$'\t'; read -a bells <<<"$DL_ANTI"
						FIRST_NAME="./temp/"${bells[0]}"_"$TARGET".bam"
						echo "Subsetting:" $FIRST_NAME
						NEG_NAME="./temp/"${bells[0]}"_"$TARGET"_anti_"$TARG_ANTI".bam"
						echo "Naming Subset: "$NEG_NAME
						if [[ "$CDNA"  == TRUE ]]
						then
							bedtools intersect -a $FIRST_NAME -b "./temp/temp_anti_bed.bed" -wa -split -v -nonamecheck > $NEG_NAME
							samtools index $NEG_NAME
						else
							bedtools intersect -a $FIRST_NAME -b "./temp/temp_anti_bed.bed" -wa -split -v -s -nonamecheck > $NEG_NAME
							samtools index $NEG_NAME
						fi
					done
				fi
			break
			fi
		done
	fi
	echo "======="
	if [[ "$MAKE_PLOT" == TRUE ]]
	then
		echo "======="
		echo "Not running R script for now to trace all bash scripts"
		BAMS=${fields[1]} #I dont know why I need this but I do
		Rscript $NANO_BLOT_RSCRIPT $BAMS $TARGET $DUP_FACTOR $TARG_NEGS
		echo "======="
	else
		echo "Skipping plot generation. If plot generation is desired remove -P flag."
	fi
done 
echo "=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~="


# This code essentially just removes all the temp files that were created by the program
if [ $CLEAN_ALL = TRUE ]
then
	echo "Clearing ./temp/"
	rm -r ./temp/*
fi
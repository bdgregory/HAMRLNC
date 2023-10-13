#!/bin/bash
set -u

# Harry Li, University of Pennsylvania & Chosen Obih, University of Arizona


usage () {
    echo ""
    echo "Usage : sh $0"
    echo ""

cat <<'EOF'
  
  ######################################### COMMAND LINE OPTIONS #############################
  REQUIRED:
    -o	<project directory>
    -t	<SRA accession list.txt or folder of raw fastq files>
    -c	<filenames for each fastq.csv>
    -g	<reference genome.fa>
    -i  <reference genome annotation.gff3>
    -l	<read length>
    -s	<genome size in bp >
    -e	<genome annotation generator>

  OPTIONAL: 
    -n  number of threads (default 4)
    -a	[use TopHat2 instead of STAR]
    -b	[Tophat library choice: fr-unstranded, fr-firststrand, fr-secondstrand]
    -f	[filter]
    -m	[HAMR model]
    -k  [activate hamrbox]
    -p  [activate evolinc_i]
    -u  [activate featurecount]
    -v  [evolinc_i_option: M or MO, default=M]
    -Q	[HAMR: minimum qualuty score, default=30]
    -C	[HAMR: minimum coveragem default=50]
    -E	[HAMR: sequencing error, default=0.01]
    -P	[HAMR: maximum p-value, default=1]
    -F	[HAMR: maximum fdr, default=0.05]
    -T  <transposable Elements file> (optional file for evolinc_i)
    -G  <CAGE RNA file> (optional file for evolinc_i)
    -D  <known lincRNA file> (optional file for evolinc_i)
    -h	[help message] 


  ################################################# END ########################################
EOF
    exit 0
}

#curdir=$(dirname "$0")
threads=4
tophat=false
quality=30
coverage=50
err=0.01
pvalue=1
fdr=0.05
evolinc_i_option="M"
tophatlib="fr-firststrand"
filter=$util/filter_SAM_number_hits.pl
model=$util/euk_trna_mods.Rdata
evolinc_i=false
featurecount=false
hamrbox=false
generator=""

#############Grabbing arguments############
while getopts ":o:t:c:g:i:z:l:b:e:v:s:n:fmhQCakTGDupEPF:" opt; do
  case $opt in
    o)
    out=$OPTARG # project output directory root
     ;;
    t)
    acc=$OPTARG # SRA accession 
     ;;
    c)
    csv=$OPTARG # SRR to filename table
     ;;
    g)
    genome=$OPTARG # reference genome directory
     ;;
    i)
    annotation=$OPTARG # reference genome annotation
    ;;
    l)
    length+=$OPTARG # read length 
     ;;
    e)
    generator=$OPTARG # organism abbreviation for annotationGenerate
    ;;
    s)
    genomelength=$OPTARG # length or size of the genome
     ;;
    f)
    filter=$OPTARG
     ;;
    m)
    model=$OPTARG
     ;;
    v)
    evolinc_i_option="$OPTARG"
    ;;
    n)
    threads=$OPTARG
    ;;
    p)
    evolinc_i=true
    ;;
    k)
    hamrbox=true
    ;;
    u)
    featurecount=true
    ;;
    Q)
    quality=$OPTARG
    ;;
    T)
    blast_file=$OPTARG
     ;;
    G)
    cage_file=$OPTARG
     ;;
    D)
    known_linc=$OPTARG
    ;;
    C)
    coverage=$OPTARG
    ;;
    b)
    tophatlib=$OPTARG
    ;;
    E)
    err=$OPTARG
    ;;
    a)
    tophat=true
    ;;
    P)
    pvalue=$OPTARG
    ;;
    F)
    fdr=$OPTARG
    ;;
    h)
    usage
    ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# reassign sample input files, genome and annotation files name and include file paths
user_dir=$(pwd)
genome="$user_dir"/"$genome"
annotation="$user_dir"/"$annotation"
out="$user_dir"/"$out"
acc="$user_dir"/"$acc"
csv="$user_dir"/"$csv"
# reassign the below only if inputed
if [ ! -z "$blast_file" ]; then blast_file="$user_dir"/"$blast_file"; fi
if [ ! -z "$cage_file" ]; then cage_file="$user_dir"/"$cage_file"; fi
if [ ! -z "$known_file" ]; then known_linc="$user_dir"/"$known_linc"; fi

# assigning additional variables
dumpout=$out/datasets
mismatch=$((length*6/100))
overhang=$((mismatch-1))
genomedir=$(dirname "$genome")
last_checkpoint=""

# Assigning the appropriate annotationGenerate.R 
if [[ $generator == "AT" ]]; then
    generator=$annotationGenerate/annotationGenerateAT.R
    echo "Model organism detected: Arabidopsis thaliana"
elif [[ $generator == "BD" ]]; then
    generator=$annotationGenerate/annotationGenerateBD.R
    echo "Model organism detected: Brachypodium distachyon"
elif [[ $generator == "ZM" ]]; then
    generator=$annotationGenerate/annotationGenerateZM.R
    echo "Model organism detected: Zea mays"
elif [[ $generator == "OSJ" ]]; then
    generator=$annotationGenerate/annotationGenerateOSJ.R
    echo "Model organism detected: Oryza sativa jadica"
elif [[ $generator == "OSIR64" ]]; then
    generator=$annotationGenerate/annotationGenerateOSIR64.R
    echo "Model organism detected: Oryza sativa IR64"
else
    echo "##################################################"
    echo "model organism code not recognized, please check your input"
    echo "HAMRLINC will proceed with limited functionalities"
    echo "##################################################"
fi

# designate log file, if exist, clear, and have all stdout written
logstart=$(date "+%Y.%m.%d-%H.%M.%S")
logfile=$out/Log_$logstart.log
exec > >(tee -a "$logfile") 2>&1
#below captures only echo...?
#exec 2>&1 1>>$logfile 3>&1

######################################################### Subprogram Definition #########################################
checkpoint () {
    echo "Checkpoint reached: $1"
    echo "$1" > "$out"/checkpoint.txt
}

fqgrab () {

  echo "begin downloading $line..." 

  fasterq-dump "$line" -O "$dumpout"/raw --verbose

  # automatically detects the suffix
  echo "$dumpout"/raw/"$line"
  if [[ -f $dumpout/raw/$line"_1.fastq" ]]; then
    suf="fastq"
    PE=true
    echo "$line is a paired-end file ending in .fastq"
  elif [[ -f $dumpout/raw/$line"_1.fq" ]]; then
    suf="fq"
    PE=true
    echo "$line is a paired-end file ending in .fq"
  elif [[ -f $dumpout/raw/$line".fastq" ]]; then
    suf="fastq"
    PE=false
    echo "$line is a single-end file ending in .fastq"
  elif [[ -f $dumpout/raw/$line".fq" ]]; then
    suf="fq"
    PE=false
    echo "$line is a single-end file ending in .fq"
  else
    echo "suffix not recognized, please check your datasets"
    exit 1
  fi

  if [[ "$PE" = false ]]; then  
    echo "[$line] performing fastqc on raw file..."
    fastqc "$dumpout"/raw/"$line"."$suf" -o "$dumpout"/fastqc_results &

    echo "[$line] trimming..."
    trim_galore -o "$dumpout"/trimmed "$dumpout"/raw/"$line"."$suf"

    echo "[$line] trimming complete, performing fastqc..."
    fastqc "$dumpout"/trimmed/"$line""_trimmed.fq" -o "$dumpout"/fastqc_results

  else 
    echo "[$line] performing fastqc on raw file..."
    fastqc "$dumpout"/raw/"$line""_1.$suf" -o "$dumpout"/fastqc_results &
    fastqc "$dumpout"/raw/"$line""_2.$suf" -o "$dumpout"/fastqc_results &

    echo "[$line] trimming..."
    trim_galore -o "$dumpout"/trimmed "$dumpout"/raw/"$line""_1.$suf"
    trim_galore -o "$dumpout"/trimmed "$dumpout"/raw/"$line""_2.$suf"

    echo "[$line] trimming complete, performing fastqc..."
    fastqc "$dumpout"/trimmed/"$line""_1_trimmed.fq" -o "$dumpout"/fastqc_results
    fastqc "$dumpout"/trimmed/"$line""_2_trimmed.fq" -o "$dumpout"/fastqc_results
  fi
  echo "[$(date '+%d/%m/%Y %H:%M:%S')] finished processing $line"
  echo ""
}

fqgrab2 () {
    sname=$(basename "$fq")
    tt=${sname%.*} 
    echo "[$sname] performing fastqc on raw file..."
    fastqc "$fq" -o "$dumpout"/fastqc_results &

    echo "[$sname] trimming..."
    trim_galore -o "$dumpout"/trimmed "$fq"

    echo "[$sname] trimming complete, performing fastqc..."
    fastqc "$dumpout"/trimmed/"$tt""_trimmed.fq" -o "$dumpout"/fastqc_results
}

fastq2hamr () {
    smpext=$(basename "$smp")
    smpdir=$(dirname "$smp")
    smpkey="${smpext%.*}"
    smpname=""
    original_ext="${smpext##*.}"

    if [[ $smpkey == *_1* ]]; then
        smpkey="${smpkey%_1*}"
        smp1="$smpdir/${smpkey}_1_trimmed.$original_ext"
        smp2="$smpdir/${smpkey}_2_trimmed.$original_ext"
        # Paired end recognized
        det=0
        echo "$smpext is a part of a paired-end sequencing file"
    elif [[ $smpkey == *_2* ]]; then
        # If _2 is in the filename, this file was processed along with its corresponding _1 so we skip
        echo "$smpext has already been processed with its _1 counter part. Skipped."
        echo ""
        exit 1
    else
        det=1
        echo "$smpext is a single-end sequencing file"
        echo ""
    fi

    # Read the CSV file into a DataFrame
    mapfile -t names < <(awk -F, '{ print $1 }' "$csv")
    mapfile -t smpf < <(awk -F, '{ print $2 }' "$csv")

    # Create a dictionary from the DataFrame
    declare -A dictionary
    for ((i=0; i<${#names[@]}; i++)); do
        dictionary[${names[i]}]=${smpf[i]}
    done

    if [[ $smpkey == *_trimmed* ]]; then
        smpkey="${smpkey%_trimmed*}"
    fi

    # Retrieve the translated value
    if [[ ${dictionary[$smpkey]+_} ]]; then
        smpname="${dictionary[$smpkey]}"
        smpname="${smpname//$'\r'}"
        echo "[$smpkey] Sample group name found: $smpname"
    else
        echo "[$smpkey] Could not locate sample group name, exiting..."
        exit 1
    fi

    # Reassign / declare pipeline file directory
    if [ ! -d "$out/pipeline/$smpkey""_temp" ]; then
        mkdir "$out/pipeline/$smpkey""_temp"
        echo "[$smpkey] created path: $out/pipeline/$smpkey""_temp"
    fi

    smpout=$out/pipeline/$smpkey"_temp"
    echo "[$smpkey] You can find all the intermediate files for $smpkey at $smpout" 


    # Reassign hamr output directory
    if [ ! -d "$out/hamr_out" ]; then
        mkdir "$out"/hamr_out
        echo "created path: $out/hamr_out"
    fi

    hamrout=$out/hamr_out
    echo "[$smpkey] You can find the HAMR output file for $smpkey at $hamrout/$smpname.mod.txt" 


    echo "$(date '+%d/%m/%Y %H:%M:%S') [$smpkey] Begin HAMR pipeline"
    cd "$smpout" || exit
    # maps the trimmed reads to provided annotated genome, can take ~1.5hr

    if [[ "$tophat" = false ]]; then  
        echo "Using STAR for mapping..."
        if [ "$det" -eq 1 ]; then
            echo "[$smpkey] Performing STAR with a single-end file."
            STAR \
            --runThreadN "$threads" \
            --genomeDir "$out"/ref/ \
            --readFilesIn "$smp" \
            --sjdbOverhang $overhang \
            --sjdbGTFfile "$annotation" \
            --sjdbGTFtagExonParentTranscript Parent \
            --outFilterMultimapNmax 10 \
            --outFilterMismatchNmax $mismatch \
            --outSAMtype BAM SortedByCoordinate
        else
            echo "[$smpkey] Performing STAR with a paired-end file."
            STAR \
            --runThreadN "$threads" \
            --genomeDir "$out"/ref/ \
            --readFilesIn "$smp1" "$smp2" \
            --sjdbOverhang $overhang \
            --sjdbGTFfile "$annotation" \
            --sjdbGTFtagExonParentTranscript Parent \
            --outFilterMultimapNmax 10 \
            --outFilterMismatchNmax $mismatch \
            --outSAMtype BAM SortedByCoordinate
        fi

    else
        echo "Using TopHat2 for mapping..."
        # set read distabce based on mistmatch num
        red=8
        if [[ $mismatch -gt 8 ]]; then red=$((mismatch +1)); fi

        if [ "$det" -eq 1 ]; then
            echo "[$smpkey] Performing TopHat2 with a single-end file."
            tophat2 \
                --library-type "$tophatlib" \
                --read-mismatches $mismatch \
                --read-edit-dist $red \
                --max-multihits 10 \
                --b2-very-sensitive \
                --transcriptome-max-hits 10 \
                --no-coverage-search \
                -G "$annotation" \
                -p "$threads" \
                "$out"/btref \
                "$smp"
        else
        echo "[$smpkey] Performing TopHat2 with a paired-end file."
            tophat2 \
                --library-type "$tophatlib" \
                --read-mismatches $mismatch \
                --read-edit-dist $red \
                --max-multihits 10 \
                --b2-very-sensitive \
                --transcriptome-max-hits 10 \
                --no-coverage-search \
                -G "$annotation" \
                -p "$threads" \
                "$out"/btref \
                "$smp1" "$smp2"
        fi
    fi
    cd || exit

    wait

    #sorts the accepted hits
    echo "[$smpkey] sorting..."
    # handles tophat or star output
    if [[ "$tophat" = false ]]; then
        samtools sort \
        -n "$smpout"/Aligned.sortedByCoord.out.bam \
        -o "$smpout"/sort_accepted.bam
    else
        samtools sort \
        -n "$smpout"/accepted_hits.bam \
        -o "$smpout"/sort_accepted.bam
    fi
    echo "[$smpkey] finished sorting"
    echo ""

    wait

    #filter the accepted hits by uniqueness
    echo "[$smpkey] filter unique..."
    samtools view \
        -h "$smpout"/sort_accepted.bam \
        | perl "$filter" 1 \
        | samtools view -bS - \
        | samtools sort \
        -o "$smpout"/unique.bam
    echo "[$smpkey] finished filtering"
    echo ""

    wait

    # translates string library prep strandedness into feature count required number
    if [[ "$tophatlib" = fr-firststrand ]]; then
        fclib=2
    elif [[ "$tophatlib" = fr-secondstrand ]]; then
        fclib=1
    else 
        fclib=0
    fi

    ###############################################
    ########evolinc_i logic here (left arm)###########
    ############################################### 
    # if user didn't suppress evolinc_i, start the pipeline, note the constitutive featurecount after evolinc
    if [[ "$evolinc_i" = true ]]; then
        echo "################################################################"
        echo "############## Entering lincRNA abundance quantification pipeline ##############"
        echo "################################################################"
        date '+%d/%m/%Y %H:%M:%S'
        if [ ! -d "$out/evolinc_out" ]; then mkdir "$out/evolinc_out"; fi
        # run stringtie accordingly, note PE and SE here is taken care of
        # output is unnamed and stored in each fastq folder
        echo "[$smpkey] producing transcript assembly using stringtie..."
        if [[ "$tophatlib" = fr-firststrand ]]; then
            echo "[$smpkey] running stringtie with --rf"
            stringtie \
                "$smpout"/unique.bam \
                -o "$smpout"/transcriptAssembly.gtf \
                -G "$annotation" \
                -p "$threads" \
                --rf
        elif [[ "$tophatlib" = fr-secondstrand ]]; then
            echo "[$smpkey] running stringtie with --fr"
            stringtie \
                "$smpout"/unique.bam \
                -o "$smpout"/transcriptAssembly.gtf \
                -G "$annotation" \
                -p "$threads" \
                --fr
        else
            echo "[$smpkey] running stringtie assuming an unstranded library"
            stringtie \
                "$smpout"/unique.bam \
                -o "$smpout/"transcriptAssembly.gtf \
                -G "$annotation" \
                -p "$threads"
        fi

        # next run cuff compare
        echo "[$smpkey] merging assemblies using cuffcompare..."
        cuffcompare \
            "$smpout"/transcriptAssembly.gtf \
            -r "$annotation" \
            -s "$genome" \
            -T \
            -o "$smpout"/cuffed

        # run evolinc
        echo "[$smpkey] annotating lincRNA using Evolinc-i..."
        if [ "$evolinc_i_option" == "M" ]; then
            echo "[$smpkey] M option identified for evolinc"
            evolinc-part-I.sh \
                -c "$smpout"/cuffed.combined.gtf \
                -g "$genome" \
                -u "$annotation" \
                -r "$annotation" \
                -n "$threads" \
                -o "$smpout"/"$smpname"_lincRNA
        elif [ "$evolinc_i_option" == "MO" ]; then
            echo "[$smpkey] MO option identified for evolinc"
            evolinc-part-I.sh \
                -c "$smpout"/cuffed.combined.gtf \
                -g "$genome" \
                -u "$annotation" \
                -r "$annotation" \
                -n "$threads" \
                -o "$smpout"/"$smpname"_lincRNA \
                -b "$blast_file" \
                -t "$cage_file" \
                -x "$known_linc"
        fi

        cd $smpout
        # house keeping for evolinc
        rm *.loci *.stats *.tracking
        mv *_lincRNA* "$out/evolinc_out"
        cd

        # run constitutive feature count within evolinc_i (left arm) if the user didn't suppress feacturecount
        if [[ "$featurecount" = true ]]; then
            echo "[$(date '+%d/%m/%Y %H:%M:%S')$smpkey] quantifying lincRNA-based transcript abundance using featurecounts..."
            if [ ! -d "$out/featurecount_out" ]; then mkdir "$out/featurecount_out"; fi
            if [ "$det" -eq 1 ]; then
                echo "[$smpkey] running featurecount with $fclib as the -s argument"
                featureCounts \
                    -T "$threads" \
                    -s $fclib \
                    -a "$smpout"/lincRNA/lincRNA.updated.gtf \
                    -o "$smpout"/"$smpname"_lincRNA_featurecount.txt \
                    "$smpout"/unique.bam
            else
                featureCounts \
                    -T "$threads" \
                    -a "$smpout"/lincRNA/lincRNA.updated.gtf \
                    -o "$smpout"/"$smpname"_lincRNA_featurecount.txt \
                    "$smpout"/unique.bam
            fi

            # housekeeping for feature counts
            cd "$smpout"
            mv *_featurecount.txt* "$out/featurecount_out"
            cd
        fi 
        echo "################################################################"
        echo "############## lincRNA abundance quantification pipeline completed ##############"
        echo "################################################################"
    fi

    wait

    ###############################################
    ########regular feature count logic here (right arm)###########
    ############################################### 
    # first create gtf file from gff3 file
    gffread \
        "$annotation" \
        -T \
        -o "$out"/ref/temp.gtf

    # run feature count for normal alignment transcript quantification if user didn't suppress
    if [[ "$featurecount" = true ]]; then
        echo "[$(date '+%d/%m/%Y %H:%M:%S')$smpkey] quantifying regular transcript abundance using featurecounts..."
        if [ ! -d "$out/featurecount_out" ]; then mkdir "$out/featurecount_out"; fi
        if [ "$det" -eq 1 ]; then
            echo "[$smpkey] running featurecount with $fclib as the -s argument"
            featureCounts \
                -T "$threads" \
                -s $fclib \
                -a "$out"/ref/temp.gtf \
                -o "$smpout"/"$smpname"_alignment_featurecount.txt \
                "$smpout"/unique.bam
        else
            featureCounts \
                -T "$threads" \
                -a "$out"/ref/temp.gtf \
                -o "$smpout"/"$smpname"_alignment_featurecount.txt \
                "$smpout"/unique.bam
        fi

        # housekeeping for regular abundance
        cd "$smpout"
        mv *_featurecount.txt* "$out/featurecount_out"
        cd
    fi
    wait

    ###############################################
    ########original continuation of fastq2hamr here###########
    ############################################### 
    # run below only if hamrbox is true
    if [[ "$hamrbox" = false ]]; then
        echo "[$(date '+%d/%m/%Y %H:%M:%S')] hamrbox functionality suppressed, $smpkey analysis completed."
    else
        #adds read groups using picard, note the RG arguments are disregarded here
        echo "[$smpkey] adding/replacing read groups..."
        gatk AddOrReplaceReadGroups \
            I="$smpout"/unique.bam \
            O="$smpout"/unique_RG.bam \
            RGID=1 \
            RGLB=xxx \
            RGPL=illumina_100se \
            RGPU=HWI-ST1395:97:d29b4acxx:8 \
            RGSM=sample
        echo "[$smpkey] finished adding/replacing read groups"
        echo ""

        wait

        #reorder the reads using picard
        echo "[$smpkey] reordering..."
        echo "$genome"
        gatk --java-options "-Xmx2g -Djava.io.tmpdir=$smpout/tmp" ReorderSam \
            I="$smpout"/unique_RG.bam \
            O="$smpout"/unique_RG_ordered.bam \
            R="$genome" \
            CREATE_INDEX=TRUE \
            SEQUENCE_DICTIONARY="$dict" \
            TMP_DIR="$smpout"/tmp
        echo "[$smpkey] finished reordering"
        echo ""

        wait

        #splitting and cigarring the reads, using genome analysis tool kit
        #note can alter arguments to allow cigar reads 
        echo "[$smpkey] getting split and cigar reads..."
        gatk --java-options "-Xmx2g -Djava.io.tmpdir=$smpout/tmp" SplitNCigarReads \
            -R "$genome" \
            -I "$smpout"/unique_RG_ordered.bam \
            -O "$smpout"/unique_RG_ordered_splitN.bam
            # -U ALLOW_N_CIGAR_READS
        echo "[$smpkey] finished splitting N cigarring"
        echo ""

        wait

        #final resorting using picard
        echo "[$smpkey] resorting..."
        gatk --java-options "-Xmx2g -Djava.io.tmpdir=$smpout/tmp" SortSam \
            I="$smpout"/unique_RG_ordered_splitN.bam \
            O="$smpout"/unique_RG_ordered_splitN.resort.bam \
            SORT_ORDER=coordinate
        echo "[$smpkey] finished resorting"
        echo ""

        wait

        #hamr step, can take ~1hr
        echo "[$smpkey] hamr..."
        #hamr_path=$(which hamr.py) 
        python /HAMR/hamr.py \
            -fe "$smpout"/unique_RG_ordered_splitN.resort.bam "$genome" "$model" "$smpout" $smpname $quality $coverage $err H4 $pvalue $fdr .05
        wait

        if [ ! -e "$smpout/${smpname}.mods.txt" ]; then 
            cd "$hamrout" || exit
            printf '%s \n' "$smpname" >> zero_mod.txt
            cd || exit
        else
        # HAMR needs separate folders to store temp for each sample, so we move at the end
            cp "$smpout"/"${smpname}".mods.txt "$hamrout"
        fi

        # Move the unique_RG_ordered.bam and unique_RG_ordered.bai to a folder for read depth analysis
        cp "$smpout"/unique_RG_ordered.bam "$out"/pipeline/depth/"$smpname".bam
        cp "$smpout"/unique_RG_ordered.bai "$out"/pipeline/depth/"$smpname".bai

        # delete more intermediate files?
        
    fi
}

consensusOverlap () {
    IFS="/" read -ra sections <<< "$smp"
    temp="${sections[-1]}"

    IFS="." read -ra templ <<< "$temp"
    smpname="${templ[0]}"

    echo "consensus file prefix: $smpname"
    echo ""

    count=$(ls -1 "$genomedir"/*_CDS.bed 2>/dev/null | wc -l)
    if [ "$count" != 0 ]; then 
        cds=$(find "$genomedir" -maxdepth 1 -name "*_CDS.bed")
        #overlap with cds
        intersectBed \
            -a "$cds" \
            -b "$smp" \
            -wa -wb \
            > "$out"/lap/"$smpname""_CDS".bed
        echo "finished finding overlap with CDS library"
    fi

    count=$(ls -1 "$genomedir"/*_fiveUTR.bed 2>/dev/null | wc -l)
    if [ "$count" != 0 ]; then 
        fiveutr=$(find "$genomedir" -maxdepth 1 -name "*_fiveUTR.bed")
        #overlap with 5utr
        intersectBed \
            -a "$fiveutr" \
            -b "$smp" \
            -wa -wb \
            > "$out"/lap/"$smpname""_fiveUTR".bed
        echo "finished finding overlap with 5UTR library"
    fi

    count=$(ls -1 "$genomedir"/*_threeUTR.bed 2>/dev/null | wc -l)
    if [ "$count" != 0  ]; then 
        threeutr=$(find "$genomedir" -maxdepth 1 -name "*_threeUTR.bed")
        #overlap with 3utr
        intersectBed \
            -a "$threeutr" \
            -b "$smp" \
            -wa -wb \
            > "$out"/lap/"$smpname""_threeUTR".bed
        echo "finished finding overlap with 3UTR library"
    fi

    count=$(ls -1 "$genomedir"/*_gene.bed 2>/dev/null | wc -l)
    if [ "$count" != 0 ]; then 
        gene=$(find "$genomedir" -maxdepth 1 -name "*_gene.bed")
        #overlap with gene
        intersectBed \
            -a "$gene" \
            -b "$smp" \
            -wa -wb \
            > "$out"/lap/"$smpname""_gene".bed
        echo "finished finding overlap with gene library"
    fi

    count=$(ls -1 "$genomedir"/*_primarymRNA.bed 2>/dev/null | wc -l)
    if [ "$count" != 0 ]; then 
        mrna=$(find "$genomedir" -maxdepth 1 -name "*_primarymRNA.bed")
        #overlap with mrna
        intersectBed \
            -a "$mrna" \
            -b "$smp" \
            -wa -wb \
            > "$out"/lap/"$smpname""_primarymRNA".bed
        echo "finished finding overlap with primary mRNA library"
    fi

    count=$(ls -1 "$genomedir"/*_exon.bed 2>/dev/null | wc -l)
    if [ "$count" != 0 ]; then 
        exon=$(find "$genomedir" -maxdepth 1 -name "*_exon.bed")
        #overlap with exon
        intersectBed \
            -a "$exon" \
            -b "$smp" \
            -wa -wb \
            > "$out"/lap/"$smpname""_exon".bed
        echo "finished finding overlap with exon library"
    fi

    count=$(ls -1 "$genomedir"/*_ncRNA.bed 2>/dev/null | wc -l)
    if [ "$count" != 0 ]; then 
        nc=$(find "$genomedir" -maxdepth 1 -name "*_ncRNA.bed")
        #overlap with nc rna
        intersectBed \
            -a "$nc" \
            -b "$smp" \
            -wa -wb \
            > "$out"/lap/"$smpname""_ncRNA".bed
        echo "finished finding overlap with ncRNA library"
    fi
}

fqgrabhouse () {
    ##########fqgrab housekeeping begins#########
    if [ ! -d "$out" ]; then mkdir "$out"; echo "created path: $out"; fi

    if [ ! -d "$out/datasets" ]; then mkdir "$out"/datasets; echo "created path: $out/datasets"; fi

    # first see what input is provided
    if [[ $acc == *.txt ]]; then
        echo "SRR accession list provided, using fasterq-dump for .fastq acquisition..."

        # Create directory to store original fastq files
        if [ ! -d "$out/datasets/raw" ]; then mkdir "$out"/datasets/raw; fi
        echo "You can find your original fastq files at $out/datasets/raw" 
        mode=1

    elif [[ -d $acc ]]; then
        echo "Directory $acc is found, assuming raw fastq files are provided..."
        mode=2
    else
        echo "Error recognizing input source, exiting..."
        exit 1
    fi

    # relocate user-provided inputs
    if [ ! -d "$out/ref" ]; then 
        mkdir "$out/ref"
        echo "created path: $out/ref"
        #cp $genome "$out/ref"
        #genome="$out/ref/$(basename $genome)"
        #cp $annotation "$out/ref"
        #annotation="$out/ref/$(basename $annotation)"
    fi

    if [ ! -d "$out/fileprep" ]; then 
        mkdir "$out/fileprep"
        echo "created path: $out/fileprep"
        #cp $acc "$out/fileprep"
        #cp $csv "$out/fileprep"
    fi

    # Create directory to store trimmed fastq files
    if [ ! -d "$out/datasets/trimmed" ]; then mkdir "$out"/datasets/trimmed; fi
    echo "You can find your trimmed fastq files at $out/datasets/trimmed"

    # Create directory to store fastqc results
    if [ ! -d "$out/datasets/fastqc_results" ]; then mkdir "$out"/datasets/fastqc_results; fi
    echo "You can find all the fastqc test results at $out/datasets/fastqc_results"

    # Run a series of command checks to ensure the entire script can run smoothly
    if ! command -v fasterq-dump > /dev/null; then
        echo "Failed to call fasterq-dump command. Please check your installation."
        exit 1
    fi

    if ! command -v fastqc > /dev/null; then
        echo "Failed to call fastqc command. Please check your installation."
        exit 1
    fi

    if ! command -v trim_galore > /dev/null; then
        echo "Failed to call trim_galore command. Please check your installation."
        exit 1
    fi

    if ! command -v gatk > /dev/null; then
        echo "Failed to call gatk command. Please check your installation."
        exit 1
    fi
    ##########fqgrab housekeeping ends#########
}

fastq2hamrhouse () {
    ############fastq2hamr housekeeping begins##############
    # Checks if the files were trimmed or cleaned, and if so, take those files for downstream
    hamrin=""
    suf=""
    # If trimmed folder present, then user specified trimming, we take trimmed files with .fq
    if [ -d "$dumpout/trimmed" ]; then 
        hamrin=$dumpout/trimmed
        suf="fq"
    else
        echo "failed to locate trimmed fastq files"
        exit 1
    fi

    # Creating some folders
    if [ ! -d "$out/pipeline" ]; then mkdir "$out"/pipeline; echo "created path: $out/pipeline"; fi

    if [ ! -d "$out/hamr_out" ]; then mkdir "$out"/hamr_out; echo "created path: $out/hamr_out"; fi

    # Check if zero_mod is present already, if not then create one
    if [ ! -e "$out/hamr_out/zero_mod.txt" ]; then
        cd "$out/hamr_out" || exit
        echo "Below samples have 0 HAMR predicted mods:" > zero_mod.txt
        cd || exit
    fi


    # create dict file using fasta genome file
    count=$(ls -1 "$genomedir"/*.dict 2>/dev/null | wc -l)
    if [ "$count" == 0 ]; then 
    gatk CreateSequenceDictionary \
        R="$genome"
    fi
    dict=$(find "$genomedir" -maxdepth 1 -name "*.dict")

    # create fai index file using fasta genome
    count=$(ls -1 "$genomedir"/*.fai 2>/dev/null | wc -l)
    if [ "$count" == 0 ]; then 
    samtools faidx "$genome"
    fi

    # Check which mapping software, and check for index
    if [[ "$tophat" = false ]]; then  
    # Check if indexed files already present for STAR
        if [ -e "$out/ref/SAindex" ]; then
            echo "STAR Genome Directory with indexed genome detected, skipping STAR indexing"
        else
            # Now, do the indexing step
            # Define the SA index number argument
            log_result=$(echo "scale=2; l($genomelength)/l(2)/2 - 1" | bc -l)
            sain=$(echo "scale=0; if ($log_result < 14) $log_result else 14" | bc)
            echo "Creating STAR genome index..."
            # Create genome index 
            STAR \
                --runThreadN $threads \
                --runMode genomeGenerate \
                --genomeDir "$out"/ref \
                --genomeFastaFiles "$genome" \
                --sjdbGTFfile "$annotation" \
                --sjdbGTFtagExonParentTranscript Parent \
                --sjdbOverhang $overhang \
                --genomeSAindexNbases $sain
        fi
    else
        # Check if bowtie index directory is already present
        if [ -e "$out/btref" ]; then
            echo "bowtie indexed directory detected, skipping generating bowtie index"
        else
        # If not, first check if ref folder is present, if not then make
            if [ ! -d "$out/btref" ]; then mkdir "$out/btref"; echo "created path: $out/btref"; fi
            echo "Creating Bowtie references..."
            bowtie2-build "$genome" "$out"/btref
        fi
    fi

    # Run a series of command checks to ensure fastq2hamr can run smoothly
    if ! command -v mapfile > /dev/null; then
        echo "Failed to call mapfile command. Please check your installation."
        exit 1
    fi

    if ! command -v STAR > /dev/null; then
        echo "Failed to call STAR command. Please check your installation."
        exit 1
    fi

    if ! command -v samtools > /dev/null; then
        echo "Failed to call samtools command. Please check your installation."
        exit 1
    fi

    if ! command -v stringtie > /dev/null; then
        echo "Failed to call stringtie command. Please check your installation."
        exit 1
    fi

    if ! command -v cuffcompare > /dev/null; then
        echo "Failed to call cuffcompare command. Please check your installation."
        exit 1
    fi

    if ! command -v featureCounts > /dev/null; then
        echo "Failed to call featureCounts command. Please check your installation."
        exit 1
    fi

    if ! command -v gatk > /dev/null; then
        echo "Failed to call gatk command. Please check your installation."
        exit 1
    fi

    if ! command -v python > /dev/null; then
        echo "Failed to call python command. Please check your installation."
        exit 1
    fi

    # Creates a folder for depth analysis
    if [ ! -d "$out/pipeline/depth" ]; then mkdir "$out"/pipeline/depth; echo "created path: $out/pipeline/depth"; fi
    #############fastq2hamr housekeeping ends#############
}

######################################################### Main Program Begins #########################################

echo ""
echo "##################################### Begin HAMRLINC #################################"
echo ""

# Check if the required arguments are provided
if [ -z "$out" ]; then 
    echo "output directory not detected, exiting..."
    exit 1
elif [ -z "$acc" ]; then
    echo "input SRR or fastq files not detected, exiting..."
    exit 1
elif [ -z "$csv" ]; then
    echo "filename dictionary csv not detected, exiting..."
    exit 1
elif [ -z "$genome" ]; then
    echo "model organism genmome fasta not detected, exiting..."
    exit 1
elif [ -z "$annotation" ]; then
    echo "model organism genmome annotation gff3 not detected, exiting..."
    exit 1
elif [ -z "$length" ]; then
    echo "read length not detected, exiting..."
    exit 1
elif [ -z "$genomelength" ]; then
    echo "genome size not detected, exiting..."
    exit 1
else
    echo "all required arguments provided, proceding..."
fi

# check that the user didn't suppress all three programs -- if so, there's no need to run anything
if [ $evolinc_i = false ] && [ $featurecount = false ] && [ $hamrbox = false ]; then
    echo "User has not activated any functionalities. Exiting..."
    exit 0
fi

# check whether checkpoint.txt is present 
if [ -e "$out"/checkpoint.txt ]; then
    last_checkpoint=$(cat "$out"/checkpoint.txt)
    echo "Resuming from checkpoint: $last_checkpoint"
else
    last_checkpoint="start"
fi

# run fqgrab when checkpoint agrees so
if [ "$last_checkpoint" = "start" ] || [ "$last_checkpoint" = "" ]; then
    fqgrabhouse
    ##########fqgrab main begins#########
    if [[ $mode -eq 1 ]]; then
        # Grabs the fastq files from acc list provided into the dir ~/datasets
        i=0
        while IFS= read -r line
        do ((i=i%threads)); ((i++==0)) && wait
        fqgrab &
        done < "$acc"

    elif [[ $mode -eq 2 ]]; then
        i=0
        for fq in "$acc"/*; do
        ((i=i%threads)); ((i++==0)) && wait
        fqgrab2 &
        done
    fi
    wait
    ##################fqgrab main ends#################

    echo ""
    echo "################ Finished downloading and processing all fastq files. Entering pipeline for HAMR analysis. ######################"
    date '+%d/%m/%Y %H:%M:%S'
    echo ""
    # obtained all processed fastq files, record down checkpoint
    last_checkpoint="checkpoint1"
    checkpoint $last_checkpoint
fi

# run fastq2hamr when checkpoint agrees
if [ "$last_checkpoint" = "checkpoint1" ]; then 
    fastq2hamrhouse
    #############fastq2hamr main begins###############
    # Pipes each fastq down the hamr pipeline, and stores out put in ~/hamr_out
    # Note there's also a hamr_out in ~/pipeline/SRRNUMBER_temp/, but that one's for temp files
    
    #mkdir trimmed_temp && mv "$hamrin"/*."$suf" trimmed_temp && chmod -R 777 trimmed_temp
    #cd trimmed_temp
    #current_dir=$(pwd)
    #cd ..

    i=0
    ttop=$((threads/2))
    for smp in "$hamrin"/*."$suf"
    do ((i=i%ttop)); ((i++==0)) && wait
    fastq2hamr #&& mv "$trimmed_temp"/*".$suf" "$hamrin"/ && rm -r trimmed_temp &
    done

    if [[ "$hamrbox" = false ]]; then
        exit 0
    fi

    wait

    # Check whether any hamr.mod.text is present, if not, halt the program here
    if [ -z "$(ls -A "$out"/hamr_out)" ]; then
    echo "No HAMR predicted mod found for any sequencing data in this project, please see log for verification"
    exit 1
    fi

    # If program didn't exit, at least 1 mod file, move zero mod record outside so it doesn't get read as a modtbl next
    mv "$out"/hamr_out/zero_mod.txt "$out"

    echo ""
    echo "################ Finished HAMR analysis. Producing consensus mod table and depth analysis. ######################"
    echo "$(date '+%d/%m/%Y %H:%M:%S')"
    echo ""

    #############fastq2hamr main ends###############

    # obtained all HAMR txts, record down checkpoint
    last_checkpoint="checkpoint2"
    checkpoint $last_checkpoint
fi

# run consensus when checkpoint agrees
if [ "$last_checkpoint" = "checkpoint2" ]; then 
    ##############consensus finding begins##############
    # Produce consensus bam files based on filename (per extracted from name.csv) and store in ~/consensus
    if [ ! -d "$out/consensus" ]; then mkdir "$out"/consensus; echo "created path: $out/consensus"; fi

    # Run a series of command checks to ensure findConsensus can run smoothly
    if ! command -v Rscript > /dev/null; then
        echo "Failed to call Rscript command. Please check your installation."
        exit 1
    fi

    echo "Producing consensus file across biological replicates..."
    # Find consensus accross all reps of a given sample group
    Rscript "$scripts"/findConsensus.R \
        "$out"/hamr_out \
        "$out"/consensus
    wait
    echo "done"

    # The case where no consensus file is found, prevents *.bed from being created
    if [ -z "$(ls -A "$out"/consensus)" ]; then
    echo "No consensus mods found within any sequencing group. Please see check individual rep for analysis. "
    exit 1
    fi

    # Add depth columns with info from each rep alignment, mutate in place
    for f in "$out"/consensus/*.bed
    do
    t=$(basename "$f")
    d=$(dirname "$f")
    n=${t%.*}
    echo "starting depth analysis on $n"
    for ff in "$out"/pipeline/depth/*.bam
        do
            if echo "$ff" | grep -q "$n"
            then
                tt=$(basename "$ff")
                nn=${tt%.*}
                echo "[$n] extracting depth information from $nn"
                for i in $(seq 1 $(wc -l < "$f"))
                do
                    chr=$(sed "${i}q;d" "$f" | sed 's/\t/\n/g' | sed '1q;d')
                    pos=$(sed "${i}q;d" "$f" | sed 's/\t/\n/g' | sed '2q;d')
                    dph=$(samtools coverage \
                        -r "$chr":"$pos"-"$pos" \
                        "$ff" \
                        | awk 'NR==2' | awk -F'\t' '{print $7}')
                    awk -v "i=$i" 'NR==i {print $0"\t"var; next} 1' var="$dph" "$f" > "$d"/"${nn}"_new.bed && mv "$d"/"${nn}"_new.bed "$f" 
                done
                echo "[$n] finished $nn"
                fi
            done &
        done
    wait

    for f in "$out"/consensus/*.bed
    do
    if [ -s "$f" ]; then
    # The file is not-empty.
        t=$(basename "$f")
        n=${t%.*}
        echo "computing depth across reps for $n"
        Rscript "$scripts"/depth_helper_average.R "$f"
    fi
    done

    wait

    #############consensus finding ends###############

    # obtained all consensus HAMR mods with depth, record down checkpoint
    last_checkpoint="checkpoint3"
    checkpoint $last_checkpoint
fi

# run overlap when checkpoint agrees
if [ "$last_checkpoint" = "checkpoint3" ]; then 
    ##############overlapping begins##############
    # Produce overlap bam files with the provided annotation library folders and store in ~/lap
    if [ ! -d "$out/lap" ]; then mkdir "$out"/lap; echo "created path: $out/lap"; fi

    # Run a series of command checks to ensure consensusOverlap can run smoothly
    if ! command -v intersectBed > /dev/null; then
        echo "Failed to call intersectBed command. Please check your installation."
        exit 1
    fi

    # checks if genomedir is populated with generated annotation files, if not, hamrbox can't run anymore, exit
    count=$(ls -1 "$genomedir"/*.bed 2>/dev/null | wc -l)
    if [ "$count" == 0 ]; then 
        if [[ -n "$generator" ]]; then
            echo "generating annotations for overlap..."
            Rscript "$generator" "$annotation"
        else
            echo "#########NOTICE###########"
            echo "##########No annotation generator or annotation files found, please check your supplied arguments##########"
            echo "##########As a result, HAMRLINC will stop here. Please provide the above files in the next run############"
            exit 1
        fi
    else 
        echo "generated annotation detected, proceeding to overlapping"
    fi

    # Overlap with provided libraries for each sample group
    for smp in "$out"/consensus/*
    do consensusOverlap
    done

    if [ -z "$(ls -A "$out"/lap)" ]; then
    echo "No overlapped mods found within any sequencing group. Please see check individual rep for analysis. "
    exit 1
    fi

    #############overlapping ends###############

    # obtained all overlapped HAMR mods, record down checkpoint
    last_checkpoint="checkpoint4"
    checkpoint $last_checkpoint
fi

# run R analysis when checkpoint agrees
if [ "$last_checkpoint" = "checkpoint4" ]; then 
    ##############R analysis begins##############
    echo ""
    echo "###############SMACK portion completed, entering EXTRACT################"
    date '+%d/%m/%Y %H:%M:%S'
    echo ""
    #######################################begins EXTRACT######################################
    dir=$out

    echo "generating long modification table..."
    # collapse all overlapped data into longdf
    Rscript "$scripts"/allLapPrep.R \
        "$dir"/lap \
        "$dir"
    echo "done"
    echo ""

    echo "plotting modification abundance..."
    # overview of modification proportion
    Rscript "$scripts"/abundByLap.R \
        "$dir"/mod_long.csv \
        "$genomedir" \
        "$dir"
    echo "done"
    echo ""

    echo "performing modification cluster analysis..."
    # analyze hamr-mediated/true clustering across project
    Rscript "$scripts"/clusterAnalysis.R \
        "$dir"/mod_long.csv \
        "$dir"
    echo "done"
    echo ""

    # if [ ! -z "${4+x}" ]; then
    #     echo "known modification landscape provided, performing relative positional analysis to known mod..."
    #     # The csv (in modtbl format) of the known mod you want analyzed in distToKnownMod
    #     antcsv=$4
    #     # analyze hamr-mediated/true clustering across project
    #     Rscript $scripts/distToKnownMod.R \
    #         $dir/mod_long.csv \
    #         $antcsv
    #     echo "done"
    #     echo ""
    # else 
    #     echo "known modification file not detected, skipping relative positional analysis"
    #     echo ""
    # fi

    echo "classifying modified RNA subtype..."
    # looking at RNA subtype for mods
    Rscript "$scripts"/RNAtype.R \
        "$dir"/mod_long.csv
    echo "done"
    echo ""

    if [ ! -d "$dir/go" ]; then mkdir "$dir"/go; echo "created path: $dir/go"; fi

    if [ ! -d "$dir/go/genelists" ]; then mkdir "$dir"/go/genelists; echo "created path: $dir/go/genelists"; fi

    if [ ! -d "$dir/go/pantherout" ]; then mkdir "$dir"/go/pantherout; echo "created path: $dir/go/pantherout"; fi

    if [ -z "/pantherapi-pyclient" ]; then
        echo "panther installation not found, skipping go analysis"
    else    
        echo "generating genelist from mod table..."
        # produce gene lists for all GMUCT (for now) groups
        Rscript "$scripts"/produceGenelist.R \
            "$dir"/mod_long.csv \
            "$dir"/go/genelists

        # proceed if genelists directory is not empty
        if [ -n "$(ls "$dir"/go/genelists)" ]; then
            echo "sending each gene list to panther for overrepresentation analysis..."
            # Send each gene list into panther API and generate a overrepresentation result file in another folter
            for f in "$dir"/go/genelists/*.txt
            do
            n=$(basename "$f")
            echo "$n"
            python /pantherapi-pyclient/pthr_go_annots.py \
                --service enrich \
                --params_file /pantherapi-pyclient/params/enrich.json \
                --seq_id_file "$f" \
                > "$dir"/go/pantherout/"$n"
            done

            echo "producing heatmap..."
            # Run the R script that scavenges through a directory for result files and produce heatmap from it
            Rscript "$scripts"/panther2heatmap.R \
                "$dir"/go/pantherout \
                "$dir"
        fi
    fi

    echo "classifying modified RNA subtype..."
    # looking at RNA subtype for mods
    Rscript "$scripts"/RNAtype.R \
        "$dir"/mod_long.csv
    echo "done"
    echo ""

    if [ -e "$genomedir"/*_CDS.bed ] && [ -e "$genomedir"/*_fiveUTR.bed ] && [ -e "$genomedir"/*_threeUTR.bed ]; then
        c=$(find "$genomedir" -type f -name "*_CDS.bed")
        f=$(find "$genomedir" -type f -name "*_fiveUTR.bed")
        t=$(find "$genomedir" -type f -name "*_threeUTR.bed")
        echo "mapping modification regional distribution landscape..."
        # looking at RNA subtype for mods
        Rscript "$scripts"/modRegionMapping.R \
            "$dir"/mod_long.csv \
            "$f" \
            "$c" \
            "$t"
        echo "done"
        echo ""
    fi

    echo ""
    echo "#################################### HAMRLINC has finished running #######################################"
    date '+%d/%m/%Y %H:%M:%S'
    echo ""
fi
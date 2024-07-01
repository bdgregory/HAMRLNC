#!/bin/bash
set -u

# Harry Li, University of Pennsylvania & Chosen Obih, University of Arizona
# man page
usage () {
    echo ""
    echo "Usage : sh $0"
    echo ""

    cat <<'EOF'
  ######################################### COMMAND LINE OPTIONS #############################
  REQUIRED:
    -o  <project directory>
    -c  <filenames with key and name>
    -g  <reference genome.fa>
    -i  <reference genome annotation.gff3>
    -l  <read length>

  OPTIONAL: 
    -d  [input a directory of fastq]
    -D  [input a directory of bam]
    -b  [sort input bam files, default=false]
    -h  [help message] 
    -n  number of threads (default=4)
    -O  [Panther: organism taxon ID, default=3702]
    -A  [Panther: annotation dataset, default="GO:0008150"]
    -Y  [Panther: test type, default="FISHER"]
    -R  [Panther: correction type, default="FDR"]
    -y  [keep intermediate bam files, default=false]
    -q  [halt program upon completion of checkpoint 2, default=false]
    -G  [attribute used for featurecount, default="gene_id"]
    -k  [activate modification annotation workflow, default=false]
    -p  [activate lncRNA annotation workflow, default=false]
    -u  [activate featurecount workflow, default=false]
    -H  [SERVER alt path for panther]
    -U  [SERVER alt path for HAMRLINC]
    -W  [SERVER alt path for GATK]
    -S  [SERVER alt path for HAMR]
    -J  [SERVER alt path for CPC2]
    -M  [SERVER alt path for Rfam]
    -f  [HAMR: filter]
    -m  [HAMR: model]
    -Q  [HAMR: minimum quality score, default=30]
    -E  [HAMR: sequencing error, default=0.01]
    -P  [HAMR: maximum p-value, default=1]
    -F  [HAMR: maximum fdr, default=0.05]
    -C  [HAMR: minimum coverage, default=10]
  ################################################# END ########################################
EOF
    exit 0
}


############################################# Define Default Variables #####################################
# hamr related defaults
quality=30
coverage=10
err=0.01
pvalue=1
fdr=0.05
path_to_HAMR="/HAMR"
path_to_rfam=""

# hamr downstream
execpthr="/pantherapi-pyclient/pthr_go_annots.py"
execcpc="/CPC2/bin/CPC2.py"

# subprogram activation boolean
run_lnc=false
run_mod=false
run_featurecount=false
mod_partial=false

# other initialization
threads=4
attribute_fc="gene_id"
clean=true
#curdir=$(dirname "$0")
hsref=""
fastq_in=""
bam_in=""
bam_sorted=true
porg=""
pterm=""
ptest=""
pcorrect=""

# might change this later
hamrlinc_dir=""
gatk_dir=""


######################################################### Grab Arguments #########################################
while getopts ":o:c:g:i:l:d:D:bhn:O:A:Y:R:yqG:kupH:U:W:S:M:J:f:m:Q:E:P:F:C:" opt; do
  case $opt in
    o)
    out=$OPTARG # project output directory root
    ;;
    c)
    csv=$OPTARG # SRR to filename table
    ;;
    g)
    genome=$OPTARG # reference genome 
    ;;
    i)
    annotation=$OPTARG # reference genome annotation
    ;;
    l)
    length+=$OPTARG # read length 
    ;;
    f)
    filter=$OPTARG
    ;;
    m)
    model=$OPTARG
    ;;
    n)
    threads=$OPTARG
    ;;
    p)
    run_lnc=true
    ;;
    k)
    run_mod=true
    ;;
    u)
    run_featurecount=true
    ;;
    q)
    mod_partial=true
    ;;
    y)
    clean=false
    ;;
    Q)
    quality=$OPTARG
    ;;
    O)
    porg=$OPTARG
    ;;
    G)
    attribute_fc=$OPTARG
    ;;
    A)
    pterm=$OPTARG
    ;;
    Y)
    ptest=$OPTARG
    ;;
    R)
    pcorrect=$OPTARG
    ;;
    d)
    fastq_in=$OPTARG
    ;;
    D)
    bam_in=$OPTARG
    ;;
    b)
    bam_sorted=false
    ;;
    C)
    coverage=$OPTARG
    ;;
    E)
    err=$OPTARG
    ;;
    P)
    pvalue=$OPTARG
    ;;
    S)
    path_to_HAMR="$OPTARG"
    ;;
    M)
    path_to_rfam="$OPTARG"
    ;;
    H)
    execpthr="$OPTARG"
    ;;
    J)
    execcpc="$OPTARG"
    ;;
    F)
    fdr=$OPTARG
    ;;
    U)
    hamrlinc_dir=$OPTARG
    ;;
    W)
    gatk_dir=$OPTARG
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


############################################### Derive Other Variables #########################################
# reassign sample input files, genome and annotation files name and include file paths
user_dir=$(pwd)
genome="$user_dir"/"$genome"
annotation="$user_dir"/"$annotation"
out="$user_dir"/"$out"
csv="$user_dir"/"$csv"

# assigning additional variables
dumpout=$out/datasets
ttop=$((threads/2))
mismatch=$((length*6/100))
overhang=$((mismatch-1))
genomedir=$(dirname "$genome")
last_checkpoint=""

# HAMR components path assignment
exechamrpy="$path_to_HAMR"/"hamr.py"
execignoreends="$path_to_HAMR"/"ignoreBamReadEnds.py"
# this might not work?
export util="$hamrlinc_dir"/"util"
export scripts="$hamrlinc_dir"/"scripts"
export PATH="$gatk_dir/:$PATH"

filter="$util"/filter_SAM_number_hits.pl
model="$util"/euk_trna_mods.Rdata
json="$util"/panther_params.json
generator="$scripts"/annotationGenerateUnified.R


################################################ Subprogram Definitions #########################################
# announces reached checkpoint and updates checkpoint file, or creates txt if it didn't exist
checkpoint () {
    echo "Checkpoint reached: $1"
    echo "$1" > "$out"/checkpoint.txt
}

# is repeated for each accession code found in csv, performs fasterq-dump, fastqc, and trimming; automatic paired-end recognition
fastqGrabSRA () {

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

    # remove unneeded raw
    rm "$dumpout"/raw/"$line"."$suf"

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

    # remove unneeded raw
    rm "$dumpout"/raw/"$line""_1.$suf"
    rm "$dumpout"/raw/"$line""_2.$suf"
  fi

  echo "[$(date '+%d/%m/%Y %H:%M:%S')] finished processing $line"
  echo ""
}

# is repeated for each local file provided, performs fastqc and trimming; automatic paired-end recognition
fastqGrabLocal () {
    sname=$(basename "$fq")
    tt=$(echo "$sname" | cut -d'.' -f1)
    echo "[$sname] performing fastqc on raw file..."
    fastqc "$fq" -o "$dumpout"/fastqc_results &

    echo "[$sname] trimming..."
    trim_galore -o "$dumpout"/trimmed "$fq" --dont_gzip

    echo "[$sname] trimming complete, performing fastqc..."
    fastqc "$dumpout"/trimmed/"$tt""_trimmed.fq" -o "$dumpout"/fastqc_results
    
    # choosing not to remove user provided raw fastq
}

# called upon completion of each sorted BAM files, takes the file through pre-processing, and performs hamr
hamrBranch () {
    if [[ $currProg_mod == "2" ]]; then
        #adds read groups using picard, note the RG arguments are disregarded here
        echo "[$smpkey] adding/replacing read groups..."
        gatk AddOrReplaceReadGroups \
            I="$smpout"/sort_accepted.bam \
            O="$smpout"/sorted_RG.bam \
            RGPU=HWI-ST1395:97:d29b4acxx:8 \
            RGID=1 \
            RGSM=xxx \
            RGLB=xxx \
            RGPL=illumina 
        echo "[$smpkey] finished adding/replacing read groups (MOD 1/7)"
        echo ""

        # RG finished without exiting
        echo "3" > "$smpout"/progress_mod.txt
        currProg_mod="3"
    fi

    wait

    if [[ $currProg_mod == "3" ]]; then
        #filter the accepted hits by uniqueness
        echo "[$smpkey] filtering uniquely mapped reads..."
        samtools view \
            -h "$smpout"/sorted_RG.bam \
            | perl "$filter" 1 \
            | samtools view -bS - \
            | samtools sort \
            -o "$smpout"/sorted_RG_unique.bam
        echo "[$smpkey] finished filtering (MOD 2/7)"
        echo ""

        # filtering unique completed without erroring out if this is reached
        echo "4" > "$smpout"/progress_mod.txt
        currProg_mod="4"
    fi

    wait

    if [[ $currProg_mod == "4" ]]; then
        echo "[$smpkey] excluding read ends..."
        # first index the input bam file with samtools
        samtools index \
            "$smpout"/sorted_RG_unique.bam \
            -o "$smpout"/sorted_RG_unique.bai

        wait

        # ignore read ends
        python $execignoreends \
            -5p 1 -3p 1 \
            "$smpout"/sorted_RG_unique.bam \
            "$smpout"/sorted_RG_unique_endsIGN.bam
        echo "[$smpkey] finished excluding (MOD 3/7)"
        echo ""

        # removing ends finished without exiting
        echo "5" > "$smpout"/progress_mod.txt
        currProg_mod="5"
    fi 

    wait


    if [[ $currProg_mod == "5" ]]; then
        #reorder the reads using picard
        # 6/26 adding allow discordance, remove ends might have messed with the length?
        echo "[$smpkey] reordering..."
        gatk --java-options "-Xmx2g -Djava.io.tmpdir=$smpout/tmp" ReorderSam \
            I="$smpout"/sorted_RG_unique_endsIGN.bam \
            O="$smpout"/sorted_RG_unique_endsIGN_reordered.bam \
            R="$genome" \
            CREATE_INDEX=TRUE \
            SEQUENCE_DICTIONARY="$dict" \
            TMP_DIR="$smpout"/tmp \
            ALLOW_CONTIG_LENGTH_DISCORDANCE=true
        echo "[$smpkey] finished reordering (MOD 4/7)"
        echo ""

        # ordering finished without exiting
        echo "6" > "$smpout"/progress_mod.txt
        currProg_mod="6"
    fi 

    wait

    if [[ $currProg_mod == "6" ]]; then
        #splitting and cigarring the reads, using genome analysis tool kit
        #note can alter arguments to allow cigar reads 
        echo "[$smpkey] getting split and cigar reads..."
        gatk --java-options "-Xmx2g -Djava.io.tmpdir=$smpout/tmp" SplitNCigarReads \
            -R "$genome" \
            -I "$smpout"/sorted_RG_unique_endsIGN_reordered.bam \
            -O "$smpout"/sorted_RG_unique_endsIGN_reordered_SNC.bam \
            # apparently this is now outdated
            # -U ALLOW_N_CIGAR_READS
        echo "[$smpkey] finished splitting N cigarring (MOD 5/7)"
        echo ""

        # cigaring and spliting finished without exiting
        echo "7" > "$smpout"/progress_mod.txt
        currProg_mod="7"
    fi 

    wait

    if [[ $currProg_mod == "7" ]]; then
        #final resorting using picard
        echo "[$smpkey] resorting..."
        gatk --java-options "-Xmx2g -Djava.io.tmpdir=$smpout/tmp" SortSam \
            I="$smpout"/sorted_RG_unique_endsIGN_reordered_SNC.bam \
            O="$smpout"/sorted_RG_unique_endsIGN_reordered_SNC_resorted.bam \
            SORT_ORDER=coordinate
        echo "[$smpkey] finished resorting (MOD 6/7)"
        echo ""

        # cigaring and spliting finished without exiting
        echo "8" > "$smpout"/progress_mod.txt
        currProg_mod="8"
    fi 

    wait

    if [[ $currProg_mod == "8" ]]; then
        #hamr step, can take ~1hr
        echo "[$smpkey] hamr..."
        #hamr_path=$(which hamr.py) 
        python $exechamrpy \
            -fe "$smpout"/sorted_RG_unique_endsIGN_reordered_SNC_resorted.bam "$genome" "$model" "$smpout" $smpname $quality $coverage $err H4 $pvalue $fdr .05
        wait

        if [ ! -e "$smpout/${smpname}.mods.txt" ]; then 
            cd "$hamrout" || exit
            printf '%s \n' "$smpname" >> zero_mod.txt
            cd || exit
        else
        # HAMR needs separate folders to store temp for each sample, so we move at the end
            cp "$smpout"/"${smpname}".mods.txt "$hamrout"
        fi
        echo "9" > "$smpout"/progress_mod.txt
        currProg_mod="9"
    fi
}

# called upon completion of each sorted BAM files, takes the sorted BAM through lncRNA prediction pipeline
lncCallBranch () {

    cd $smpout

    if [[ $currProg_lnc == "2" ]]; then
        echo "[$smpkey] running stringtie for gtf conversion..."

        # turn bam into gtf
        stringtie \
            sort_accepted.bam \
            -G $annotation -o stringtie_out.gtf \
            -f 0.05 -j 9 -c 7 -s 20

        echo "[$smpkey] finished converting bam to gtf (LNC 1/15)"
        echo ""

        echo "3" > "$smpout"/progress_lnc.txt
        currProg_lnc="3"
    fi
    
    if [[ $currProg_lnc == "3" ]]; then
        echo "[$smpkey] merging sample gtf with reference gtf..."
        
        # merge gtf from bam with ref gtf
        stringtie --merge -G $annotation \
            -o stringtie_merge_out.gtf \
            stringtie_out.gtf

        echo "[$smpkey] finished merging (LNC 2/15)"
        echo ""

        echo "4" > "$smpout"/progress_lnc.txt
        currProg_lnc="4"
    fi


    if [[ $currProg_lnc == "4" ]]; then
        echo "[$smpkey] running gffcompare on merged gtf..."
        
        # gffcmp merged gtf
        gffcompare -r $annotation stringtie_merge_out.gtf

        echo "[$smpkey] finished gffcompare (LNC 3/15)"
        echo ""

        echo "5" > "$smpout"/progress_lnc.txt
        currProg_lnc="5"
    fi

    if [[ $currProg_lnc == "5" ]]; then
        echo "[$smpkey] filtering..."
        
        # filter with awk
        awk '$7 != "." {print}' gffcmp.annotated.gtf > filtered_gffcmp_annotated.gtf

        echo "[$smpkey] finished filtering (LNC 4/15)"
        echo ""

        echo "6" > "$smpout"/progress_lnc.txt
        currProg_lnc="6"
    fi


    if [[ $currProg_lnc == "6" ]]; then
        echo "[$smpkey] filtering for UX..."
        
        # filter with grep for UX class codes
        grep -E 'class_code "u";|class_code "x";' filtered_gffcmp_annotated.gtf > UXfiltered_gffcmp_annotated.gtf

        echo "[$smpkey] finished filtering (LNC 6/15)"
        echo ""

        echo "7" > "$smpout"/progress_lnc.txt
        currProg_lnc="7"
    fi


    if [[ $currProg_lnc == "7" ]] && [[ -s UXfiltered_gffcmp_annotated.gtf ]]; then
        echo "[$smpkey] creating index file..."
        
        # I could copy over the already made fai but... eh
        samtools faidx $genome

        echo "[$smpkey] finished indexing (LNC 7/15)"
        echo ""

        echo "8" > "$smpout"/progress_lnc.txt
        currProg_lnc="8"
    fi

    
    if [[ $currProg_lnc == "8" ]]; then
        echo "[$smpkey] converting filtered gtf to gff3..."
        
        # covnvert to gff3
        gffread UXfiltered_gffcmp_annotated.gtf -T -o UXfiltered_gffcmp_annotated.gff3

        echo "[$smpkey] finished conversion (LNC 8/15)"
        echo ""

        echo "9" > "$smpout"/progress_lnc.txt
        currProg_lnc="9"
    fi


    if [[ $currProg_lnc == "9" ]]; then
        echo "[$smpkey] writing fa file from filtered gtf..."
        
        # write gtf to fasta
        gffread UXfiltered_gffcmp_annotated.gtf -g $genome -w transcripts.fa

        echo "[$smpkey] finished writing (LNC 9/15)"
        echo ""

        echo "10" > "$smpout"/progress_lnc.txt
        currProg_lnc="10"
    fi


    if [[ $currProg_lnc == "10" ]]; then
        echo "[$smpkey] analyzing for transcript coding probability with cpc2..."
        
        # use cpc2 to analyze for coding probablity
        python $execcpc -i transcripts.fa -o cpc2_output

        echo "[$smpkey] finished analysis (LNC 10/15)"
        echo ""

        echo "11" > "$smpout"/progress_lnc.txt
        currProg_lnc="11"
    fi


    if [[ $currProg_lnc == "11" ]]; then
        echo "[$smpkey] extracting transcripts with probability less than 0.5..."
        
        # use awk to extract result
        awk '$7 < 0.5' cpc2_output.txt > filtered_transcripts.txt

        echo "[$smpkey] finished extraction (LNC 11/15)"
        echo ""

        echo "12" > "$smpout"/progress_lnc.txt
        currProg_lnc="12"
    fi


    if [[ $currProg_lnc == "12" ]]; then
        echo "[$smpkey] using cpc2 results to filter gtf..."
        
        # loop to fetch entries from gtf
        inputFile="$smpout/filtered_transcripts.txt"
        gtfFile="$smpout/UXfiltered_gffcmp_annotated.gtf"
        outputFile="$smpout/cpc_filtered_transcripts.txt"
        while IFS= read -r line; do
            pattern=$(echo "$line" | cut -f1)
            grep "$pattern" "$gtfFile" >> "$outputFile"
        done < "$inputFile"

        echo "[$smpkey] finished filtering (LNC 12/15)"
        echo ""

        echo "13" > "$smpout"/progress_lnc.txt
        currProg_lnc="13"
    fi

    
    if [[ $currProg_lnc == "13" ]]; then
        echo "[$smpkey] creating fa file from cpc2 filtered gtf..."
        
        # create fa from cpc filtered gtf
        gffread cpc_filtered_transcripts.txt -g $genome -w rfam_in.fa

        echo "[$smpkey] finished writing (LNC 13/15)"
        echo ""

        echo "14" > "$smpout"/progress_lnc.txt
        currProg_lnc="14"
    fi


    if [[ $currProg_lnc == "14" ]]; then
        echo "[$smpkey] performing cmscan..."
        
        # create fa from cpc filtered gtf
        cmscan --nohmmonly \
            --rfam --cut_ga --fmt 2 --oclan --oskip \
            --clanin "$path_to_rfam"/Rfam.clanin -o "$smpout"/my.cmscan.out \
            --tblout "$smpout"/my.cmscan.tblout "$path_to_rfam"/Rfam.cm "$smpout"/rfam_in.fa

        echo "[$smpkey] finished (LNC 14/15)"
        echo ""

        echo "15" > "$smpout"/progress_lnc.txt
        currProg_lnc="15"
    fi


    if [[ $currProg_lnc == "15" ]]; then
        echo "[$smpkey] sacnning against rfam and finishing lncRNA annotation..."
        
        # tblout info extraction 
        inputFile="$smpout/my.cmscan.tblout"
        gtfFile="$smpout/cpc_filtered_transcripts.txt"
        outputFile="$smpout/rfam_filtered_transcripts.txt"
        tail -n +3 "$inputFile" >> parsed_rfam_out.tblout

        # created a python script to deal with infernal's space delimited file
        cp "cpc_filtered_transcripts.txt" "rfam_filtered_transcripts.txt"
        while IFS= read -r line; do
            if [[ $line =~ ^#$ ]]; then 
                break
            else
                pattern=$(python "$scripts"/parser.py "$line")
                sed -i "/$pattern/d" "$outputFile"
            fi
        done < "$smpout"/parsed_rfam_out.tblout

        # rename, combine
        cat $annotation "$smpout"/rfam_filtered_transcripts.txt > "$smpout"/final_combined.gtf
        mv rfam_filtered_transcripts.txt "${smpname}".lnc.gtf

        echo "[$smpkey] finished lncRNA annotation (LNC 15/15)"
        echo ""

        echo "16" > "$smpout"/progress_lnc.txt
        currProg_lnc="16"
    fi

    cd
    
    if [[ $currProg_lnc == "16" ]]; then
        echo "[$smpkey] processing identified lncRNA into GTF..."
        Rscript "$scripts"/lnc_processing.R \
            "$smpout"/rfam_filtered_transcripts.txt \
            "$smpout"/"${smpname}"

        cp "$smpout"/"${smpname}".lnc.gtf "$lncout"
    elif [[ $currProg_lnc == "7" ]]; then
        echo "[$smpkey] filtering by class codes -u -x in gtf yielded no entries, no lncRNA can be annotated"
    fi

    echo "[$smpkey] done (LNC)"
    echo ""
}

# called upon completion of lncCall, performs abundance analysis for each BAM dependning on lnc arm
featureCountBranch () {
    echo "[$(date '+%d/%m/%Y %H:%M:%S')$smpkey] quantifying regular transcript abundance using featurecounts..."
    if [ ! -d "$out/featurecount_out" ]; then mkdir "$out/featurecount_out"; fi

    if [[ "$run_lnc" = true ]]; then
        # if lncRNA annotated we also feature count with the combined gtf, separate by PE det
        echo "[$smpkey] quantifying transcripts found in reads with lncRNA..."
        echo "[$smpkey] running featurecount..."
        featureCounts \
            -T 2 \
            -t transcript \
            -g $attribute_fc \
            -a "$out"/final_combined.gtf \
            -o "$smpout"/"$smpname"_transcript_abundance_lncRNA-included.txt \
            "$smpout"/sort_accepted.bam

        # housekeeping for lnc abundance
        mv "$smpout"/"$smpname"_transcript_abundance_lncRNA-included.txt "$out/featurecount_out"
        echo "[$smpkey] finished quantifying read features"
    fi

    # always do feature count with the regular gtf
    # first create gtf file from gff3 file
    gffread \
        "$annotation" \
        -T \
        -o "$out"/temp.gtf
    
    echo "[$smpkey] quantifying exons found in reads..."
    echo "[$smpkey] running featurecount..."
    featureCounts \
        -T 2 \
        -t exon \
        -g $attribute_fc \
        -a "$out"/temp.gtf \
        -o "$smpout"/"$smpname"_exon_abundance.txt \
        "$smpout"/sort_accepted.bam
    echo "[$smpkey] finished quantifying read features"

    # housekeeping for regular abundance
    mv "$smpout"/"$smpname"_exon_abundance.txt "$out/featurecount_out"
}

# the wrapper around hamrBranch and lncCallBranch, is called once for each rep (or input file)
fastq2raw () {

    # Read the CSV file into a DataFrame
    mapfile -t names < <(awk -F, '{ print $1 }' "$csv")
    mapfile -t smpf < <(awk -F, '{ print $2 }' "$csv")

    # Create a dictionary from the DataFrame
    declare -A dictionary
    for ((i=0; i<${#names[@]}; i++)); 
    do
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

    # shared info for all 3 arms
    smpout=$out/pipeline/$smpkey"_temp"
    echo "[$smpkey] You can find all the intermediate files for $smpkey at $smpout" 

    # sort out directories and progress info for mod
    currProg_mod="0"
    if [[ "$run_mod" = true ]]; then
        # Reassign hamr output directory
        if [ ! -d "$out/hamr_out" ]; then
            mkdir "$out"/hamr_out
            echo "created path: $out/hamr_out"
        fi
        hamrout=$out/hamr_out
        echo "[$smpkey] You can find the HAMR output file for $smpkey at $hamrout/$smpname.mod.txt"

        # check if progress_mod.txt exists, if not, create it with 0
        if [[ ! -e "$smpout"/progress_mod.txt ]]; then
            echo "0" > "$smpout"/progress_mod.txt
        fi

        currProg_mod=$(cat "$smpout"/progress_mod.txt)
        echo "-------------Folder $smpkey is at progress number $currProg_mod for this run of modification annotation--------------"
    fi

    # sort out directories and progress info for lnc
    currProg_lnc="0"
    if [[ "$run_lnc" = true ]]; then
        # Reassign lnc output directory
        if [ ! -d "$out/lnc_out" ]; then
            mkdir "$out"/lnc_out
            echo "created path: $out/lnc_out"
        fi
        lncout=$out/lnc_out
        echo "[$smpkey] You can find the lncRNA output file for $smpkey at $lncout/$smpname.mod.txt" 

        # check if progress_lnc.txt exists, if not, create it with 0
        if [[ ! -e "$smpout"/progress_lnc.txt ]]; then
            echo "0" > "$smpout"/progress_lnc.txt
        fi

        currProg_lnc=$(cat "$smpout"/progress_lnc.txt)
        echo "-------------Folder $smpkey is at progress number $currProg_lnc for this run of lncRNA annotation--------------"
    fi

    echo "$(date '+%d/%m/%Y %H:%M:%S') [$smpkey] Begin preprocessing pipeline"

    
    # if 0, then either this run failed before mapping completion or this run just started
    if [[ $currProg_mod == "0" || $currProg_lnc == "0" ]]; then
        if [[ -z $bam_in ]]; then
            cd "$smpout" || exit
            # maps the trimmed reads to provided annotated genome, can take ~1.5hr
            echo "--------Entering mapping step--------"
            if [ "$det" -eq 1 ]; then
                echo "[$smpkey] Using STAR for mapping with a single-end file."
                STAR \
                --runThreadN 2 \
                --genomeDir "$out"/STARref \
                --readFilesIn "$smp" \
                --sjdbOverhang $overhang \
                --sjdbGTFfile "$annotation" \
                --sjdbGTFtagExonParentTranscript Parent \
                --outFilterMultimapNmax 10 \
                --outFilterMismatchNmax $mismatch \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMstrandField intronMotif \
                --outFilterIntronMotifs RemoveNoncanonical
            else
                echo "[$smpkey] Using STAR for mapping with a paired-end file."
                STAR \
                --runThreadN 2 \
                --genomeDir "$out"/STARref \
                --readFilesIn "$smp1" "$smp2" \
                --sjdbOverhang $overhang \
                --sjdbGTFfile "$annotation" \
                --sjdbGTFtagExonParentTranscript Parent \
                --outFilterMultimapNmax 10 \
                --outFilterMismatchNmax $mismatch \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMstrandField intronMotif \
                --outFilterIntronMotifs RemoveNoncanonical
            fi
            cd || exit
        else 
            echo "relocating provided bam file..."
            if [[ "$bam_sorted" = true ]]; then
                mv "$smp" "$smpout/sort_accepted.bam"
            else
                echo "user requested sorting provided bam files..."
                samtools sort \
                    -n "$smp" \
                    -o "$smpout"/sort_accepted.bam
                echo "[$smpkey] finished sorting"
            fi
        fi

        # mapping completed without erroring out if this is reached
        echo "1" > "$smpout"/progress_mod.txt
        echo "1" > "$smpout"/progress_lnc.txt

        # update directly for normal progression
        currProg_mod="1"
        currProg_lnc="1"
    fi

    wait

    # if 1, then either last run failed before sorting completion or this run just came out of mapping
    if [[ $currProg_mod == "1" || $currProg_lnc == "1" ]]; then
        if [[ -z $bam_in ]]; then
            echo "[$smpkey] renaming file"
            # star already sorts its output
            mv "$smpout"/Aligned.sortedByCoord.out.bam "$smpout"/sort_accepted.bam
            echo "[$smpkey] finished"
            echo ""
        else
            echo "@@@@ bam input mode has now merged track with fastq mode"
            echo "@@@@ developer can go get a coffee now"
        fi

        # sorting completed without erroring out if this is reached
        echo "2" > "$smpout"/progress_mod.txt
        echo "2" > "$smpout"/progress_lnc.txt
        currProg_mod="2"
        currProg_lnc="2"
    fi

    wait

    # if both lnc and mod are true, then run them parallelized
    if [[ "$run_lnc" = true ]] && [[ "$run_mod" = true ]]; then
        hamrBranch &
        lncCallBranch
    # this means only lnc runs, mod doesn't
    elif [[ "$run_lnc" = true ]]; then
        lncCallBranch
    # this means only mod runs, lnc doesn't
    elif [[ "$run_mod" = true ]]; then
        hamrBranch
    fi

    wait

    if [[ "$run_featurecount" = true ]]; then
        featureCountBranch
    fi    

    wait

    # intermediate file clean up
    if [[ $currProg_mod == "9" ]]; then
        # Move the unique_RG_ordered.bam and unique_RG_ordered.bai to a folder for read depth analysis
        cp "$smpout"/sorted_RG_unique_endsIGN_reordered.bam "$out"/pipeline/depth/"$smpname".bam
        cp "$smpout"/sorted_RG_unique_endsIGN_reordered.bai "$out"/pipeline/depth/"$smpname".bai

        # only delete intermediates if user didn't disable it
        if [[ "$clean" == true ]]; then
            # delete more intermediate files
            echo "[$smpkey] removing large intermediate files..."
            rm "$smpout"/sorted_RG.bam
            rm "$smpout"/sorted_RG_unique.bam
            rm "$smpout"/sorted_RG_unique_endsIGN.bam
            rm "$smpout"/sorted_RG_unique_endsIGN_reordered.bam
            rm "$smpout"/sorted_RG_unique_endsIGN_reordered_SNC.bam
            rm "$smpout"/sorted_RG_unique_endsIGN_reordered_SNC_resorted.bam
        fi

        echo "[$smpkey] finished cleaning (MOD 7/7)"
        echo "10" > "$smpout"/progress_mod.txt
    fi
}

# the wrapper around fastq2raw, iterates through existing files so fastq2raw can be called once for each file
parallelWrap () {
    smpext=$(basename "$smp")
    smpdir=$(dirname "$smp")
    smpkey="${smpext%.*}"
    smpname=""
    original_ext="${smpext##*.}"

    # this function will be split into bam route and fastq route, 6/30/24
    if [[ -z $bam_in ]]; then
        # always run the below to ensure necessary variables are assigned
        if [[ $smpkey == *_1 ]]; then
            smpkey="${smpkey%_1*}"
            smp1="$smpdir/${smpkey}_1_trimmed.$original_ext"
            smp2="$smpdir/${smpkey}_2_trimmed.$original_ext"
            # Paired end recognized
            det=0
            echo "$smpext is a part of a paired-end sequencing file"
            fastq2raw
        elif [[ $smpkey == *_2 ]]; then
            # If _2 is in the filename, this file was processed along with its corresponding _1 so we skip
            echo "$smpext has already been processed with its _1 counter part. Skipped."
            echo ""
        else
            det=1
            echo "$smpext is a single-end sequencing file"
            fastq2raw
        fi
    else
        echo "$smpext is a bam file, HAMRLINC will shortcut into appropriate steps"
        fastq2raw
    fi
}

# de novo function to deal with consensus modifications
consensusOverlap () {
    IFS="/" read -ra sections <<< "$smp"
    temp="${sections[-1]}"

    IFS="." read -ra templ <<< "$temp"
    smpname="${templ[0]}"

    echo "consensus file prefix: $smpname"

    count=$(ls -1 "$out"/annotBeds/*_CDS.bed 2>/dev/null | wc -l)
    if [ "$count" != 0 ]; then 
        cds=$(find "$out"/annotBeds -maxdepth 1 -name "*_CDS.bed")
        #overlap with cds
        intersectBed \
            -a "$cds" \
            -b "$smp" \
            -wa -wb \
            > "$out"/lap/"$smpname""_CDS".bed
        echo "finished finding overlap with CDS library"
    fi

    count=$(ls -1 "$out"/annotBeds/*_fiveUTR.bed 2>/dev/null | wc -l)
    if [ "$count" != 0 ]; then 
        fiveutr=$(find "$out"/annotBeds -maxdepth 1 -name "*_fiveUTR.bed")
        #overlap with 5utr
        intersectBed \
            -a "$fiveutr" \
            -b "$smp" \
            -wa -wb \
            > "$out"/lap/"$smpname""_fiveUTR".bed
        echo "finished finding overlap with 5UTR library"
    fi

    count=$(ls -1 "$out"/annotBeds/*_threeUTR.bed 2>/dev/null | wc -l)
    if [ "$count" != 0  ]; then 
        threeutr=$(find "$out"/annotBeds -maxdepth 1 -name "*_threeUTR.bed")
        #overlap with 3utr
        intersectBed \
            -a "$threeutr" \
            -b "$smp" \
            -wa -wb \
            > "$out"/lap/"$smpname""_threeUTR".bed
        echo "finished finding overlap with 3UTR library"
    fi

    count=$(ls -1 "$out"/annotBeds/*_gene.bed 2>/dev/null | wc -l)
    if [ "$count" != 0 ]; then 
        gene=$(find "$out"/annotBeds -maxdepth 1 -name "*_gene.bed")
        #overlap with gene
        intersectBed \
            -a "$gene" \
            -b "$smp" \
            -wa -wb \
            > "$out"/lap/"$smpname""_gene".bed
        echo "finished finding overlap with gene library"
    fi

    count=$(ls -1 "$out"/annotBeds/*_primarymRNA.bed 2>/dev/null | wc -l)
    if [ "$count" != 0 ]; then 
        mrna=$(find "$out"/annotBeds -maxdepth 1 -name "*_primarymRNA.bed")
        #overlap with mrna
        intersectBed \
            -a "$mrna" \
            -b "$smp" \
            -wa -wb \
            > "$out"/lap/"$smpname""_primarymRNA".bed
        echo "finished finding overlap with primary mRNA library"
    fi

    count=$(ls -1 "$out"/annotBeds/*_exon.bed 2>/dev/null | wc -l)
    if [ "$count" != 0 ]; then 
        exon=$(find "$out"/annotBeds -maxdepth 1 -name "*_exon.bed")
        #overlap with exon
        intersectBed \
            -a "$exon" \
            -b "$smp" \
            -wa -wb \
            > "$out"/lap/"$smpname""_exon".bed
        echo "finished finding overlap with exon library"
    fi

    count=$(ls -1 "$out"/annotBeds/*_ncRNA.bed 2>/dev/null | wc -l)
    if [ "$count" != 0 ]; then 
        nc=$(find "$out"/annotBeds -maxdepth 1 -name "*_ncRNA.bed")
        #overlap with nc rna
        intersectBed \
            -a "$nc" \
            -b "$smp" \
            -wa -wb \
            > "$out"/lap/"$smpname""_ncRNA".bed
        echo "finished finding overlap with ncRNA library"
    fi

    ######## this is from the lncRNA identification steps ##########
    # given lines of lncRNA region in gtf, see if any mod can be found there
    if [[ "$run_lnc" = true ]]; then
        # intersectBed \
        #     -a "$smpout"/"$smpname".lnc.gtf \
        #     -b "$smp" \
        #     -wa -wb \
        #     > "$out"/lap/"$smpname"_overlapped_lnc.bed
        echo "finished finding overlap with lncRNA predictions"
    fi

    echo ""
}

# house keeping steps for fastqGrab functions, mostly creating folders and checking function calls
fastqGrabHouseKeeping () {
    ##########fastqGrab housekeeping begins#########
    if [ ! -d "$out/datasets" ]; then mkdir "$out"/datasets; echo "created path: $out/datasets"; fi

    # first see whether input folder is provided
    if [[ ! -z $fastq_in ]]; then
        fastq_in="$user_dir"/"$fastq_in"
        echo "Directory $fastq_in is found, assuming raw fastq files are provided..."
        mode=2
    else
        # Create directory to store original fastq files
        if [ ! -d "$out/datasets/raw" ]; then mkdir "$out"/datasets/raw; fi
        echo "You can find your original fastq files at $out/datasets/raw" 
        mode=1
        # grab txt from csv
        awk -F "," '{print $1}' $csv > "$user_dir"/accession.txt
        acc="$user_dir"/accession.txt
    fi

    if [ ! -d "$out/filein" ]; then 
        mkdir "$out/filein"
        echo "created path: $out/filein"
        # keeps a reference for the user inputted required files
        cp $genome "$out/filein"
        cp $annotation "$out/filein"
        cp $csv "$out/filein"
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
    ##########fastqGrab housekeeping ends#########
}

# house keeping steps for fastq2raw, mostly creating folders, some indices, and checking function calls
fastq2rawHouseKeeping () {
    ############fastq2raw housekeeping begins##############
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
    # previous code didn't work when the genome dir contained another genome with its dict
    fafilename=$(basename "$genome")
    fafilestem="${fafilename%.*}"
    dict="$genomedir"/"$fafilestem".dict
    if [[ ! -f "$dict" ]]; then
        gatk CreateSequenceDictionary \
        R="$genome"
    fi

    # create fai index file using fasta genome
    # same as above, note fai by convention has .fa in file name
    fai="$genomedir"/"$fafilestem".fa.fai
    if [[ ! -f "$dict" ]]; then
        samtools faidx "$genome"
    fi

    # Check if indexed files already present for STAR
    if [ -e "$out/STARref/SAindex" ]; then
        echo "STAR Genome Directory with indexed genome detected, skipping STAR indexing"
    else
        # get genome length
        genomelength=$(bioawk -c fastx '{ print length($seq) }' < $genome | awk '{sum += $1} END {print sum}')
        echo "For reference, your provided genome length is $genomelength long"

        # Define the SA index number argument
        log_result=$(echo "scale=2; l($genomelength)/l(2)/2 - 1" | bc -l)
        sain=$(echo "scale=0; if ($log_result < 14) $log_result else 14" | bc)
        echo "Creating STAR genome index..."
        # Create genome index 
        STAR \
            --runThreadN $threads \
            --runMode genomeGenerate \
            --genomeDir "$out"/STARref \
            --genomeFastaFiles "$genome" \
            --sjdbGTFfile "$annotation" \
            --sjdbGTFtagExonParentTranscript Parent \
            --sjdbOverhang $overhang \
            --genomeSAindexNbases $sain
    fi

    # Run a series of command checks to ensure fastq2raw can run smoothly
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
    #############fastq2raw housekeeping ends#############
}

# house keeping steps for when the user inputs bam files instead of fastq
bamEntranceHouseKeeping () {

    bam_in="$user_dir"/"$bam_in"
    if [[ -d "$bam_in" ]]; then
        echo "Directory $bam_in is found, assuming mapped bam files are provided..."
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

    # Creates a folder for depth analysis
    if [ ! -d "$out/pipeline/depth" ]; then mkdir "$out"/pipeline/depth; echo "created path: $out/pipeline/depth"; fi

    # create dict file using fasta genome file
    # previous code didn't work when the genome dir contained another genome with its dict
    fafilename=$(basename "$genome")
    fafilestem="${fafilename%.*}"
    dict="$genomedir"/"$fafilestem".dict
    if [[ ! -f "$dict" ]]; then
        gatk CreateSequenceDictionary \
        R="$genome"
    fi

    # create fai index file using fasta genome
    # same as above, note fai by convention has .fa in file name
    fai="$genomedir"/"$fafilestem".fa.fai
    if [[ ! -f "$dict" ]]; then
        samtools faidx "$genome"
    fi


    if ! command -v gatk > /dev/null; then
        echo "Failed to call gatk command. Please check your installation."
        exit 1
    fi

    if ! command -v mapfile > /dev/null; then
        echo "Failed to call mapfile command. Please check your installation."
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
}

# house keeping steps before starting the main program, checks key arguments, set checkpoints, etc
mainHouseKeeping () {
    # Check if the required arguments are provided
    if [ -z "$out" ]; then 
        echo "output directory not detected, exiting..."
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
    else
        echo "all required arguments provided, proceding..."
    fi

    # check that the user didn't suppress all three programs -- if so, there's no need to run anything
    if [[ "$run_lnc" = false ]] && [[ "$run_featurecount" = false ]] && [[ "$run_mod" = false ]]; then
        echo "User has not activated any functionalities. Exiting..."
        exit 0
    fi

    # check if run checkpoint.txt exists, if not, create it with start
    if [[ ! -e "$out"/checkpoint.txt ]]; then
        touch "$out"/checkpoint.txt
        echo "start" > "$out"/checkpoint.txt
    fi

    # determine stage of progress for this sample folder at this run
    last_checkpoint=$(cat "$out"/checkpoint.txt)
    echo "-------------------------------------------"
    echo "Starting from checkpoint: $last_checkpoint"
    echo "-------------------------------------------"
}

######################################################### Main Program #########################################
echo ""
echo "##################################### Begin HAMRLINC #################################"
echo ""

if [ ! -d "$out" ]; then mkdir "$out"; echo "created path: $out"; fi

# designate log file, if exist, clear, and have all stdout written
logstart=$(date "+%Y.%m.%d-%H.%M.%S")
logfile=$out/Log_$logstart.log
exec > >(tee -a "$logfile") 2>&1
#below captures only echo...?
#exec 2>&1 1>>$logfile 3>&1

# perform house keeping steps
mainHouseKeeping

if [[ -z $bam_in ]]; then
    # run fastqGrab when checkpoint is at start
    if [ "$last_checkpoint" = "start" ] || [ "$last_checkpoint" = "" ]; then
        fastqGrabHouseKeeping
        ##########fastqGrab main begins#########
        if [[ $mode -eq 1 ]]; then
            # Grabs the fastq files from acc list provided into the dir ~/datasets
            i=0
            while IFS= read -r line
            do ((i=i%threads)); ((i++==0)) && wait
            fastqGrabSRA &
            done < "$acc"

        elif [[ $mode -eq 2 ]]; then
            i=0
            for fq in "$fastq_in"/*; 
            do
                ((i=i%threads)); ((i++==0)) && wait
                fastqGrabLocal &
            done
        fi
        wait
        ##########fastqGrab main ends############
        echo ""
        echo "################ Finished downloading and processing all fastq files. Entering pipeline for HAMR analysis. ######################"
        date '+%d/%m/%Y %H:%M:%S'
        echo ""

        # obtained all processed fastq files, record down checkpoint
        last_checkpoint="checkpoint1"
        checkpoint $last_checkpoint
    fi

    # run fastq2raw if program is at checkpoint 1
    if [ "$last_checkpoint" = "checkpoint1" ]; then 
        fastq2rawHouseKeeping
        #############fastq2raw main begins###############
        # Pipes each fastq down the hamr pipeline, and stores out put in ~/hamr_out
        # Note there's also a hamr_out in ~/pipeline/SRRNUMBER_temp/, but that one's for temp files
        
        #mkdir trimmed_temp && mv "$hamrin"/*."$suf" trimmed_temp && chmod -R 777 trimmed_temp
        #cd trimmed_temp
        #current_dir=$(pwd)
        #cd ..

        i=0
        for smp in "$hamrin"/*."$suf"; 
        do
            ((i=i%ttop)); ((i++==0)) && wait   
            parallelWrap &
        done

        wait

        # these checks apply only if mod arm was on
        if [[ "$run_mod" = true ]]; then
            # Check whether any hamr.mod.text is present, if not, halt the program here
            if [[ -z "$(ls -A "$out"/hamr_out)" ]]; then
                echo "No HAMR predicted mod found for any sequencing data in this project, please see log for verification"
                exit 1
            else
                #at least 1 mod file, move zero mod record outside so it doesn't get read as a modtbl next
                mv "$out"/hamr_out/zero_mod.txt "$out"
            fi
        fi

        echo ""
        echo "################ Finished the requested analysis for each fastq file. Now producing consensus files and depth analysis. ######################"
        echo "$(date '+%d/%m/%Y %H:%M:%S')"
        echo ""

        #############fastq2raw main ends###############

        # obtained all HAMR / lnc results, record down checkpoint
        last_checkpoint="checkpoint2"
        checkpoint $last_checkpoint
    fi
else
    # set up project folder
    bamEntranceHouseKeeping

    # declare cycle patterns
    suf="bam"

    # main cycle through supplied bam folder
    i=0
    for smp in "$bam_in"/*."$suf"; 
    do
        ((i=i%ttop)); ((i++==0)) && wait   
        parallelWrap &
    done

    wait

    # these checks apply only if mod arm was on
    if [[ "$run_mod" = true ]]; then
        # Check whether any hamr.mod.text is present, if not, halt the program here
        if [[ -z "$(ls -A "$out"/hamr_out)" ]]; then
            echo "No HAMR predicted mod found for any sequencing data in this project, please see log for verification"
            exit 1
        else
            #at least 1 mod file, move zero mod record outside so it doesn't get read as a modtbl next
            mv "$out"/hamr_out/zero_mod.txt "$out"
        fi
    fi

    echo ""
    echo "################ Finished the requested analysis for each fastq file. Now producing consensus files and depth analysis. ######################"
    echo "$(date '+%d/%m/%Y %H:%M:%S')"
    echo ""

    #############fastq2raw main ends###############

    # obtained all HAMR / lnc results, record down checkpoint
    last_checkpoint="checkpoint2"
    checkpoint $last_checkpoint
fi

# check for mod_partial flag, if not, run consensus when checkpoint is at 2
if [[ "$mod_partial" = true ]]; then 
    echo "User has selected to perform only partial HAMRLINC functions"
    echo ""
    echo "#################################### HAMRLINC has finished running #######################################"
    date '+%d/%m/%Y %H:%M:%S'
    echo ""
elif [ "$last_checkpoint" = "checkpoint2" ]; then 
    ##############consensus finding begins##############
    # Produce consensus bam files based on filename (per extracted from name.csv) and store in ~/consensus
    if [ ! -d "$out/hamr_consensus" ]; then mkdir "$out"/hamr_consensus; echo "created path: $out/hamr_consensus"; fi
    if [ ! -d "$out/lnc_consensus" ]; then mkdir "$out"/lnc_consensus; echo "created path: $out/lnc_consensus"; fi

    # Run a series of command checks to ensure findConsensus can run smoothly
    if ! command -v Rscript > /dev/null; then
        echo "Failed to call Rscript command. Please check your installation."
        exit 1
    fi

    echo "Producing consensus file across biological replicates..."
    # Find consensus accross all reps of a given sample group
    if [[ "$run_mod" = true ]]; then
        Rscript "$scripts"/findConsensus.R \
            "$out"/hamr_out \
            "$out"/hamr_consensus
    fi

    if [[ "$run_lnc" = true ]]; then
        Rscript "$scripts"/findConsensus_lnc.R \
            "$out"/lnc_out \
            "$out"/lnc_consensus
    fi

    wait
    echo "done"

    # The case where no consensus file is found, prevents *.bed from being created
    if [ -z "$(ls -A "$out"/hamr_consensus)" ]; then
    echo "No consensus mods found within any sequencing group. Please see check individual rep for analysis. "
    exit 1
    fi

    # Add depth columns with info from each rep alignment, mutate in place
    for f in "$out"/hamr_consensus/*.bed
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

    for f in "$out"/hamr_consensus/*.bed
    do
        if [ -s "$f" ]; then
        # The file is not-empty.
            t=$(basename "$f")
            n=${t%.*}
            echo "computing depth across reps for $n"
            Rscript "$scripts"/depthHelperAverage.R "$f"
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
    # Produce overlap bam files with the provided annotation library folder and store in ~/lap
    if [ ! -d "$out/lap" ]; then mkdir "$out"/lap; echo "created path: $out/lap"; fi

    # Run a series of command checks to ensure consensusOverlap can run smoothly
    if ! command -v intersectBed > /dev/null; then
        echo "Failed to call intersectBed command. Please check your installation."
        exit 1
    fi

    # annot bed target directory
    if [ ! -d "$out/annotBeds" ]; then mkdir "$out"/annotBeds; echo "created path: $out/annotBeds"; fi

    # checks if genomedir is populated with generated annotation files, if not, hamrbox can't run anymore, exit
    count=$(ls -1 "$out/annotBeds"/*.bed 2>/dev/null | wc -l)
    if [ "$count" == 0 ]; then 
        if [[ -e "$generator" ]]; then
            echo "generating annotations for overlap..."
            # 11/17 redirect annotation generate output to out/annotBeds, second arg added
            Rscript "$generator" "$annotation" "$out/annotBeds"
            echo ""
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
    for smp in "$out"/hamr_consensus/*
    do 
        consensusOverlap
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
    echo "############### Obtained all necessary pipeline files, generating dynamic epitranscriptomic summary ################"
    date '+%d/%m/%Y %H:%M:%S'
    echo ""
    #######################################begins EXTRACT######################################
    if [ ! -d "$out/results" ]; then mkdir "$out"/results; echo "created path: $out/results"; fi
    dir="$out/results"

    echo "generating long modification table..."
    # collapse all overlapped data into longdf
    Rscript "$scripts"/concatenate4R.R \
        "$out"/lap \
        "$out/results"
    echo "done"
    echo ""

    # note mod_long.csv is now in dir/results, update

    echo "plotting modification abundance per sample group..."
    # overview of modification proportion
    Rscript "$scripts"/countPerGroup.R \
        "$dir"/mod_long.csv \
        "$out"/annotBeds \
        "$dir"
    echo "done"
    echo ""

    echo "plotting modification abundance per mod type..."
    # overview of modification proportion
    Rscript "$scripts"/countPerMod.R \
        "$dir"/mod_long.csv \
        "$out"/annotBeds \
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

        echo "editing panther param file with user input..."
        # edit params/enrich.json with user's input
        cd $util
        if [[ ! -z $porg ]]; then
            mv panther_params.json temp.json
            jq --arg jq_in $porg -r '.organism |= $jq_in' temp.json > panther_params.json
            rm temp.json

            mv panther_params.json temp.json
            jq --arg jq_in $porg -r '.refOrganism |= $jq_in' temp.json > panther_params.json
            rm temp.json
        fi

        if [[ ! -z $pterm ]]; then
            mv panther_params.json temp.json
            jq --arg jq_in $pterm -r '.annotDataSet |= $jq_in' temp.json > panther_params.json
            rm temp.json
        fi

        if [[ ! -z $ptest ]]; then
            mv panther_params.json temp.json
            jq --arg jq_in $ptest -r '.enrichmentTestType |= $jq_in' temp.json > panther_params.json
            rm temp.json
        fi

        if [[ ! -z $pcorrect ]]; then
            mv panther_params.json temp.json
            jq --arg jq_in $pcorrect -r '.correction |= $jq_in' temp.json > panther_params.json
            rm temp.json
        fi
        cd

        # proceed if genelists directory is not empty
        if [ -n "$(ls "$dir"/go/genelists)" ]; then
            echo "sending each gene list to panther for overrepresentation analysis..."
            # Send each gene list into panther API and generate a overrepresentation result file in another folter
            for f in "$dir"/go/genelists/*.txt
            do
                n=$(basename "$f")
                echo "$n"
                python $execpthr \
                    --service enrich \
                    --params_file $json \
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
    echo "done"
    echo ""

    echo "classifying modified RNA subtype..."
    # looking at RNA subtype for mods
    Rscript "$scripts"/RNAtype.R \
        "$dir"/mod_long.csv
    echo "done"
    echo ""

    if [ -e "$out"/annotBeds/*_CDS.bed ] && [ -e "$out"/annotBeds/*_fiveUTR.bed ] && [ -e "$out"/annotBeds/*_threeUTR.bed ]; then
        c=$(find "$out"/annotBeds -type f -name "*_CDS.bed")
        f=$(find "$out"/annotBeds -type f -name "*_fiveUTR.bed")
        t=$(find "$out"/annotBeds -type f -name "*_threeUTR.bed")
        echo "mapping modification regional distribution landscape..."
        # improved region mapping
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

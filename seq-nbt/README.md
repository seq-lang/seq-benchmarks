# Seq Experimental Notebook

The following datafiles were used to conduct the experiments:

- Human genome v37 (hs37d5): 
  - [`genome.fa` and `chr1.fa`](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz)
  - [`Homo_sapiens_assembly19_1000genomes_decoy.fa`](gs://gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.fasta)
- BWA-MEM FASTQs: 
  - [`sample.fastq`](http://cb.csail.mit.edu/cb/seq/nbt/bwa-data.bz2)
- mrsFAST FASTQs: 
  - [`FIN1_1.fastq`](http://cb.csail.mit.edu/cb/seq/nbt/mrsfast-data.bz2)
- Smith-Waterman sequences: 
  - [`queries*`, `targets*`](http://cb.csail.mit.edu/cb/seq/nbt/sw-data.tar.bz2)
- Segmental duplications for AVID: 
  - [`wgac-data` and `sedef-data`](http://cb.csail.mit.edu/cb/seq/nbt/avid-data.tar.bz2)
- UMI-tools FASTQs: 
  - [`hgmm_100_R1.fastq`](http://cb.csail.mit.edu/cb/seq/nbt/umi-data.bz2)
- GATK BAMs and HapTree-X BAMs:
  - Details are available [here](https://github.com/0xTCG/haptreex/blob/master/paper/experiments.ipynb). 
    For 10X sample, we used NA12878 sample BAM restricted to chromosomes 19–22 (inclusive).

The following commands will set up the benchmarking environment:
```bash
# Location of tools' sources and executables
export TOOLS=`pwd`/tools
# Location of test datasets
export DATA=`pwd`/data
# Locatiion of the local Seq installation
export SEQ=`pwd`/seq
# Make sure to point LD_LIBRARY_PATH, PATH and SEQ_LIBRARY_PATH to your Seq installation directory
export LD_LIBRARY_PATH=$SEQ/build:$SEQ/deps/lib:${LD_LIBRARY_PATH}
export PATH=$SEQ/build:${PATH}
export SEQ_LIBRARY_PATH=$SEQ/build
# Use 1 thread
export OMP_NUM_THREADS=1
# Current working directory
DIR=`pwd`
# Make the necessary directories
mkdir -p $DIR/runs/time
mkdir -p $DIR/runs/stderr
mkdir -p $DIR/runs/stdout
mkdir -p $TOOLS/seq
export OUT=$DIR/runs/out
mkdir -p $OUT
```

The following function is used to collect the performance metrics of each experiment.
```bash
function tm {
  TIME=`which time`  # use gtime on macOS
  NAME=$1
  if grep -qs 'Exit status: 0' $DIR/runs/time/$NAME ; then
    echo -ne "- $NAME already done "
    grep  'Exit status:' $DIR/runs/time/$NAME
  else
    cmd="$TIME -v -o $DIR/runs/time/$NAME ${@:2} >$DIR/runs/stdout/$NAME 2>$DIR/runs/stderr/$NAME"
    eval "$cmd"
    es=$?
    if [ $es -ne 0 ] ; then
      echo "- $NAME failed with $es: $cmd"
      cat $DIR/runs/stderr/$NAME
    else
      echo "- $NAME completed"
    fi
  fi
}
```

## CORA (Homology table construction)

We have used [CORA 1.1.5b](https://github.com/denizy/cora/tree/6a53c7b9b25eb8c33d3de92f00481310c16fdd5c) to evaluate CORA performance.

CORA first needs to be run with default parameters to generate `__temporary_CORA_files` directory.
You can kill CORA once this directory is generated and proceed by running the commands below:

```bash
cd $OUT/cora
# Generate inexact homology table
tm cora-inexact $TOOLS/cora/homTable_setup $DATA/genome.fa 64 __temporary_CORA_files/__TEMP__Aux_Physical_Splits_File cora_hom_exact_k64 cora_hom_inexact_k64 __temporary_CORA_files/__TEMP__Aux_Parallel_Splits_File EXACT 1 FULL
# Generate exact homology table (needs 2 steps)
tm cora-exact-1 $TOOLS/cora/homTable_setup $DATA/genome.fa 64 __temporary_CORA_files/__TEMP__Aux_Physical_Splits_File cora_hom_exact_k64 cora_hom_inexact_k64 __temporary_CORA_files/__TEMP__Aux_Parallel_Splits_File EXACT 1 ONLY_EXACT_COMPACT
tm cora-exact-2 $TOOLS/cora/homTable_setup $DATA/genome.fa 64 __temporary_CORA_files/__TEMP__Aux_Physical_Splits_File cora_hom_exact_k64 cora_hom_inexact_k64 __temporary_CORA_files/__TEMP__Aux_Parallel_Splits_File BOTH  1 ONLY_CONSTRUCT
cd ..
```

For Seq, we used:
```bash
# Build
tm seq-build-cora-exact   seqc build -release -o $TOOLS/seq/cora-exact   $TOOLS/seq/cora/hom_exact.seq
tm seq-build-cora-inexact seqc build -release -o $TOOLS/seq/cora-inexact $TOOLS/seq/cora/hom_inexact.seq

# Use 24 threads
export OMP_NUM_THREADS=24
# Generate exact homology table
tm seq-cora-exact $TOOLS/seq/cora-exact   $DATA/genome.fa   $OUT/seq-cora-k64
# Generate inexact homology table
tm seq-cora-inexact $TOOLS/seq/cora-inexact $DATA/genome.fa 1 $OUT/seq-cora-k64
export OMP_NUM_THREADS=1
```

Seq reported 438k homologies that were not found by CORA; CORA, on the other hand, reported 252k homologies that were not found by Seq.
These homologies are available for inspection [here](http://cb.csail.mit.edu/cb/seq/nbt/cora-homologies.tbz).
After inspection, we discovered that the homologies not found by Seq are not true homologies due to a bug in the original version
(for example, a reported homology between `chr6:78300426-78300505` and `chr6:78300466-78300545` is incorrect as can be seen by comparing [the first](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr6%3A78300426%2D78300505&hgsid=1079285129_4USpZq0gDfCkCIMWcmtMyuKdrAMJ) and [the second](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr6%3A78300466%2D78300545&hgsid=1079285129_4USpZq0gDfCkCIMWcmtMyuKdrAMJ) region— note the sequence differences at the end).

## BWA-MEM (fastmap)

We used [BWA 0.7.17](https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2) with the [following patch](tools/bwa-0.7.17.patch) that reports the time needed to load an index. 

BWA was run as follows:
```bash
# Prepare the index
tm bwa-index-chr1 $TOOLS/bwa-0.7.17/bwa index   $DATA/chr1.fa
tm bwa-chr1       $TOOLS/bwa-0.7.17/bwa fastmap $DATA/chr1.fa $DATA/sample.fastq
```

Seq was run as follows:
```bash
# Build prefetch and no-prefetch versions of fastmap
tm seq-build-fastmap-bn seqc build -release -o $TOOLS/seq/fastmap-bn $TOOLS/seq/bwa/fastmap_build.seq
tm seq-build-fastmap-bp "sed 's/^#\@prefetch/\@prefetch/' $TOOLS/seq/bwa/fastmap_build.seq | seqc build -release -o $TOOLS/seq/fastmap-bp -"
# Run both versions
tm seq-fastmap-bn-chr1 $TOOLS/seq/fastmap-bn search $DATA/chr1.fa $DATA/sample.fastq $OUT/seq-fastmap-bn
tm seq-fastmap-bp-chr1 $TOOLS/seq/fastmap-bp search $DATA/chr1.fa $DATA/sample.fastq $OUT/seq-fastmap-bp
```


Rust code was built with [rust-bio v0.32.0](https://github.com/rust-bio/rust-bio/tree/v0.32.0) with the following [patch](tools/rust/rust-bio.patch) that allows proper SMEM reconstruction:
```bash
# Build the tool
(cd $TOOLS/rust; cargo build --release)
# Run Rust fastmap
tm rust-fastmap-chr1 $TOOLS/rust/target/release/biorust fastmap $DATA/chr1.fa $DATA/sample.fastq $OUT/rust-fastmap-chr1
```

SeqAn/C++ code was built with [SeqAn v3.0.2](https://github.com/seqan/seqan3/releases/download/3.0.2/seqan3-3.0.2-Source.tar.xz) with g++ 10:
```bash
# Build the tool
( cd $TOOLS/seqan; 
  g++ -std=c++2a -fconcepts -O3 smem_seqan.cc -o fastmap \
    -I seqan/3.0.2/include -I seqan/3.0.2/include/seqan3/submodules/range-v3/include -I seqan/3.0.2/include/seqan3/submodules/sdsl-lite/include )
# Run SeqAn fastmap
tm seqan-fastmap-chr1 $TOOLS/seqan/fastmap $DATA/chr1.fa $DATA/sample.fastq
```


## mrsFAST

We used [mrsFAST v3.4.1](https://github.com/sfu-compbio/mrsfast/tree/75b49304992302ed8754da139ccb464297fda930) with the following [patch](tools/mrsfast-3.4.1.patch)
that turns off the calculation of MD tag for fair comparison, and that force-enabes SSE4 version to be built.

```bash
# Index the genomes
tm mrsfast-index-hg19 $TOOLS/mrsfast/mrsfast --index  $DATA/genome.fa
tm mrsfast-index-chr1 $TOOLS/mrsfast/mrsfast --index  $DATA/chr1.fa
# Inexact version (e=2)
tm mrsfast-chr1 $TOOLS/mrsfast/mrsfast --search $DATA/chr1.fa --seq $DATA/FIN1_1.fastq --crop 96 -e 2 -o $OUT/mrsfast-chr1 -u /dev/null --disable-sam-header --mem 10
# Exact version
tm mrsfast-hg19-all $TOOLS/mrsfast/mrsfast --search $DATA/genome.fa --seq $DATA/FIN1_1.fastq -e 0 -o $OUT/mrsfast-hg19-all -u /dev/null --disable-sam-header --mem 200
```

Seq version was run as follows:

```bash
# Build prefetch and no-prefetch versions of exact match mapper (FM-index)
tm seq-build-mrsfast-exact-n seqc build -release -o $TOOLS/seq/mrsfast-exact-n $TOOLS/seq/mrsfast/exact.seq
tm seq-build-mrsfast-exact-p "sed 's/^#@prefetch/@prefetch/' $TOOLS/seq/mrsfast/exact.seq | seqc build -release -o $TOOLS/seq/mrsfast-exact-p -"
# Build inexact mapper (k-mer index)
tm seq-build-mrsfast-inexact seqc build -release -o $TOOLS/seq/mrsfast-inexact $TOOLS/seq/mrsfast/mrsfast.seq
# Index the genomes
tm seq-mrsfast-index-hg19 $TOOLS/seq/mrsfast-exact-n index  $DATA/genome.fa
tm seq-mrsfast-index-chr1 $TOOLS/seq/mrsfast-inexact index  $DATA/chr1.fa
# Inexact version
tm seq-mrsfast-chr1 $TOOLS/seq/mrsfast-inexact search $DATA/chr1.fa $DATA/FIN1_1.fastq $OUT/seq-mrsfast-chr1
# Exact version
tm seq-mrsfast-hg19-n $TOOLS/seq/mrsfast-exact-n search $DATA/genome.fa $DATA/FIN1_1.fastq $OUT/seq-mrsfast-hg19-n
tm seq-mrsfast-hg19-p $TOOLS/seq/mrsfast-exact-p search $DATA/genome.fa $DATA/FIN1_1.fastq $OUT/seq-mrsfast-hg19-p
```

## minimap2 (Smith-Waterman alignment)

The original [minimap2](https://github.com/lh3/minimap2) uses [KSW2](https://github.com/lh3/ksw2) library that Seq uses internally for `seq.align` function.
Thus, we just run unoptimized `seq.align` to measure the performance of minimap2 / KSW2 alignment.

```bash
# Build the tool
tm seq-build-sw seqc build -release -o $TOOLS/seq/sw $TOOLS/seq/minimap2/sw.seq
# Run the benchmark. It reports the average time for: 
#  - non-optimized alignment, 
#  - SSE4-optimized inter-alignment, 
#  - AVX2-optimized inter-alignment, and 
#  - AVX512F-optimized inter-alignment (if available).
tm seq-sw $TOOLS/seq/sw $DATA/queries_m $DATA/targets_m
```

We compared our alignment with [rust-bio v0.32.0](https://github.com/rust-bio/rust-bio/tree/v0.32.0), [SeqAn v3.0.2](https://github.com/seqan/seqan3/releases/download/3.0.2/seqan3-3.0.2-Source.tar.xz) and [Bio.jl and BioAlignments.jl v1.0.1](https://biojulia.net/Bio.jl/stable) in Julia v1.6.0.
Other tools (SeqAn, BioJulia and rust-bio) do not support consistently advanced alignment options (e.g. setting a band width, Z-drop score and so on).
Thus, we compared a simple version of Seq's `seq.align` to these tools to maintain the consistency:

```bash
# Build and run a simple version of the Seq alignment benchmark
tm seq-build-sw-simple seqc build -release -o $TOOLS/seq/sw-simple $TOOLS/seq/minimap2/sw_simple.seq
tm seq-sw-simple $TOOLS/seq/sw-simple $DATA/queries.small $DATA/targets.small
# Build and run rust-bio alignment benchmark
(cd $TOOLS/rust; cargo build --release)
tm rust-sw-simple $TOOLS/rust/target/release/biorust align $DATA/queries.small $DATA/targets.small
# Build and run SeqAn v3 alignment benchmark
( cd $TOOLS/seqan; 
  g++ -std=c++2a -fconcepts -O3 align_seqan.cc -o sw \
    -I seqan/3.0.2/include -I seqan/3.0.2/include/seqan3/submodules/range-v3/include -I seqan/3.0.2/include/seqan3/submodules/sdsl-lite/include )
tm seqan-sw-simple $TOOLS/seqan/sw $DATA/queries.small $DATA/targets.small
# Run Julia/BioJulia alignment benchmark
tm julia-sw-simple julia --check-bounds=no -O3 $TOOLS/julia/align.jl $DATA/queries.small $DATA/targets.small
```

## AVID (global alignment)

We used a Linux version of [AVID v2.1 build 0](http://genome.lbl.gov/vista/mvista/download.shtml). (Note: you need to register to download a statically compiled binary. Source code is not available).

As AVID only supports reading a single pair of sequences, we used a wrapper script [runavid.py](tools/runavid.py) to run AVID and its Seq counterpart on many pairs of sequences and collect the correct time.

```bash
# Run AVID on WGAC and SEDEF-generated segmental duplications
tm avid-wgac  python3 $TOOLS/runavid.py $DATA/wgac-segdups  $OUT/avid-wgac  avid
tm avid-sedef python3 $TOOLS/runavid.py $DATA/sedef-segdups $OUT/avid-sedef avid
```

Seq was run as follows:
```bash
# Build Seq-AVID
tm seq-build-sw seqc build -release -o $TOOLS/seq/sw $TOOLS/seq/avid/avid.seq
# Run Seq-AVID on WGAC and SEDEF-generated segmental duplications
tm seq-avid-wgac  python3 $TOOLS/runavid.py $DATA/wgac-segdups  $OUT/seq-avid-wgac  seq
tm seq-avid-sedef python3 $TOOLS/runavid.py $DATA/sedef-segdups $OUT/seq-avid-sedef seq
```

## UMI-tools (barcode whitelisting)

We used [UMI-tools v1.1.1](https://github.com/CGATOxford/UMI-tools/tree/1.1.1):

```bash
tm umi-tools umi_tools whitelist --stdin $DATA/hgmm_100_R1.fastq --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN --log2stderr
```

Seq verson was run as follows:
```bash
# Build
tm seq-build-umi seqc build -release -o $TOOLS/seq/umi $TOOLS/seq/umi-tools/whitelist.seq
# Run
tm seq-umi $TOOLS/seq/umi $DATA/hgmm_100_R1.fastq
```

## GATK (SplitNCigar)

To ensure fairness in comparison, all files were converted to single-ended BAM files with no optional tags:
```bash
for i in giab k562 cytosol; do 
  ( samtools view $DATA/$i.bam -H; \
    samtools view $DATA/$i.bam \
    | awk '{$1="read"NR; $2=and($2,compl(235)); $7="*"; $8=0; $9=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' 'OFS=\t' ) \
  | samtools view - -1 -o $DATA/$i.gatk.bam -@ 4 
done
```

We used [GATK 4.1.4.1](https://github.com/broadinstitute/gatk/releases/download/4.1.4.1/gatk-4.1.4.1.zip) to run `SplitNCigar` as follows:
```bash
for i in giab k562 cytosol; do
  tm gatk-$i $TOOLS/gatk-4.1.4.1/gatk SplitNCigarReads -I $DATA/$i.gatk.bam -O $OUT/gatk-$i.sam -R $DATA/Homo_sapiens_assembly19_1000genomes_decoy.fasta
done
```

Seq was run as follows:
```bash
tm seq-build-gatk seqc build -release -o $TOOLS/seq/gatk $TOOLS/seq/gatk-splitncigar/gatk-split.seq
for i in giab k562 cytosol; do
  tm seq-gatk-$i $TOOLS/seq/gatk $DATA/$i.gatk.bam $DATA/Homo_sapiens_assembly19_1000genomes_decoy.fasta $OUT/gatk-$i.sam
done
```

## HapTree-X (Haplotype phasing)

We compared [HapTree-X fd3bef](https://github.com/0xTCG/haptreex/tree/fd3befbebfca0c2134bf8a642e5294dca6af5e8e) to [HapCUT2 v1.2](https://github.com/vibansal/HapCUT2/tree/de33b57106898fc46751d55bb338456dfd5dbbe5) and [phASER v1.1.1](https://github.com/secastel/phaser/tree/20eb3b44394fe0f2f3e002ffbcf7b0db6ff01966).

HapCUT2 was run as follows:
```bash
# RNA-seq and exome data:
for i in na12878 exome ; do
  # - Extract fragments
  tm hapcut2-extract-$i $TOOLS/hapcut2/build/extractHAIRS --bam $DATA/$i.bam --VCF $DATA/$i.vcf --out $OUT/${i}-frag-hapcut2
  # - Run phasing
  tm hapcut2-$i $TOOLS/hapcut2/build/HAPCUT2 --fragments $OUT/${i}-frag-hapcut2 --VCF $DATA/$i.vcf --output $OUT/${i}-hapcut2
done

# 10X data:
# - Extract fragments
tm hapcut2-extract-10x $TOOLS/hapcut2/build/extractHAIRS --bam $DATA/10x.bam --VCF $DATA/10x.vcf --out $OUT/10x-frag-hapcut2 --10X 1
# - Link barcode fragments
tm hapcut2-link-10x python3 $TOOLS/hapcut2/utilities/LinkFragments.py --bam $DATA/10x.bam --VCF $DATA/10x.vcf --fragments $OUT/10x-frag-hapcut2 --out $OUT/10x-link-hapcut2
# - Run phasing
tm hapcut2-10x $TOOLS/hapcut2/build/HAPCUT2 --nf 1 --VCF $DATA/10x.vcf --fragments $OUT/10x-link-hapcut2 --output $OUT/10x-hapcut2
```

phASER was run as follows:
```bash
tm phaser python2 $TOOLS/phaser/phaser/phaser.py --vcf $DATA/na12878.vcf.gz --bam $DATA/na12878.bam --paired_end 1 --mapq 255 --baseq 0 --sample NA12878 --o $DATA/na12878.phaser
```

phASER was run only on RNA-seq sample as it is only designed to phase RNA-seq data.

HapTree-X was run as follows:
```
# Build
tm seq-build-haptreex seqc build -release -o $TOOLS/seq/haptreex $TOOLS/haptreex/src/main.seq
# Run
tm haptreex-rna   $TOOLS/seq/haptreex -v $DATA/na12878.vcf -r $DATA/na12878.bam -o $OUT/na12878.haptreex -g $DATA/Homo_sapiens.GRCh37.87.gtf
tm haptreex-exome $TOOLS/seq/haptreex -v $DATA/exome.vcf -d $DATA/exome.bam -o $OUT/exome.haptreex
tm haptreex-10x   $TOOLS/seq/haptreex -v $DATA/10x.vcf -d $DATA/10x.bam -o $OUT/10x.haptreex --10x

# Run 10X phasing with 4 threads
export OMP_NUM_THREADS=4
tm haptreex-10x-4 $TOOLS/seq/haptreex -v $DATA/10x.vcf -d $DATA/10x.bam -o $OUT/10x.haptreex-4 --10x
```

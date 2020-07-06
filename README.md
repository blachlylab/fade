![alt text](https://github.com/blachlylab/fade/raw/master/logo/fade_logo.png "FADE")

# **F**ragmentase **A**rtifact **D**etection and **E**limination

DNA shearing is a crucial first step in most NGS protocols for Illumina. Enzymatic fragmentation has shown in recent years to be a cost and time effective alternative to physical shearing (i.e. sonication). We discovered that enzymatic fragmentation leads to unexpected alteration of the original DNA source material. We provide fade as a method of identification and removal of enymatic fragmentation artifacts.

Our [documentation](https://github.com/blachlylab/fade/blob/master/INSTALL.md) has information on installing fade and its prerequisites.

### Running FADE 
```
fade annotate -b sam1.bam ref.fa > sam1.anno.bam
samtools sort -n sam1.anno.bam > sam1.anno.qsort.bam #recommended but not neccessary
fade out -b sam1.anno.qsort.bam > sam1.filtered.bam
```

### Running FADE via Docker
```
docker run -v `pwd`:/data charesgregory/fade annotate -b /data/sam1.bam /data/ref.fa > sam1.anno.bam
docker run -v `pwd`:/data charesgregory/fade out -b -c /data/sam1.anno.qsort.bam > sam1.filtered.bam
```

## Program Details

### Command-line parameters

fade
```
Fragmentase Artifact Detection and Elimination
usage: ./fade [subcommand]
        annotate: marks artifact reads in bam tags (must be done first)
        out: eliminates artifact from reads(may require queryname sorted bam)
-h --help This help information.
```
fade annotate
```
Fragmentase Artifact Detection and Elimination
annotate: performs re-alignment of soft-clips and annotates bam records with bitflag (rs) and realignment tags (am)
usage: ./fade annotate [BAM/SAM input]

-t      --threads extra threads for parsing the bam file
     --min-length Minimum number of bases for a soft-clip to be considered for artifact detection
   --short-length Minimum number of bases for a soft-clip to not be considered short
-w  --window-size Number of bases considered outside of read or mate region for re-alignment
-q      --q-score Minimum average base-quality score of a soft-clip to be considered an artifact
-m     --mate-est Read Mate size estimate in bases
-b          --bam output bam
-u         --ubam output uncompressed bam
-h         --help This help information.
```
fade out
```
Fragmentase Artifact Detection and Elimination
out: removes all read and mates for reads contain the artifact (used after annotate and requires queryname sorted bam) or, with the -c flag, hard clips out artifact sequence from reads
usage: ./fade out [BAM/SAM input]

-c         --clip clip reads instead of filtering them
-t      --threads extra threads for parsing the bam file
-a --artifact-bam filename to extract artifact reads to (BAM/SAM)
-b          --bam output bam
-u         --ubam output uncompressed bam
-h         --help This help information.
```

## Algorithm
FADE is written in D and uses the [htslib](http://www.htslib.org/download/) library via [dhtslib](https://github.com/blachlylab/dhtslib.git), and the [parasail](https://github.com/jeffdaily/parasail) library via [dparasail](https://github.com/blachlylab/dparasail). FADE accepts SAM/BAM/CRAM files containing reads that have been mapped to a reference genome and filters or cleans up artifact-containing reads according to the following procedure. 

FADE is designed to determine a sequencing read’s enzymatic artifact status by employing aligner soft-clipping. Soft-clipping is an action performed by the aligner to improve the alignment score of a read to the reference by ignoring a portion on one end of the read. Soft-clipping can help an aligner correctly align a read that has sequencing error on one end of the read or has adapter contamination. FADE employs soft-clipping to identify potentially enzymatic artifact containing reads. 
1. It considers only those reads aligned with soft-clipping; alignments without soft-clips. 
2. The reference sequence that the alignment is mapped to is extracted such that there exists 300 nucleotides (nt) of padding on each end of the mapped read.\* 
3. The read is reverse-complemented and then aligned via a Smith-Waterman local alignment to the extracted region of reference sequence. We use a scoring matrix with a gap open penalty of 10, a gap extension penalty of 2, a mismatch penalty of 3, and a match score of 2.\*\* 

FADE makes available three subcommands that all rely on the algorithm described above. The **annotate** subcommand performs the initial analysis and adds BAM tags encoding information concerning artifact status to the alignments, used during filtration to remove the artifacts. The **out** subcommand eliminates the artifact. It either removes reads from the output BAM/SAM file completely if they or their mate contain an identified fragmentation artifact or artifact-containing reads are trimmed to remove extraneous sequence originating from the opposite strand. After elimination, FADE reports statistics describing the total number of alignments, the percentage of soft-clipped alignments, and the percentage of enzymatic artifacts found. The filtration step must be run on a queryname-sorted BAM file in order to fully filter out the read, its mate, and any other supplementary or secondary alignments. (Stats description)

<sub><sup>\* The 300 nt padding on each end of the mapped region provides ample search space for artifact alignment search without being too computationally expensive; most artifacts originate very close to the mapped region and 300 nt was chosen as an optimal tradeoff, but could be adjusted.</sub></sup>
<br/><sub><sup>\*\* Harsher gap penalties allows the algorithm to be strict in allowing gaps, since we expect the artifact sequences to directly match the reference, except for soft-clipped regions derived from sequencing error. A soft-clipped region is considered to be an artifact if there is a 90% or greater match to the opposite strand sequence. </sub></sup>

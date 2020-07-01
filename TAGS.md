### BAM tags
FADE uses BAM/SAM tags to annotate artifact reads. These tags are then used to identify and remove/clip reads.

#### rs tag
The rs tag is a 6 bit flag that indicates read artifact status.

| Bit | Description                         |
|-----|-------------------------------------|
|0    |Read is Softclipped                  | 
|1    |Read's left sc is an artifact        |
|2    |Read's right sc is an artifact       |
|3    |Left Artifact aligns to mate region  |
|4    |Right Artifact aligns to mate region |
|5    |Read has supplementary alignment     |

#### am tag
Contains mapping information of the artifact sequence as a string:
Format is CHROM, POS, CIGAR with left and right artifacts delimited by ;
i.e.   chr1,20,120S20M;chr1,120,15M120S

#### as, ap, and ar tags
As described in our manuscript we describe an expected stem loop structure:


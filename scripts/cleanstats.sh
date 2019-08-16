awk 'BEGIN{OFS="\t"}{print $2,$5,$6,$1}' $1 > read_targets.bed
awk 'BEGIN{OFS="\t"}{print $7,$8,$9,$1}' $1 > artifact_targets.bed
sed -i '1d' read_targets.bed
sed -i '1d' artifact_targets.bed
bedtools sort -i read_targets.bed > read_targets.sorted.bed
bedtools sort -i artifact_targets.bed > artifact_targets.sorted.bed
bedtools subtract -A -a read_targets.sorted.bed -b hg19.rpmasker.sorted.bed > read_targets.nrp.bed
bedtools subtract -A -a artifact_targets.sorted.bed -b hg19.rpmasker.sorted.bed > artifact_targets.nrp.bed
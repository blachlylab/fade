wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.out.gz
gzip -d hg19.fa.out.gz
awk 'BEGIN{OFS="\t"}{if(NR>3) {if($9=="C"){strand="-"}else{strand="+"};print $5,$6-1,$7,$10,".",strand}}' hg19.fa.out > hg19.rpmasker.bed
bedtools sort -i hg19.rpmasker.bed > hg19.rpmasker.sorted.bed

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz
gzip -d hg38.fa.out.gz
awk 'BEGIN{OFS="\t"}{if(NR>3) {if($9=="C"){strand="-"}else{strand="+"};print $5,$6-1,$7,$10,".",strand}}' hg38.fa.out > hg38.rpmasker.bed
bedtools sort -i hg38.rpmasker.bed > hg38.rpmasker.sorted.bed

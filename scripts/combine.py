import pandas as pd

main=pd.read_csv("test.out",delimiter="\t")
other=pd.read_csv("read_targets.nrp.bed",delimiter="\t",header=None,names=["rname","art_start","art_end","qname"])
main.set_index(["qname","rname","art_start"],inplace=True)
other.set_index(["qname","rname","art_start"],inplace=True)
main["read_rpm"]=1
main["art_rpm"]=1
main.loc[other.index,"read_rpm"]=0

other=pd.read_csv("artifact_targets.nrp.bed",delimiter="\t",header=None,names=["aln_rname","aln_start","aln_end","qname"])
main.reset_index(inplace=True)
main.set_index(["qname","aln_rname","aln_start"],inplace=True)
other.set_index(["qname","aln_rname","aln_start"],inplace=True)
main.loc[other.index,"art_rpm"]=0
main.reset_index(inplace=True)

main.to_csv("master.tsv",sep="\t",index=False)

main[["aln_rname","aln_start","aln_end"]].to_csv("master.bed",sep="\t",index=False,header=False)
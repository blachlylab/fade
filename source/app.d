import std.stdio;
import std.algorithm.iteration:splitter;
import std.algorithm:map,reverse;
import std.regex:splitter,regex;
import std.range:drop,array;
import std.conv:to;
import std.traits:ReturnType;
import bio.bam.reader;
import bio.bam.cigar;
import bio.bam.writer;
import std.algorithm.iteration:filter;
import std.algorithm:count;
import dparasail;
import dhtslib.faidx:IndexedFastaFile;

string rc(Range)(Range seq){
	seq.reverse;
	return seq.map!(x=>cast(char)x.complement).array.idup;
}

unittest{
	auto seq="ACGATGATCGATNGT".dup;
	writeln(rc(seq));
}



void main(string[] args){
	auto bam = new BamReader(args[1]);
	auto fai=IndexedFastaFile(args[2]);
	auto sc_bam= new BamWriter("sc.bam");
	sc_bam.writeSamHeader(bam.header());
	sc_bam.writeReferenceSequenceInfo(bam.reference_sequences());
	//auto db_bam=new BamWriter("db.bam");
	//db_bam.writeSamHeader(bam.header());
	//db_bam.writeReferenceSequenceInfo(bam.reference_sequences());
	auto art_bam=new BamWriter("art.bam");
	art_bam.writeSamHeader(bam.header());
	art_bam.writeReferenceSequenceInfo(bam.reference_sequences());
	//auto non_art_bam=new BamWriter("non_art.bam");
	//non_art_bam.writeSamHeader(bam.header());
	//non_art_bam.writeReferenceSequenceInfo(bam.reference_sequences());
	auto out_bam=new BamWriter("out.bam");
	out_bam.writeSamHeader(bam.header());
	out_bam.writeReferenceSequenceInfo(bam.reference_sequences());
	//string[2] strands=["+","-"];

	ubyte[string] reads;
	/+
	ubyte read status encoding
	0	Read is Normal
	1	Read is Softclipped
	2	Read has Supp Alignment
	3   Read is Artifact
	4	Artifact aligns to read region
	5	Artifact aligns to mate region
	6   Mate is on Diff Chrom than read
	7	5' or 3' artifact
	+/
	auto p=Parasail("ACTGN",1,-1,1,3);
	foreach(BamRead rec;bam.allReads()){
		ubyte read_class=1;
		if(rec.mate_ref_name()!=rec.ref_name())
			read_class|=0b100_0000;
		if(rec.is_supplementary()||rec.is_secondary_alignment()||rec.is_unmapped()){
			//reads[rec.name]=read_class;
			continue;
		}
		if (rec.cigar.filter!(x=>x.type=='S').count()==0){
			//reads[rec.name]=read_class;
			continue;
		}
		read_class|=0b10;
		if(!rec["SA"].is_nothing()){
			read_class|=0b100;
			string[] sup=rec["SA"].toString.splitter(",").array;
			if(sup[0]==rec.ref_name()){
				if (
					(sup[1].to!int>rec.position()-300)&&
					(sup[1].to!int<rec.position()+300)//&&
					//(strands[rec.is_reverse_strand]!=sup[2])
				){
					read_class|=0b1000;
					read_class|=0b1_0000;
					reads[rec.name]=read_class;
					continue;
				}
			}else if(sup[0]==rec.mate_ref_name()){
				if (
					(sup[1].to!int>rec.mate_position()-300)&&
					(sup[1].to!int<rec.mate_position()+300)//&&
					//(strands[rec.mate_is_reverse_strand]!=sup[2])
				){
					read_class|=0b1000;
					read_class|=0b10_0000;
					reads[rec.name]=read_class;
					continue;
				}
			}
		}
		CigarOperation[2] clips;
		bool first=true;
		foreach(CigarOperation op;rec.cigar){
			auto is_sc=op.is_query_consuming()&&op.is_clipping();
			if(first&&!is_sc){
				first=false;
			}else if(first&&is_sc){
				clips[0]=op;
			}else if(is_sc){
				clips[1]=op;
			}
		}
		//SW parasail needed
		string q_seq;
		string ref_seq;
		float cutoff;
		int start,end,score_read,score_mate;
		parasail_query res;
		//left soft-clip (left on reference not 5' neccesarily)
		if((clips[0].length!=0)&&(clips[1].length()==0)){
			q_seq=rec.sequence[0..clips[0].length()].rc;
			cutoff=q_seq.length*0.75;
			start=rec.position()-300;
			if(start<0){
				start=0;
			}
			end=rec.position()+rec.basesCovered()+300;
			if(end>bam.header.getSequence(rec.ref_id).length){
				end=bam.header.getSequence(rec.ref_id).length;
			}
			ref_seq=fai.fetchSequence(rec.ref_name(),start,end);
			res=p.sw_striped(q_seq,ref_seq);
			score_read=res.result.score;
			res.close();
			start=rec.mate_position()-300;
			if(start<0){
				start=0;
			}
			end=rec.mate_position()+300;
			if(end>bam.header.getSequence(rec.mate_ref_id).length){
				end=bam.header.getSequence(rec.mate_ref_id).length;
			}
			ref_seq=fai.fetchSequence(rec.ref_name(),start,end);
			res=p.sw_striped(q_seq,ref_seq);
			score_mate=res.result.score;
			res.close();
			if(score_read>cutoff||score_mate>cutoff){
				if(score_read>score_mate){
					read_class|=0b1000;
					read_class|=0b1_0000;
				}else{
					read_class|=0b1000;
					read_class|=0b10_0000;
				}
			}
		}
		//right soft-clip
		else if((clips[0].length==0)&&(clips[1].length()!=0)){
			q_seq=rec.sequence[0..clips[1].length()].rc;
			cutoff=q_seq.length*0.75;
			start=rec.position()-300;
			if(start<0){
				start=0;
			}
			end=rec.position()+rec.basesCovered()+300;
			if(end>bam.header.getSequence(rec.ref_id).length){
				end=bam.header.getSequence(rec.ref_id).length;
			}
			ref_seq=fai.fetchSequence(rec.ref_name(),start,end);
			res=p.sw_striped(q_seq,ref_seq);
			score_read=res.result.score;
			res.close();
			start=rec.mate_position()-300;
			if(start<0){
				start=0;
			}
			end=rec.mate_position()+300;
			if(end>bam.header.getSequence(rec.mate_ref_id).length){
				end=bam.header.getSequence(rec.mate_ref_id).length;
			}
			ref_seq=fai.fetchSequence(rec.ref_name(),start,end);
			res=p.sw_striped(q_seq,ref_seq);
			score_mate=res.result.score;
			res.close();
			if(score_read>cutoff||score_mate>cutoff){
				if(score_read>score_mate){
					read_class|=0b1000;
					read_class|=0b1_0000;
				}else{
					read_class|=0b1000;
					read_class|=0b10_0000;
				}
			}
		}
		//both
		else if((clips[0].length!=0)&&(clips[1].length()!=0)){
			q_seq=rec.sequence[0..clips[0].length()].rc;
			cutoff=q_seq.length*0.75;
			start=rec.position()-300;
			if(start<0){
				start=0;
			}
			end=rec.position()+rec.basesCovered()+300;
			if(end>bam.header.getSequence(rec.ref_id).length){
				end=bam.header.getSequence(rec.ref_id).length;
			}
			ref_seq=fai.fetchSequence(rec.ref_name(),start,end);
			res=p.sw_striped(q_seq,ref_seq);
			score_read=res.result.score;
			res.close();
			start=rec.mate_position()-300;
			if(start<0){
				start=0;
			}
			end=rec.mate_position()+300;
			if(end>bam.header.getSequence(rec.mate_ref_id).length){
				end=bam.header.getSequence(rec.mate_ref_id).length;
			}
			ref_seq=fai.fetchSequence(rec.ref_name(),start,end);
			res=p.sw_striped(q_seq,ref_seq);
			score_mate=res.result.score;
			res.close();
			if(score_read>cutoff||score_mate>cutoff){
				if(score_read>score_mate){
					read_class|=0b1000;
					read_class|=0b1_0000;
				}else{
					read_class|=0b1000;
					read_class|=0b10_0000;
				}
			}
			q_seq=rec.sequence[0..clips[1].length()].rc;
			cutoff=q_seq.length*0.75;
			start=rec.position()-300;
			if(start<0){
				start=0;
			}
			end=rec.position()+rec.basesCovered()+300;
			if(end>bam.header.getSequence(rec.ref_id).length){
				end=bam.header.getSequence(rec.ref_id).length;
			}
			ref_seq=fai.fetchSequence(rec.ref_name(),start,end);
			res=p.sw_striped(q_seq,ref_seq);
			score_read=res.result.score;
			res.close();
			start=rec.mate_position()-300;
			if(start<0){
				start=0;
			}
			end=rec.mate_position()+300;
			if(end>bam.header.getSequence(rec.mate_ref_id).length){
				end=bam.header.getSequence(rec.mate_ref_id).length;
			}
			ref_seq=fai.fetchSequence(rec.ref_name(),start,end);
			res=p.sw_striped(q_seq,ref_seq);
			score_mate=res.result.score;
			res.close();
			if(score_read>cutoff||score_mate>cutoff){
				if(score_read>score_mate){
					read_class|=0b1000;
					read_class|=0b1_0000;
				}else{
					read_class|=0b1000;
					read_class|=0b10_0000;
				}
			}
		}
		reads[rec.name]=read_class;
	}
	int read_count;
	int clipped;
	int sup;
	int art;
	int aln_r;
	int aln_m;
	int diff_chrom;
	int aln_r_df;
	int aln_m_df;
	foreach(BamRead rec;bam.allReads()){
		read_count++;
		if((rec.name in reads)is null)
			continue;
		uint val=reads[rec.name];
		//writefln("%b",val);
		/+
		ubyte read status encoding
		0	Read is Normal
		1	Read is Softclipped
		2	Read has Supp Alignment
		3   Read is Artifact
		4	Artifact aligns to read region
		5	Artifact aligns to mate region
		6   Mate is on Diff Chrom than read
		7
		+/
		if(val&1){
			out_bam.writeRecord(rec);
		}if((val&0b10)>>1){
			clipped++;
			sc_bam.writeRecord(rec);
		}if((val&0b100)>>2){
			sup++;
		}if((val&0b1000)>>3){
			art++;
			art_bam.writeRecord(rec);
		}if((val&0b10000)>>4){
			aln_r++;
		}if((val&0b100000)>>5){
			aln_m++;
		}if(((val&0b1000000)>>6)&&((val&0b10000)>>4)){
			diff_chrom++;
			aln_r_df++;
		}if(((val&0b1000000)>>6)&&((val&0b100000)>>5)){
			diff_chrom++;
			aln_m_df++;
		}
		//if(((val&0b10000000)>>7)==1){
		//	art_sp++;
		//	art_bam.writeRecord(rec);
		//}

	}
	writeln(read_count);
	writeln(clipped/float(read_count));
	writeln(sup/float(clipped));
	writeln(art/float(clipped));
	writeln(aln_r/float(art));
	writeln(aln_m/float(art));
	writeln(diff_chrom/float(art));
	writeln(aln_r_df/float(diff_chrom));
	writeln(aln_m_df/float(diff_chrom));
	//writeln(art_sep_chr/float(read_count));
}
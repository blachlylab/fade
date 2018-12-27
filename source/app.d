import std.stdio;
import std.algorithm.iteration:splitter;
import std.algorithm:map,reverse;
import std.regex:splitter,regex,match;
import std.range:drop,array;
import std.conv:to;
import std.uni:toUpper;
import std.traits:ReturnType;
import bio.std.hts.bam.reader;
import bio.std.hts.bam.cigar;
import bio.std.hts.bam.writer;
import std.algorithm.iteration:filter;
import std.algorithm:count;
import std.bitmanip;
import dparasail;
import dhtslib.faidx:IndexedFastaFile;

string rc(Range)(Range seq){
	//seq.array.reverse;
	return seq.array.reverse.map!(x=>cast(char)x.complement).array.idup;
}

CigarOperations cigarFromString(string cigar) {
    return match(cigar, regex(`(\d+)([A-Z=])`, "g"))
            .map!(m => CigarOperation(m[1].to!uint, m[2].to!char))
            .array;
}

unittest{
	auto seq="ACGATGATCGATNGT".dup;
	writeln(rc(seq));
}

union ReadStatus{
    ubyte raw;
    //ubyte read status encoding
    mixin(bitfields!(
        //0	Read is Softclipped
        bool,"sc",1,
        //1	Read has Supp Alignment
        bool,"sup",1,
        //2   Read is Artifact
        bool,"art",1,
        //3	Artifact aligns to mate region and not read
        bool,"mate",1,
        //4   Mate is on Diff Chrom than read
        bool,"mate_diff",1,
        //5	5' artifact not 3'
        bool,"five_prime",1,
        //6	same strand
        bool,"same_strand",1,
        //7 supp alignment not close to read or mate
        bool,"far",1,
    ));
}

CigarOperation[2] parse_clips(const CigarOperations cigar){
    CigarOperation[2] clips;
    bool first=true;
    foreach(CigarOperation op;cigar){
        auto is_sc=op.is_query_consuming()&&op.is_clipping();
        if(first&&!is_sc){
            first=false;
        }else if(first&&is_sc){
            clips[0]=op;
        }else if(is_sc){
            clips[1]=op;
        }
    }
    return clips;
}

void align_clip(BamReader * bam,IndexedFastaFile * fai,Parasail * p,BamRead * rec,
        ReadStatus * status, uint clip_len,bool left){
    //SW parasail needed
    string q_seq;
    string ref_seq;
    float cutoff;
    int start,end,score_read,score_mate;
    parasail_query res;
    //if left sofclip ? remove from left : else remove from right
    q_seq=left?rec.sequence[0..clip_len].rc:rec.sequence[$-clip_len..$].rc;
    //set cutoff
    cutoff=q_seq.length*0.75;

    start=rec.position()-300;
    //if start<0: start is zero
    if(start<0){
        start=0;
    }

    end=rec.position()+rec.basesCovered()+300;
    //if end>length of chrom: end is length of chrom
    if(end>bam.header.getSequence(rec.ref_id).length){
        end=bam.header.getSequence(rec.ref_id).length;
    }
    //get read region seq
    ref_seq=fai.fetchSequence(rec.ref_name(),start,end).toUpper;
    //align
    res=p.sw_striped(q_seq,ref_seq);
    //get read region score
    score_read=res.result.score;

    debug{
        writeln(rec.ref_name()," ",start," ",end);
        writeln(rec.name());
        writeln(q_seq);
        writeln(ref_seq);
        writeln(score_read);
    }

    res.close();
    if(rec.is_paired()&&!rec.mate_is_unmapped()){
        //rinse and repeat for mate region
        start=rec.mate_position()-300;
        if(start<0){
            start=0;
        }
        end=rec.mate_position()+300+151; //here we can't know the bases covered so estimate 151 bases
        if(end>bam.header.getSequence(rec.mate_ref_id).length){
            end=bam.header.getSequence(rec.mate_ref_id).length;
        }
        ref_seq=fai.fetchSequence(rec.mate_ref_name(),start,end).toUpper;
        res=p.sw_striped(q_seq,ref_seq);
        score_mate=res.result.score;

        debug{
            writeln(rec.mate_ref_name()," ",start," ",end);
            writeln(rec.name());
            writeln(q_seq);
            writeln(ref_seq);
            writeln(score_mate);
        }

        res.close();
        //choose score from alignments
        if(score_read>cutoff||score_mate>cutoff){
            if(score_read>=score_mate){
                status.art=true;
                status.mate=false;
            }else{
                status.art=true;
                status.mate=true;
            }
        }
    }else{
        if(score_read>cutoff){
            status.art=true;
            status.mate=false;
        }
    }
}
void main(string[] args){
    if(args[1]=="annotate"){
        annotate(args[1..$]);
    }else if(args[1]=="filter"){
        filter(args[1..$],false);
    }else if(args[1]=="clip"){
        filter(args[1..$],true);
    }
}

void annotate(string[] args){
	auto bam = new BamReader(args[1]);
	auto fai=IndexedFastaFile(args[2]);
	auto out_bam=new BamWriter(args[3]);
	out_bam.writeSamHeader(bam.header());
	out_bam.writeReferenceSequenceInfo(bam.reference_sequences());
	//string[2] strands=["+","-"];

	//ubyte[string] reads;
    //ReadStatus[string] reads;

	auto p=Parasail("ACTGN",1,-1,1,3);
	foreach(BamRead rec;bam.allReads){
        ReadStatus status;
		if(!rec.mate_is_unmapped()&&rec.mate_ref_name()!=rec.ref_name())
			//read_class|=0b100_0000;
            status.mate_diff=true;
		if(rec.mate_is_reverse_strand()==rec.is_reverse_strand()){
            status.same_strand=true;
        }
        if(rec.is_supplementary()||
            rec.is_secondary_alignment()||
            rec.is_unmapped()||
            rec.cigar.filter!(x=>x.type=='S').count()==0
        ){
            rec["RS"]=status.raw;
            out_bam.writeRecord(rec);
            continue;
        }
		//read_class|=0b10;
        status.sc=true;
        CigarOperation[2] clips=parse_clips(rec.cigar);
        if(clips[0].length!=0){
            status.five_prime=true;
        }
		if(!rec["SA"].is_nothing()){
			//read_class|=0b100;
            status.sup=true;
			string[] sup=rec["SA"].toString.splitter(",").array;
			if(sup[0]==rec.ref_name()){
				if (
					(sup[1].to!int>rec.position()-300)&&
					(sup[1].to!int<rec.position()+300)//&&
					//(strands[rec.is_reverse_strand]!=sup[2])
				){
                    status.art=true;
					//read_class|=0b1000;
                    status.mate=false;
					//read_class|=0b1_0000;
                    rec["RS"]=status.raw;
                    out_bam.writeRecord(rec);
					continue;
				}
			}else if(sup[0]==rec.mate_ref_name()){
				if (
					(sup[1].to!int>rec.mate_position()-300)&&
					(sup[1].to!int<rec.mate_position()+300)//&&
					//(strands[rec.mate_is_reverse_strand]!=sup[2])
				){
                    status.art=true;
					//read_class|=0b1000;
                    status.mate=true;
					//read_class|=0b10_0000;
                    rec["RS"]=status.raw;
                    out_bam.writeRecord(rec);
					continue;
				}
			}else{
                if(rec.strand()!=sup[2][0]){
                    auto sa_clips=parse_clips(cigarFromString(sup[3]));
                    if((clips[0].length!=0&&sa_clips[1].length<=clips[0].length)||
                        (clips[1].length!=0&&sa_clips[0].length<=clips[1].length)){
                        status.art=true;
                        //read_class|=0b1000;
                        status.mate=false;
                        status.far=true;
                        //read_class|=0b10_0000;
                        rec["RS"]=status.raw;
                        out_bam.writeRecord(rec);
                        continue;
                    }
                }
            }
		}
		//left soft-clip (left on reference not 5' neccesarily)
		if(clips[0].length!=0){
            align_clip(&bam,&fai,&p,&rec,&status,clips[0].length(),true);
		}
		//right soft-clip
		if(clips[1].length()!=0){
            align_clip(&bam,&fai,&p,&rec,&status,clips[1].length(),false);
		}
        rec["RS"]=status.raw;
        assert(rec["RS"].bam_typeid=='C');
		//reads[rec.name]=read_class;
        //if(status.art){
        //    rec["RS"]=ReadSt.rawatus;
        //}
        out_bam.writeRecord(rec);
	}
    out_bam.finish();
}

void clipRead(BamRead * rec,ReadStatus * status){
    auto new_cigar=rec.cigar.dup;
    auto qual=rec.base_qualities.dup;
    if(status.five_prime){
        if(new_cigar[0].type=='H'&&new_cigar[1].type=='S'){
            rec.sequence=rec.sequence[rec.cigar[1].length..$].map!(x=>x.asCharacter).array.idup;
            rec.base_qualities=qual[rec.cigar[1].length..$];
            new_cigar[1]=CigarOperation(new_cigar[0].length+new_cigar[1].length,new_cigar[0].type);
            rec.cigar=new_cigar[1..$];
        }else{
            rec.sequence=rec.sequence[rec.cigar[0].length..$].map!(x=>x.asCharacter).array.idup;
            rec.base_qualities=qual[rec.cigar[0].length..$];
            new_cigar[0]=CigarOperation(new_cigar[0].length,'H');
            rec.cigar=new_cigar;
        }
    }else{
        if(new_cigar[$-1].type=='H'&&new_cigar[$-2].type=='S'){
            rec.sequence=rec.sequence[0..$-rec.cigar[$-2].length].map!(x=>x.asCharacter).array.idup;
            rec.base_qualities=qual[0..$-rec.cigar[$-2].length];
            new_cigar[$-2]=CigarOperation(new_cigar[$-1].length+new_cigar[$-2].length,new_cigar[$-1].type);
            rec.cigar=new_cigar[0..$-1];
        }else{
            rec.sequence=rec.sequence[0..$-rec.cigar[$-1].length].map!(x=>x.asCharacter).array.idup;
            rec.base_qualities=qual[0..$-rec.cigar[$-1].length];
            new_cigar[$-1]=CigarOperation(new_cigar[$-1].length,'H');
            rec.cigar=new_cigar;
        }
    }

}

void filter(string[] args,bool clip){
    auto bam = new BamReader(args[1]);
    auto out_bam=new BamWriter(args[2]);
    out_bam.writeSamHeader(bam.header());
    out_bam.writeReferenceSequenceInfo(bam.reference_sequences());
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
    int read_count;
    int clipped;
    int sup;
    int art;
    int aln_r;
    int aln_m;
    int diff_chrom;
    int aln_r_df;
    int aln_m_df;
    int art_strand;
    int art_far;
    int art_5;
    foreach(BamRead rec;bam.allReads()){
        read_count++;
        ReadStatus val;
        val.raw=cast(ubyte)rec["RS"];
        //if(val.raw==0)
        //    continue;
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
        if(!val.art){
            out_bam.writeRecord(rec);
        }else{
            art++;
            art_bam.writeRecord(rec);
            if(val.mate){
                aln_m++;
            }else if(!val.far){
                aln_r++;
            }else{
                art_far++;
            }
            if(val.mate_diff){
                diff_chrom++;
                if(val.mate){
                    aln_m_df++;
                }else if(!val.far){
                    aln_r_df++;
                }
            }
            if(val.same_strand){
                art_strand++;
            }
            if(val.five_prime){
                art_5++;
            }
            if(clip){
                clipRead(&rec,&val);
                out_bam.writeRecord(rec);
            }
        }
        if(val.sc){
            clipped++;
            sc_bam.writeRecord(rec);
        }if(val.sup){
            sup++;
        }

    }
    out_bam.finish();
    sc_bam.finish();
    art_bam.finish();
    stderr.write("read count:\t");
    stderr.writeln(read_count);
    stderr.write("Clipped %:\t");
    stderr.writeln(clipped/float(read_count));
    stderr.writeln("Of those clipped:");
    stderr.write("\tWith Supplementary alns:\t");
    stderr.writeln(sup/float(clipped));
    stderr.write("\tArtifact:\t");
    stderr.writeln(art/float(clipped));
    stderr.writeln("Of those artifact:");
    stderr.write("\tArtifact on 5' end:\t");
    stderr.writeln(art_5/float(art));
    stderr.write("\tAligned near read:\t");
    stderr.writeln(aln_r/float(art));
    stderr.write("\tAligned near mate:\t");
    stderr.writeln(aln_m/float(art));
    stderr.write("\tAligned near niether:\t");
    stderr.writeln(art_far/float(art));
    stderr.write("\tMate and read on same strand:\t");
    stderr.writeln(art_strand/float(art));
    stderr.write("\tMate and read on different chromosomes:\t");
    stderr.writeln(diff_chrom/float(art));
    stderr.writeln("Of those on different chromosomes:");
    stderr.write("\tArtifact near read:\t");
    stderr.writeln(aln_r_df/float(diff_chrom));
    stderr.write("\tArtifact near mate:\t");
    stderr.writeln(aln_m_df/float(diff_chrom));
    //writeln(art_sep_chr/float(read_count));
}
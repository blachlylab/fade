import std.stdio;
import std.algorithm.iteration:splitter;
import std.algorithm:map,reverse;
import std.regex:splitter,regex,match;
import std.range:drop,array;
import std.conv:to;
import std.uni:toUpper;
import std.traits:ReturnType;
import std.algorithm.iteration:filter;
import std.algorithm:count;
import std.bitmanip;
import std.getopt;
import dparasail;
import dhtslib;
import dhtslib.htslib.sam;
import std.parallelism:defaultPoolThreads;

// string rc(Range)(Range seq){
	//seq.array.reverse;
	// return seq.array.reverse.map!(x=>cast(char)x.complement).array.idup;
// }

const(char)[16] seq_comp_table = [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15];

// extract and reverse-complement soft-clipped portion
char[] extract_soft_clip(SAMRecord * rec, int start, int end){
    ubyte * seq_ptr = (rec.b.data + (rec.b.core.n_cigar<<2) + rec.b.core.l_qname);
    char[] ret;
    ret.length=end-start;
    auto j = end-start-1;
    for(int i = start;i<end;i++){
        ret[j--]=seq_nt16_str[seq_comp_table[((seq_ptr)[(i)>>1] >> ((~(i)&1)<<2) & 0xf)]];
    }
    return ret;
}

unittest{
	auto seq="ACGATGATCGATNGT".dup;
	writeln(rc(seq));
}

/// ubyte bitflag for indicating artifact status
union ReadStatus{
    /// raw ubyte
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

/// report soft clips of a read using a cigar
CigarOp[2] parse_clips(const Cigar cigar){
    CigarOp[2] clips;
    bool first=true;
    foreach(CigarOp op;cigar.ops){
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

//quick and dirty qscore average
ushort avg_qscore(const(char)[] q){
    ushort score=q[0];
    foreach(c;q){
        score+=c;
        score>>=1;
    }
    return score;
}

/// Align the sofclip to the read region or the mate region
void align_clip(SAMReader * bam,IndexedFastaFile * fai,Parasail * p,SAMRecord * rec,
        ReadStatus * status, uint clip_len,bool left){
    string q_seq;
    const(char)[] qual_seq;
    string ref_seq;
    float cutoff;
    int start,end,score_read,score_mate;
    parasail_query res;
    if(clip_len < 6){
        return;
    }
    //if left sofclip ? remove from left : else remove from right
    q_seq=left?extract_soft_clip(rec,0,clip_len).idup:extract_soft_clip(rec,rec.b.core.l_qseq-clip_len,rec.b.core.l_qseq).idup;
    qual_seq=left?rec.qscores!false()[0..clip_len]:rec.qscores!false()[$-clip_len..$];
    // writeln(qual_seq);
    if(avg_qscore(qual_seq)<20) return;
    //set cutoff
    cutoff=q_seq.length*0.75;

    start=rec.pos()-300;
    //if start<0: start is zero
    if(start<0){
        start=0;
    }

    end=rec.pos()+rec.cigar.ref_bases_covered()+300;
    //if end>length of chrom: end is length of chrom
    if(end>bam.target_lens[rec.tid]){
        end=bam.target_lens[rec.tid];
    }
    //get read region seq
    ref_seq=fai.fetchSequence(bam.target_names[rec.tid],start,end).toUpper;
    //align
    res=p.sw_striped(q_seq,ref_seq);
    //get read region score
    score_read=res.result.score;

    debug{
        writeln(rec.tid()," ",start," ",end);
        writeln(rec.queryName());
        writeln(q_seq);
        writeln(ref_seq);
        writeln(score_read);
    }

    res.close();
    if(rec.isPaired()&&!rec.isMateMapped()){
        //rinse and repeat for mate region
        start=rec.matePos()-300;
        if(start<0){
            start=0;
        }
        end=rec.matePos()+300+151; //here we can't know the bases covered so estimate 151 bases
        if(end>bam.target_lens[rec.mateTID]){
            end=bam.target_lens[rec.mateTID];
        }
        ref_seq=fai.fetchSequence(bam.target_names[rec.mateTID],start,end).toUpper;
        res=p.sw_striped(q_seq,ref_seq);
        score_mate=res.result.score;

        debug{
            writeln(rec.mateTID()," ",start," ",end);
            writeln(rec.queryName());
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

int threads;

void main(string[] args){
    auto res=getopt(args,config.bundling,
	"threads|t","threads for parsing the bam file",&threads);
	if (res.helpWanted) {
		defaultGetoptPrinter(
            "usage: ./artifact [annotate] [bam] [reference fasta with fai] [out bam]\n"~
            "usage: ./artifact [filter or clip] [bam] [out bam]\n"~
            "annotate: marks artifact reads in bam tags (must be done first)\n"~
            "filter: removes fragments (read and mate) with artifact (requires queryname sorted bam)\n"~
            "clip: removes artifact region only", res.options);
		stderr.writeln();
		return;
	}
	if(args.length<3){
		writeln("usage: ./artifact [annotate] [bam] [reference fasta with fai] [out bam]");
		return;
	}else{
		if(threads!=0){
			defaultPoolThreads(threads);
		}
    }
    if(args[1]=="annotate"){
        annotate(args[1..$]);
    }else if(args[1]=="filter"){
        filter!false(args[1..$]);
    }else if(args[1]=="clip"){
        filter!true(args[1..$]);
    }
}

void annotate(string[] args){
	auto bam = SAMReader(args[1]);
	auto fai=IndexedFastaFile(args[2]);
	auto out_bam=SAMWriter(args[3],bam.header);
	//string[2] strands=["+","-"];

	//ubyte[string] reads;
    //ReadStatus[string] reads;

	auto p=Parasail("ACTGN",1,-1,1,3);
	foreach(SAMRecord rec;bam.all_records){
        ReadStatus status;
		if(rec.isMateMapped()&&rec.mateTID()!=rec.tid())
			//read_class|=0b100_0000;
            status.mate_diff=true;
		if(rec.mateReversed()==rec.isReversed()){
            status.same_strand=true;
        }
        if(rec.isSupplementary()||
            rec.isSecondary()||
            !rec.isMapped()||
            rec.cigar.ops.filter!(x=>x.op==Ops.SOFT_CLIP).count()==0
        ){
            rec["rs"]=status.raw;
            out_bam.write(&rec);
            continue;
        }
		//read_class|=0b10;
        status.sc=true;
        CigarOp[2] clips=parse_clips(rec.cigar);
        if(clips[0].length!=0){
            status.five_prime=true;
        }
		if(!(rec["SA"].data==null)){
			//read_class|=0b100;
            status.sup=true;
			string[] sup=rec["SA"].toString.splitter(",").array;
			if(sup[0]==bam.target_names[rec.tid]){
				if (
					(sup[1].to!int>rec.pos()-300)&&
					(sup[1].to!int<rec.pos()+300)//&&
					//(strands[rec.is_reverse_strand]!=sup[2])
				){
                    status.art=true;
					//read_class|=0b1000;
                    status.mate=false;
					//read_class|=0b1_0000;
                    rec["rs"]=status.raw;
                    out_bam.write(&rec);
					continue;
				}
			}else if(sup[0]==bam.target_names[rec.mateTID]){
				if (
					(sup[1].to!int>rec.matePos()-300)&&
					(sup[1].to!int<rec.matePos()+300)//&&
					//(strands[rec.mate_is_reverse_strand]!=sup[2])
				){
                    status.art=true;
					//read_class|=0b1000;
                    status.mate=true;
					//read_class|=0b10_0000;
                    rec["rs"]=status.raw;
                    out_bam.write(&rec);
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
                        rec["rs"]=status.raw;
                        out_bam.write(&rec);
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
        rec["rs"]=status.raw;
        assert(rec["rs"].check!ubyte);
		//reads[rec.name]=read_class;
        //if(status.art){
        //    rec["rs"]=ReadSt.rawatus;
        //}
        out_bam.write(&rec);
	}
}

void clipRead(SAMRecord * rec,ReadStatus * status){
    auto new_cigar=rec.cigar.ops.dup;
    auto qual=rec.qscores!false();
    if(status.five_prime){
        if(new_cigar[0].op==Ops.HARD_CLIP && new_cigar[1].op==Ops.SOFT_CLIP){
            rec.sequence=rec.sequence[rec.cigar.ops[1].length..$];
            rec.q_scores!false(qual[rec.cigar.ops[1].length..$]);
            new_cigar[1]=CigarOp(new_cigar[0].length+new_cigar[1].length,Ops.HARD_CLIP);
            rec.cigar=Cigar(new_cigar[1..$]);
        }else{
            rec.sequence=rec.sequence[rec.cigar.ops[0].length..$];
            rec.q_scores!false(qual[rec.cigar.ops[0].length..$]);
            new_cigar[0]=CigarOp(new_cigar[0].length,Ops.HARD_CLIP);
            rec.cigar=Cigar(new_cigar);
        }
    }else{
        if(new_cigar[$-1].op==Ops.HARD_CLIP&&new_cigar[$-2].op==Ops.SOFT_CLIP){
            rec.sequence=rec.sequence[0..$-rec.cigar.ops[$-2].length];
            rec.q_scores!false(qual[0..$-rec.cigar.ops[$-2].length]);
            new_cigar[$-2]=CigarOp(new_cigar[$-1].length+new_cigar[$-2].length,Ops.HARD_CLIP);
            rec.cigar=Cigar(new_cigar[0..$-1]);
        }else{
            rec.sequence=rec.sequence[0..$-rec.cigar.ops[$-1].length];
            rec.q_scores!false(qual[0..$-rec.cigar.ops[$-1].length]);
            new_cigar[$-1]=CigarOp(new_cigar[$-1].length,Ops.HARD_CLIP);
            rec.cigar=Cigar(new_cigar);
        }
    }

}

void filter(bool clip)(string[] args){
    auto bam = SAMReader(args[1]);
    auto out_bam=SAMWriter(args[2],bam.header);
    //auto sc_bam= SAMWriter("sc.bam");
    //sc_bam.writeSamHeader(bam.header());
    //sc_bam.writeReferenceSequenceInfo(bam.reference_sequences());
    //auto db_bam=SAMWriter("db.bam");
    //db_bam.writeSamHeader(bam.header());
    //db_bam.writeReferenceSequenceInfo(bam.reference_sequences());
    auto art_bam=SAMWriter(args[2]~".art.bam",bam.header);
    //auto non_art_bam=SAMWriter("non_art.bam");
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
    static if(clip==true){
        foreach(SAMRecord rec;bam.all_records()){
            read_count++;
            ReadStatus val;
            auto tag=rec["rs"];
            if(tag.data==null){
                out_bam.write(&rec);
                continue;
            }
            val.raw=tag.to!ubyte;
            //if(val.raw==0)
            //    continue;
            //writefln("%b",val);
            /+
            //0,sc	Read is Softclipped
            //1,sup	Read has Supp Alignment        
            //2,art   Read is Artifact
            //3,mate	Artifact aligns to mate region and not read
            //4,mate_diff   Mate is on Diff Chrom than read
            //5,five_prime artifact not 3'
            //6,same_strand	same strand
            //7,far supp alignment not close to read or mate
            +/
            if(!val.art){
                out_bam.write(&rec);
            }else{
                art++;
                art_bam.write(&rec);
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
                clipRead(&rec,&val);
                out_bam.write(&rec);
            }
            if(val.sc){
                clipped++;
                //sc_bam.writeRecord(rec);
            }if(val.sup){
                sup++;
            }
        }
    }
    static if(clip==false){
        import std.algorithm.iteration:chunkBy;
        foreach(recs;bam.all_records.chunkBy!((a,b)=>a.queryName==b.queryName)){
            // recs.map!(x=>x.queryName).each!writeln;
            // bam.fp.fp.bgzf.is_write.writeln;
            auto grouped_reads=recs.array;
            bool art_found=false;
            foreach(rec;grouped_reads){
                read_count++;
                ReadStatus val;
                auto tag=rec["rs"];
                if(tag.data==null){
                    // out_bam.writeRecord(rec);
                    continue;
                }
                val.raw=tag.to!ubyte;
                //if(val.raw==0)
                //    continue;
                //writefln("%b",val);
                /+
                //0,sc	Read is Softclipped
                //1,sup	Read has Supp Alignment        
                //2,art   Read is Artifact
                //3,mate	Artifact aligns to mate region and not read
                //4,mate_diff   Mate is on Diff Chrom than read
                //5,five_prime artifact not 3'
                //6,same_strand	same strand
                //7,far supp alignment not close to read or mate
                +/
                if(val.art){
                    art_found=true;
                    art++;
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
                }
                if(val.sc){
                    clipped++;
                    //sc_bam.writeRecord(rec);
                }if(val.sup){
                    sup++;
                }
            }
            if(art_found){
                foreach(rec;grouped_reads){
                    art_bam.write(&rec);
                }
            }else{
                foreach(rec;grouped_reads){
                    out_bam.write(&rec);
                }
            }
        }
    }
    
    //sc_bam.finish();
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

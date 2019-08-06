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
import std.math:abs;

// string rc(Range)(Range seq){
	//seq.array.reverse;
	// return seq.array.reverse.map!(x=>cast(char)x.complement).array.idup;
// }

const(char)[16] seq_comp_table = [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15];

// extract and reverse-complement soft-clipped portion
pragma(inline,true)
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


/// ushort bitflag for clip status
union ReadStatus{
    ushort raw;
    struct{
        ClipStatus left;
        ClipStatus right;
    }
}
/// ubyte bitflag for indicating artifact status
union ClipStatus {
    /// raw ubyte
    ubyte raw;
    //ubyte read status encoding
    mixin(bitfields!(
        //0	Read is Softclipped
        bool,"sc",1,
        //1	Read has Supp Alignment
        bool,"sup",1,
        //2   Supp is on opposite strand from read
        bool,"sup_opp_strand",1,
        //3   Sc doesn't meet qual cutoff
        bool,"qual",1,
        //4   Read is artifact
        bool,"art",1,
        //5	Artifact aligns to mate region and not read
        bool,"art_mate",1,
        //6   Artifact is greater than 5 bp long but shorter than 15 (TODO: set empirically)
        bool,"art_short",1,
        //7 supp alignment not close to read or mate
        bool,"far",1,
    ));
}

/// report soft clips of a read using a cigar
pragma(inline,true)
CigarOp[2] parse_clips(const Cigar cigar){
    CigarOp[2] clips;
    bool first=true;
    foreach(CigarOp op;cigar.ops){
        //skip hard clips
        if(op.op==Ops.HARD_CLIP) continue;
        auto is_sc=op.op==Ops.SOFT_CLIP;
        // if left soft-clip
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
pragma(inline,true)
ushort avg_qscore(const(char)[] q){
    ushort score=q[0];
    foreach(c;q){
        score+=c;
        score>>=1;
    }
    return score;
}

/// Align the sofclip to the read region or the mate region
string align_clip(SAMReader * bam,IndexedFastaFile * fai,Parasail * p,SAMRecord * rec,
        ReadStatus * status, uint clip_len,bool left){
    string q_seq;
    const(char)[] qual_seq;
    string ref_seq;
    float cutoff;
    int start,end,start_mate,end_mate;
    parasail_query res,res_mate;

    //if clip too short
    if(clip_len <= artifact_floor_length){
        return "";
    }

    //if left sofclip ? remove from left : else remove from right
    q_seq=left?extract_soft_clip(rec,0,clip_len).idup:extract_soft_clip(rec,rec.b.core.l_qseq-clip_len,rec.b.core.l_qseq).idup;
    qual_seq=left?rec.qscores!false()[0..clip_len]:rec.qscores!false()[$-clip_len..$];

    // if quality low
    if(avg_qscore(qual_seq)<qscore_cutoff) {

        if(left) status.left.qual=true;
        else status.right.qual=true;

        return "";
    }

    //set sw score cutoff
    cutoff=q_seq.length*0.75;

    start=rec.pos()-align_buffer_size;

    //if start<0: start is zero
    if(start<0){
        start=0;
    }

    end=rec.pos()+rec.cigar.ref_bases_covered()+align_buffer_size;

    //if end>length of chrom: end is length of chrom
    if(end>bam.target_lens[rec.tid]){
        end=bam.target_lens[rec.tid];
    }

    //get read region seq
    ref_seq=fai.fetchSequence(bam.target_names[rec.tid],start,end).toUpper;
    
    //align
    res=p.sw_striped(q_seq,ref_seq);

    ClipStatus clip = left ? status.left : status.right;

    string align_string;
    if(rec.isPaired() && rec.isMateMapped()){

        //rinse and repeat for mate region
        start_mate=rec.matePos()-align_buffer_size;
        if(start_mate<0){
            start_mate=0;
        }
        //TODO: check for mate cigar to get mate region size 
        end_mate=rec.matePos()+align_buffer_size+mate_size_est; //here we can't know the bases covered so estimate 151 bases
        if(end_mate>bam.target_lens[rec.mateTID]){
            end_mate=bam.target_lens[rec.mateTID];
        }
        ref_seq=fai.fetchSequence(bam.target_names[rec.mateTID],start_mate,end_mate).toUpper;
        res_mate=p.sw_striped(q_seq,ref_seq);

        //choose score from alignments
        if(res.result.score>cutoff||res_mate.result.score>cutoff){
            if(res.result.score>=res_mate.result.score){
                clip.art=true;
                clip.art_mate=false;
                align_string = bam.target_names[rec.tid]~","~
                    (start+res.beg_ref).to!string~","~
                    (start+res.result.end_ref).to!string~","~
                    res.result.score.to!string;
            }else{
                clip.art=true;
                clip.art_mate=true;
                align_string = bam.target_names[rec.mateTID]~","~
                    (start_mate+res_mate.beg_ref).to!string~","~
                    (start_mate+res_mate.result.end_ref).to!string~","~
                    res_mate.result.score.to!string;
            }
            if(clip_len < artifact_short_cutoff){
                clip.art_short=true;
            }
        }
        res_mate.close();
    }else{
        if(res.result.score>cutoff){
            clip.art=true;
            clip.art_mate=false;
            if(clip_len < artifact_short_cutoff){
                clip.art_short=true;
            }
            align_string = bam.target_names[rec.tid]~","~
                (start+res.beg_ref).to!string~","~
                (start+res.result.end_ref).to!string~","~
                res.result.score.to!string;
        }
    }
    if(left) status.left=clip;
    else status.right=clip;
    res.close();
    return align_string;
}

int artifact_floor_length=5;
int artifact_short_cutoff=15;
int align_buffer_size=300;
int mate_size_est=151;
int qscore_cutoff=20;
int sa_size_wiggle=5;
int threads;

void main(string[] args){
    auto res=getopt(args,config.bundling,
	"threads|t","threads for parsing the bam file",&threads);
	if (res.helpWanted) {
		defaultGetoptPrinter(
            "usage: ./fade [annotate] [bam] [reference fasta with fai] [out bam]\n"~
            "usage: ./fade [filter or clip] [bam] [out bam]\n"~
            "annotate: marks artifact reads in bam tags (must be done first)\n"~
            "filter: removes fragments (read and mate) with artifact (requires queryname sorted bam)\n"~
            "clip: removes artifact region only", res.options);
		stderr.writeln();
		return;
	}
	if(args.length<3){
		writeln("usage: ./fade [annotate] [bam] [reference fasta with fai] [out bam]\nusage: ./fade [filter or clip] [bam] [out bam]");
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
        //0	Read is Softclipped
        // sc
        //1	Read has Supp Alignment
        // sup
        //2   Supp is on opposite strand from read
        // sup_opp_strand
        //3   Sc doesn't meet qual cutoff
        // qual
        //4   Read is artifact
        // art
        //5	Artifact aligns to mate region and not read
        // art_mate
        //6   Artifact is greater than 5 bp long but shorter than 15 (TODO: set empirically)
        // art_short
        //7 supp alignment not close to read or mate
        // far
	auto p=Parasail("ACTGN",1,-1,1,3);
	foreach(SAMRecord rec;bam.all_records){
        ReadStatus status;
        if(rec.isSupplementary()||
            rec.isSecondary()||
            !rec.isMapped()||
            rec.cigar.ops.filter!(x=>x.op==Ops.SOFT_CLIP).count()==0
        ){
            rec["rs"]=status.raw;
            out_bam.write(&rec);
            continue;
        }
        CigarOp[2] clips=parse_clips(rec.cigar);
        if(clips[0].length!=0) status.left.sc=true;
        if(clips[1].length!=0) status.right.sc=true;
		if(rec["SA"].data!=null){
			string[] sup=rec["SA"].toString.splitter(",").array;
            if(sup[2][0]!=rec.strand){
                status.left.sup_opp_strand=true;
                status.right.sup_opp_strand=true;
                auto sa_cigar =cigarFromString(sup[3]);
                if(clips[0].length!=0){
                    if(sa_cigar.ops[0].op!=Ops.SOFT_CLIP){
                        if(abs(sa_cigar.ops[0].length-clips[0].length)<=sa_size_wiggle){
                            status.left.sup=true;
                            status.left.art=true;
                            rec["rs"]=status.raw;
                            rec["am"]=sup[0]~","~sup[1]~","~(sup[1].to!int+sa_cigar.ops[0].length).to!string~","~sa_cigar.ops[0].length.to!string~";";
                            out_bam.write(&rec);
                            continue;
                        }
                    }
                }
                if(clips[1].length!=0){
                    if(sa_cigar.ops[$-1].op!=Ops.SOFT_CLIP){
                        if(abs(sa_cigar.ops[$-1].length-clips[1].length)<=sa_size_wiggle){
                            status.right.sup=true;
                            status.right.art=true;
                            rec["rs"]=status.raw;
                            rec["am"]=";"~sup[0]~","~sup[1]~","~(sup[1].to!int+sa_cigar.ops[0].length).to!string~","~sa_cigar.ops[0].length.to!string;
                            out_bam.write(&rec);
                            continue;
                        }
                    }
                }
            }
		}
		//left soft-clip (left on reference not 5' neccesarily)
        string align_str;
		if(clips[0].length!=0){
            align_str~=align_clip(&bam,&fai,&p,&rec,&status,clips[0].length(),true);
		}
        align_str~=";";
		//right soft-clip
		if(clips[1].length()!=0){
            align_str~=align_clip(&bam,&fai,&p,&rec,&status,clips[1].length(),false);
		}
        // writeln(status.raw);
        rec["rs"]=status.raw;
        rec["am"]=align_str;
        assert(rec["rs"].check!ubyte);
        assert(rec["am"].check!string);
        out_bam.write(&rec);
	}
}

void clipRead(SAMRecord * rec,ReadStatus * status){
    auto new_cigar=rec.cigar.ops.dup;
    auto qual=rec.qscores!false();
    if(status.left.art){
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
    }else if(status.right.art){
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

struct Stats {
    int read_count;
    int clipped;
    int sup;
    int sup_opp;
    int art;
    int qual;
    int art_mate;
    int art_short;

    int aln_l;
    int aln_r;
    //0	Read is Softclipped
    // sc
    //1	Read has Supp Alignment
    // sup
    //2   Supp is on opposite strand from read
    // sup_opp_strand
    //3   Sc doesn't meet qual cutoff
    // qual
    //4   Read is artifact
    // art
    //5	Artifact aligns to mate region and not read
    // art_mate
    //6   Artifact is greater than 5 bp long but shorter than 15 (TODO: set empirically)
    // art_short
    //7 supp alignment not close to read or mate
    // far
    void parse(ReadStatus rs){
        read_count++;
        clipped+=(rs.left.raw|rs.right.raw)&1;
        sup+=(rs.left.raw>>1|rs.right.raw>>1)&1;
        sup_opp+=(rs.left.raw>>2|rs.right.raw>>2)&1;
        qual+=(rs.left.raw>>3|rs.right.raw>>3)&1;
        art+=(rs.left.raw>>4|rs.right.raw>>4)&1;
        art_mate+=(rs.left.raw>>5|rs.right.raw>>5)&1;
        art_short+=(rs.left.raw>>6|rs.right.raw>>6)&1;

        aln_l+=(rs.left.raw>>4)&1;
        aln_r+=(rs.right.raw>>4)&1;
    }
    void print(){
        stderr.write("read count:\t");
        stderr.writeln(read_count);
        stderr.write("Clipped %:\t");
        stderr.writeln(clipped/float(read_count));
        stderr.write("With Supplementary alns:\t");
        stderr.writeln(sup/float(read_count));
        stderr.write("With Supplementary alns on opposite strand:\t");
        stderr.writeln(sup_opp/float(read_count));
        stderr.write("Artifact rate:\t");
        stderr.writeln(art/float(read_count));
        stderr.write("Artifact rate left only:\t");
        stderr.writeln(aln_l/float(read_count));
        stderr.write("Artifact rate right only:\t");
        stderr.writeln(aln_r/float(read_count));
        stderr.write("Artifact rate short (<15bp):\t");
        stderr.writeln(art_short/float(read_count));
        stderr.write("Artifact rate mate:\t");
        stderr.writeln(art_mate/float(read_count));
    }
}

void filter(bool clip)(string[] args){
    auto bam = SAMReader(args[1]);
    auto out_bam=SAMWriter(args[2],bam.header);
    auto art_bam=SAMWriter(args[2]~".art.bam",bam.header);
    Stats stats;
    static if(clip==true){
        foreach(SAMRecord rec;bam.all_records()){
            ReadStatus val;
            auto tag=rec["rs"];
            if(tag.data==null){
                out_bam.write(&rec);
                continue;
            }
            val.raw=tag.to!ushort;
            stats.parse(val);
            if(!(val.left.art | val.right.art)){
                out_bam.write(&rec);
            }else{
                art_bam.write(&rec);
                clipRead(&rec,&val);
                out_bam.write(&rec);
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
                ReadStatus val;
                auto tag=rec["rs"];
                if(tag.data==null){
                    // out_bam.writeRecord(rec);
                    continue;
                }
                val.raw=tag.to!ushort;
                stats.parse(val);
                if(val.left.art|val.right.art){
                    art_found=true;
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
        stats.print;
    }
}

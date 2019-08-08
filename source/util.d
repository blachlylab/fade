module util;
import dhtslib;
import dhtslib.htslib.hts:seq_nt16_str;

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
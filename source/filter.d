module filter;
import std.array:array;
import dhtslib;
import readstatus;
import stats;


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

void filter(bool clip)(string[] args){
    auto bam = SAMReader(args[1]);
    auto out_bam=SAMWriter(args[2],bam.header);
    auto art_bam=SAMWriter(args[2]~".art.bam",bam.header);
    Stats stats;
    static if(clip==true){
        foreach(SAMRecord rec;bam.all_records()){
            stats.read_count++;
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
                stats.read_count++;
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

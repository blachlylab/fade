module filter;
import std.array:array;
import std.algorithm.iteration:splitter;
import std.range:drop,retro;
import std.conv:to;
import std.stdio;
import dhtslib;
import dhtslib.htslib.sam;
import readstatus;
import stats;
import util;


void clipRead(SAMRecord * rec,ReadStatus * status){
    auto new_cigar=rec.cigar.ops.dup;
    auto qual=rec.qscores!false();
    if(status.art_left){
        //get artifact cigar
        writeln((*rec)["am"].toString.splitter(";").front.splitter(",").drop(2).front);
        auto art_cigar = cigarFromString((*rec)["am"].toString.splitter(";").front.splitter(",").drop(2).front);

        //assert left side is soft-clipped 
        if(art_cigar.ops[0].op!=Ops.SOFT_CLIP) return;
        if(art_cigar.ref_bases_covered>rec.length) return;

        //trim sequence
        rec.sequence=rec.sequence[art_cigar.ref_bases_covered..$];
        rec.q_scores!false(qual[art_cigar.ref_bases_covered..$]);

        if(new_cigar[0].op==Ops.HARD_CLIP) new_cigar=new_cigar[1..$];
        assert(new_cigar[0].op==Ops.SOFT_CLIP);
        auto overlap_len =  art_cigar.ref_bases_covered-new_cigar[0].length;
        //remove soft clip
        new_cigar=new_cigar[1..$];
        while(overlap_len>0){
            if(new_cigar[0].is_query_consuming){
                if(new_cigar[0].length >= overlap_len){
                    new_cigar[0].length = new_cigar[0].length - overlap_len;
                    overlap_len = 0;
                }else{
                    overlap_len -= new_cigar[0].length;
                    new_cigar=new_cigar[1..$];
                }
            }else{
                new_cigar=new_cigar[1..$];
            }
            if(new_cigar.length==0) break;
        } 
    }
    if(status.art_right){
        //get artifact cigar
        auto art_cigar = cigarFromString((*rec)["am"].toString.splitter(";").drop(1).front.splitter(",").drop(2).front);

        //assert right side is soft-clipped 
        if(art_cigar.ops[$-1].op!=Ops.SOFT_CLIP) return;
        if(art_cigar.ref_bases_covered>rec.length) return;

        //trim sequence
        rec.sequence=rec.sequence[0..$-art_cigar.ref_bases_covered];
        rec.q_scores!false(qual[0..$-art_cigar.ref_bases_covered]);

        if(new_cigar[$-1].op==Ops.HARD_CLIP) new_cigar=new_cigar[0..$-1];
        assert(new_cigar[$-1].op==Ops.SOFT_CLIP);
        auto overlap_len =  art_cigar.ref_bases_covered-new_cigar[$-1].length;
        //remove soft clip
        new_cigar=new_cigar[0..$-1];
        while(overlap_len>0){
            if(new_cigar[$-1].is_query_consuming){
                if(new_cigar[$-1].length >=overlap_len){
                    new_cigar[$-1].length = new_cigar[$-1].length - overlap_len;
                    overlap_len = 0;
                }else{
                    overlap_len -= new_cigar[$-1].length;
                    new_cigar=new_cigar[0..$-1];
                }
            }else{
                new_cigar=new_cigar[0..$-1];
            }
            if(new_cigar.length==0) break;
        }
    }
    rec.cigar=Cigar(new_cigar);
}

SAMRecord makeArtifactRecord(SAMRecord * original,bool left, bool mate){
    auto rec =  new SAMRecord(bam_dup1(original.b));
    rec.sequence = reverse_complement_sam_record(&rec);
    writeln(original.queryName);
    writeln(rec["am"].toString);
    // rec.q_scores!false(cast(char[])(cast(ubyte[])((*original).qscores!false).retro.array));
    if(left){
        rec.cigar=cigarFromString(rec["am"].toString.splitter(";").front.splitter(",").drop(2).front);
    }else{
        rec.cigar=cigarFromString(rec["am"].toString.splitter(";").drop(1).front.splitter(",").drop(2).front);
    }
    //set unpaired change strand and supplementary
    rec.b.core.flag&=0b1111_1111_0011_1100;
    rec.b.core.flag^=0b0000_0000_0001_0000;
    rec.b.core.flag|=0b0000_1001_0000_0000;
    if(mate){
        rec.tid=rec.mateTID;
    }
    if(left){
        rec.pos=rec["am"].to!string.splitter(";").front.splitter(",").drop(1).front.to!int;
    }else{
        rec.pos=rec["am"].to!string.splitter(";").drop(1).front.splitter(",").drop(1).front.to!int;
    }
    return rec;
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
            val.raw=tag.to!ubyte;
            stats.parse(val);
            if(!(val.art_left | val.art_right)){
                out_bam.write(&rec);
            }else{
                art_bam.write(&rec);
                if(val.art_left){
                    auto art_rec =makeArtifactRecord(&rec,true,val.mate_left);
                    art_bam.write(&art_rec);
                }
                if(val.art_right){
                    auto art_rec = makeArtifactRecord(&rec,false,val.mate_right);
                    art_bam.write(&art_rec);
                }
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
                val.raw=tag.to!ubyte;
                stats.parse(val);
                if(val.art_left|val.art_right){
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
    }
    stats.print;
}

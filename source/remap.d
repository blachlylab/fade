module remap;
import dhtslib.sam;
import dhtslib.cigar;
import dhtslib.htslib.sam: BAM_FREVERSE;
import std.algorithm.mutation: reverse;
import std.conv: to;
import readstatus;
import util;

void remapArtifacts(string[] args, ubyte con){
    import std.array:join,split;
    import std.format:format;
    auto bam = SAMReader(args[1]);
    auto out_bam=getWriter(con,bam.header);
    foreach(SAMRecord rec;bam.all_records()){
        auto tag=rec["rs"];
        if(!tag.exists) continue;
        ReadStatus rs;
        rs.raw=tag.to!ubyte;
        if(!(rs.art_left|rs.art_right)) continue;
        tag=rec["am"];
        if(!tag.exists) continue;
        auto am = tag.toString;
        auto am_split = am.split(';');
        if(rs.art_left){
            SAMRecord newRec = new SAMRecord();
            auto am_fields = am_split[0].split(',');
            
            newRec.queryName = rec.queryName.idup;
            newRec.tid = bam.target_id(am_fields[0]);
            newRec.pos = am_fields[1].to!int;
            if(rec.isReversed()){
                newRec.flag = 0;
            }else{
                newRec.flag = 0 | BAM_FREVERSE;
            }
            newRec.sequence = reverse_complement_sam_record(&rec);
            newRec.q_scores!false(rec.qscores!false.dup.reverse);
            newRec.cigar = cigarFromString(am_fields[2]);
            out_bam.write(&newRec);
        }
        if(rs.art_right){
            SAMRecord newRec = new SAMRecord();
            auto am_fields = am_split[1].split(',');
            
            newRec.queryName = rec.queryName.idup;
            newRec.tid = bam.target_id(am_fields[0]);
            newRec.pos = am_fields[1].to!int;
            if(rec.isReversed()){
                newRec.flag = 0;
            }else{
                newRec.flag = 0 | BAM_FREVERSE;
            }
            newRec.sequence = reverse_complement_sam_record(&rec);
            newRec.q_scores!false(rec.qscores!false.dup.reverse);
            newRec.cigar = cigarFromString(am_fields[2]);
            out_bam.write(&newRec);
        }
        out_bam.close;
    }
}
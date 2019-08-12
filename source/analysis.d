module analysis;
import std.stdio;
import std.uni:toUpper;
import std.conv:to;
import dhtslib;
import dparasail;
import readstatus;
import util;
import main:artifact_floor_length,artifact_short_cutoff,qscore_cutoff,align_buffer_size,mate_size_est;

/// Align the sofclip to the read region or the mate region
string align_clip(bool left)(SAMReader * bam,string fai_f,Parasail * p,SAMRecord * rec,
        ReadStatus * status, uint clip_len){
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
    q_seq=reverse_complement_sam_record(rec).idup;

    //set sw score cutoff
    cutoff= clip_len*0.9;

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
    auto fai = IndexedFastaFile(fai_f);
    //get read region seq
    ref_seq=fai.fetchSequence(bam.target_names[rec.tid],start,end).toUpper;
    
    //align
    res=p.sw_striped(q_seq,ref_seq);

    // ClipStatus clip = left ? status.left : status.right;
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
        static if(left){
            if(res.cigar.ops[$-1].op==Ops.EQUAL||res_mate.cigar.ops[$-1].op==Ops.EQUAL){
                //choose score from alignments
                if(res.result.score>cutoff||res_mate.result.score>cutoff){
                    if(res.result.score>=res_mate.result.score){
                        status.art_left=true;
                        status.mate_left=false;
                        align_string = bam.target_names[rec.tid]~","~
                            (start+res.beg_ref).to!string~","~
                            (start+res.result.end_ref).to!string~","~
                            res.result.score.to!string~","~
                            res.cigar.toString;
                    }else{
                        status.art_left=true;
                        status.mate_left=true;
                        align_string = bam.target_names[rec.mateTID]~","~
                            (start_mate+res_mate.beg_ref).to!string~","~
                            (start_mate+res_mate.result.end_ref).to!string~","~
                            res_mate.result.score.to!string~","~
                            res.cigar.toString;
                    }
                }
            }
        }else{
            if(res.cigar.ops[0].op==Ops.EQUAL||res_mate.cigar.ops[0].op==Ops.EQUAL){
                //choose score from alignments
                if(res.result.score>cutoff||res_mate.result.score>cutoff){
                    if(res.result.score>=res_mate.result.score){
                        status.art_right=true;
                        status.mate_right=false;
                        align_string = bam.target_names[rec.tid]~","~
                            (start+res.beg_ref).to!string~","~
                            (start+res.result.end_ref).to!string~","~
                            clip_len.to!string~","~
                            res.result.score.to!string~","~
                            res.cigar.toString;
                    }else{
                        status.art_right=true;
                        status.mate_right=true;
                        align_string = bam.target_names[rec.mateTID]~","~
                            (start_mate+res_mate.beg_ref).to!string~","~
                            (start_mate+res_mate.result.end_ref).to!string~","~
                            clip_len.to!string~","~
                            res_mate.result.score.to!string~","~
                            res.cigar.toString;
                    }
                }
            }
        }
        res_mate.close();
    }else{
        static if(left){
            if(res.cigar.ops[$-1].op==Ops.EQUAL){
                if(res.result.score>cutoff){
                    status.art_left=true;
                    status.mate_left=false;
                    align_string = bam.target_names[rec.tid]~","~
                        (start+res.beg_ref).to!string~","~
                        (start+res.result.end_ref).to!string~","~
                        res.result.score.to!string;
                }
            }
        }else{
            if(res.cigar.ops[0].op==Ops.EQUAL){
                if(res.result.score>cutoff){
                    status.art_right=true;
                    status.mate_right=false;
                    align_string = bam.target_names[rec.tid]~","~
                        (start+res.beg_ref).to!string~","~
                        (start+res.result.end_ref).to!string~","~
                        res.result.score.to!string;
                }
            }
        }
    }
    res.close();
    return align_string;
}


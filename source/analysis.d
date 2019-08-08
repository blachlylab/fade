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
                if(clip_len < artifact_short_cutoff){
                    clip.art_short=true;
                }
            }else{
                clip.art=true;
                clip.art_mate=true;
                align_string = bam.target_names[rec.mateTID]~","~
                    (start_mate+res_mate.beg_ref).to!string~","~
                    (start_mate+res_mate.result.end_ref).to!string~","~
                    res_mate.result.score.to!string;
                if(clip_len < artifact_short_cutoff){
                    clip.art_short=true;
                }
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


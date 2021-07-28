module analysis;
import std.stdio;
import std.uni : toUpper;
import std.conv : to;
import core.sync.mutex : Mutex;
import dhtslib;
import dparasail;
import readstatus;
import util;
import std.stdio;

struct Align_Result
{
    string alignment;
    const(char)[] bq;
    string stem_loop;
    string stem_loop_rc;

}

/// Align the sofclip to the read region or the mate region
Align_Result align_clip(bool left)(SAMReader* bam, IndexedFastaFile* fai, Parasail* p,
        SAMRecord rec, ReadStatus* status, uint clip_len, Mutex* m,
        int artifact_floor_length, int align_buffer_size)
{
    string q_seq;
    string ref_seq;
    float cutoff;
    long start, end;
    // parasail_query res;
    Align_Result alignment;

    //if clip too short
    if (clip_len <= artifact_floor_length)
    {
        return alignment;
    }

    //if left sofclip ? remove from left : else remove from right
    q_seq = reverse_complement_sam_record(rec).idup;

    //set sw score cutoff
    cutoff = clip_len * 0.9 * 2;

    start = rec.pos() - align_buffer_size;

    //if start<0: start is zero
    if (start < 0)
    {
        start = 0;
    }

    end = rec.pos() + rec.cigar.alignedLength() + align_buffer_size;

    //if end>length of chrom: end is length of chrom
    if (end > bam.header.targetLength(rec.tid))
    {
        end = bam.header.targetLength(rec.tid);
    }

    m.lock();
    //get read region seq
    auto coords = ChromCoordinates!(CoordSystem.zbho)(bam.header.targetName(rec.tid).idup,ZBHO(start, end));
    ref_seq = fai.fetchSequence(coords).toUpper;
    m.unlock();

    //align
    auto res = p.sw_striped(q_seq, ref_seq);
    // ClipStatus clip = left ? status.left : status.right;
    if ((res.cigar.length == 0) | (res.cigar.length > 10))
        return alignment;

    static if (left)
    {
        if (res.cigar[$ - 1].op == Ops.EQUAL)
        {
            if (res.score > cutoff)
            {
                auto clips = parse_clips(res.cigar);
                if (clips[1].length != 0 || clips[0].length == 0)
                    return alignment;

                status.art_left = true;
                status.mate_left = false;
                alignment.alignment = bam.header.targetName(rec.tid).idup ~ "," ~ (start + res.position)
                    .to!string ~ "," ~ res.cigar.toString;
                auto overlap = start + res.position >= rec.pos - clip_len
                    ? start + res.position - (rec.pos - clip_len) : 0;
                auto plen = (rec.length - clips[0].length) + (overlap);
                plen = plen > rec.length ? rec.length : plen;
                alignment.stem_loop = rec.sequence[0 .. plen].idup;
                alignment.stem_loop_rc = q_seq[$ - plen .. $];
                alignment.bq = rec.qscoresPhredScaled[0 .. plen];
            }
        }
    }
    else
    {
        if (res.cigar[0].op == Ops.EQUAL)
        {
            if (res.score > cutoff)
            {
                auto clips = parse_clips(res.cigar);
                if (clips[0].length != 0 || clips[1].length == 0)
                    return alignment;

                status.art_right = true;
                status.mate_right = false;
                alignment.alignment = bam.header.targetName(rec.tid).idup ~ "," ~ (start + res.position)
                    .to!string ~ "," ~ res.cigar.toString;
                auto overlap = rec.pos + rec.cigar.alignedLength + clip_len >= start
                    + res.position + res.cigar.alignedLength
                    ? (rec.pos + rec.cigar.alignedLength + clip_len) - (
                            start + res.position + res.cigar.alignedLength) : 0;
                auto plen = (rec.length - clips[1].length) + (overlap);
                plen = plen > rec.length ? rec.length : plen;
                alignment.stem_loop = rec.sequence[$ - plen .. $].idup;
                alignment.stem_loop_rc = q_seq[0 .. plen];
                alignment.bq = rec.qscoresPhredScaled[$ - plen .. $];
            }
        }
    }

    return alignment;
}

/// Align the sofclip to the read region or the mate region
string self_align(bool left)(SAMReader* bam, string fai_f, Parasail* p,
        SAMRecord* rec, ReadStatus* status, uint clip_len)
{
    string q_seq;
    const(char)[] qual_seq;
    string ref_seq;
    float cutoff;
    int start, end, start_mate, end_mate;
    parasail_query res, res_mate;

    //if clip too short
    if (clip_len <= artifact_floor_length)
    {
        return "";
    }

    //if left sofclip ? remove from left : else remove from right
    q_seq = reverse_complement_sam_record(rec).idup;

    //set sw score cutoff
    cutoff = clip_len * 0.9;

    //align
    res = p.sw_striped(q_seq, rec.sequence.idup);

    // ClipStatus clip = left ? status.left : status.right;
    string align_string;
    if ((res.cigar.ops.length == 0) | (res.cigar.ops.length > 10))
        return "";
    static if (left)
    {
        align_string = rec.queryName.idup ~ "," ~ (res.position).to!string ~ "," ~ res
            .cigar.toString;
    }
    else
    {
        align_string = rec.queryName.idup ~ "," ~ (res.position).to!string ~ "," ~ res
            .cigar.toString;
    }
    res.close();
    return align_string;
}

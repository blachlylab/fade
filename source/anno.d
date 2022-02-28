module anno;
import std.algorithm : filter, map, each;
import std.conv : to;
import std.algorithm.searching : count;
import std.algorithm.iteration : splitter;
import std.array : array;
import core.sync.mutex : Mutex;
import std.parallelism : parallel;
import dhtslib;
import htslib.hts_log;
import dparasail;
import readstatus;
import analysis;
import util;

int annotate(string cl,string[] args, ubyte con, int artifact_floor_length, int align_buffer_size)
{
    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    hts_log_warning("fade annotate","Output SAM/BAM will not be sorted (regardless of prior sorting)");
    // open bam read and writer
    // also modify header
    auto bam = SAMReader(args[1]);
    auto fai = IndexedFastaFile(args[2]);
    auto header = bam.header.dup;
    header.addLine(
        RecordType.PG, 
        "ID", "fade-annotate",
        "PN", "fade",
        "VN", VERSION,
        "PP", header.valueByPos(RecordType.PG, header.numRecords(RecordType.PG) - 1, "ID"), 
        "CL", cl
        );
    auto out_bam = getWriter(con, header);
    
    // initialize parasail profile
    auto p = Parasail("ACTGN", 10, 2, 2, -3);

    // need a mutex for non-thread-safe bam writing
    auto mwrite = new Mutex();
    // need a mutex for non-thread-safe fai reading
    auto mfai = new Mutex();

    // process reads in parallel
    foreach(rec; parallel(bam.allRecords))
    {
        auto newRec = annotateTask(rec, &p, fai, mfai, artifact_floor_length, align_buffer_size);
        mwrite.lock;
        out_bam.write(newRec);
        mwrite.unlock;
    }
    return 0;
}

/// Performs enz-frag artifact detection/annotation per SAMRecord
static SAMRecord annotateTask(SAMRecord rec, Parasail * p, IndexedFastaFile fai, Mutex mfai, int artifact_floor_length, int align_buffer_size)
{
    ReadStatus status;

    // if read is supp, sec, or not mapped
    // set status and write out
    if (rec.isSupplementary() || rec.isSecondary() || !rec.isMapped()
            || rec.cigar[].filter!(x => x.op == Ops.SOFT_CLIP).count() == 0)
    {
        rec["rs"] = status.raw;
        return rec;
    }

    // check clips
    CigarOp[2] clips = parse_clips(rec.cigar);
    if (clips[0].length != 0 || clips[1].length != 0)
        status.sc = true;

    // check sup alignment
    if (rec["SA"].exists)
        status.sup = true;

    // if left soft-clip (left on reference not 5' neccesarily)
    // perform realignment of rc'd read
    Align_Result align_1, align_2;
    if (clips[0].length != 0)
    {
        align_1 = align_clip!true(fai, mfai, p, rec, &status,
                clips[0].length(), artifact_floor_length, align_buffer_size);
    }

    // if right soft-clip (right on reference not 3' neccesarily)
    // perform realignment of rc'd read
    if (clips[1].length() != 0)
    {
        align_2 = align_clip!false(fai, mfai, p, rec, &status,
                clips[1].length(), artifact_floor_length, align_buffer_size);
    }
    
    // report read status
    rec["rs"] = status.raw;

    // if artifact indicated
    // write other data
    if(status.art_left | status.art_right){
        // report both artifact cigars
        rec["am"] = align_1.alignment ~ ";" ~ align_2.alignment;
        // report both predicted stem loops
        rec["as"] = align_1.stem_loop ~ ";" ~ align_2.stem_loop;
        // report both rc'd predicted stem loops
        rec["ar"] = align_1.stem_loop_rc ~ ";" ~ align_2.stem_loop_rc;
        // report both artifact quality seqs
        rec["ab"] = (align_1.bq ~ ";" ~ align_2.bq).idup;
    }

    return rec;
}
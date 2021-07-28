module anno;
import std.algorithm : filter;
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

void annotate(string cl,string[] args, ubyte con, int artifact_floor_length, int align_buffer_size)
{
    hts_log_warning("fade annotate","Output SAM/BAM will not be sorted (reguardless of prior sorting)");
    // open bam read and writer
    // also modify header
    auto bam = SAMReader(args[1]);
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
    auto p = Parasail("ACTGN", 2, -3, 10, 2);
    auto m = new Mutex();
    auto fai = IndexedFastaFile(args[2]);

    // process reads in parallel
    foreach (SAMRecord rec; parallel(bam.allRecords))
    {
        ReadStatus status;

        // if read is supp, sec, or not mapped
        // set status and write out
        if (rec.isSupplementary() || rec.isSecondary() || !rec.isMapped()
                || rec.cigar[].filter!(x => x.op == Ops.SOFT_CLIP).count() == 0)
        {
            rec["rs"] = status.raw;
            m.lock;
            out_bam.write(rec);
            m.unlock;
            continue;
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
            align_1 = align_clip!true(&bam, &fai, &p, rec, &status,
                    clips[0].length(), &m, artifact_floor_length, align_buffer_size);
        }

        // if right soft-clip (right on reference not 3' neccesarily)
        // perform realignment of rc'd read
        if (clips[1].length() != 0)
        {
            align_2 = align_clip!false(&bam, &fai, &p, rec, &status,
                    clips[1].length(), &m, artifact_floor_length, align_buffer_size);
        }
        
        // report read status
        rec["rs"] = status.raw;

        // report both artifact cigars
        rec["am"] = align_1.alignment ~ ";" ~ align_2.alignment;
        // report both predicted stem loops
        rec["as"] = align_1.stem_loop ~ ";" ~ align_2.stem_loop;
        // report both rc'd predicted stem loops
        rec["ar"] = align_1.stem_loop_rc ~ ";" ~ align_2.stem_loop_rc;
        // report both artifact quality seqs
        rec["ab"] = (cast(char[])(align_1.bq ~ ";" ~ align_2.bq)).idup;

        // write record
        m.lock;
        out_bam.write(rec);
        m.unlock;
    }
    out_bam.close;
}

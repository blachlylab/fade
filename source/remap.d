module remap;
import dhtslib.sam;
import dhtslib.coordinates;
import htslib.sam : BAM_FREVERSE;
import htslib.hts_log;
import std.algorithm.mutation : reverse;
import std.conv : to;
import readstatus;
import util;

void remapArtifacts(string cl, string[] args, ubyte con)
{
    hts_log_warning("fade extract","Output SAM/BAM will not be sorted");
    import std.array : join, split;
    import std.format : format;

    auto bam = SAMReader(args[1]);
    auto header = bam.header.dup;
    header.addLine(
        RecordType.PG, 
        "ID", "fade-extract",
        "PN", "fade",
        "VN", VERSION,
        "PP", header.valueByPos(RecordType.PG, header.numRecords(RecordType.PG) - 1, "ID"), 
        "CL", cl
        );
    auto out_bam = getWriter(con, header);
    
    foreach (SAMRecord rec; bam.allRecords())
    {
        auto tag = rec["rs"];
        if (!tag.exists)
            continue;
        ReadStatus rs;
        rs.raw = tag.to!ubyte;
        if (!(rs.art_left | rs.art_right))
            continue;
        tag = rec["am"];
        if (!tag.exists)
            continue;
        auto am = tag.toString;
        auto am_split = am.split(';');
        if (rs.art_left)
        {
            SAMRecord newRec = SAMRecord(out_bam.header);
            auto am_fields = am_split[0].split(',');

            newRec.queryName = rec.queryName.idup;
            newRec.tid = bam.header.targetId(am_fields[0]);
            newRec.pos = ZB(am_fields[1].to!long);
            if (rec.isReversed())
            {
                newRec.flag = 0;
            }
            else
            {
                newRec.flag = 0 | BAM_FREVERSE;
            }
            newRec.sequence = reverse_complement_sam_record(rec);
            newRec.qscores(rec.qscores.dup.reverse);
            newRec.cigar = cigarFromString(am_fields[2]);
            out_bam.write(newRec);
        }
        if (rs.art_right)
        {
            SAMRecord newRec = SAMRecord(out_bam.header);
            auto am_fields = am_split[1].split(',');

            newRec.queryName = rec.queryName.idup;
            newRec.tid = bam.header.targetId(am_fields[0]);
            newRec.pos = ZB(am_fields[1].to!long);
            if (rec.isReversed())
            {
                newRec.flag = 0;
            }
            else
            {
                newRec.flag = 0 | BAM_FREVERSE;
            }
            newRec.sequence = reverse_complement_sam_record(rec);
            newRec.qscores(rec.qscores.dup.reverse);
            newRec.cigar = cigarFromString(am_fields[2]);
            out_bam.write(newRec);
        }
    }
}

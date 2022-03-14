module filter;
import std.array : array;
import std.algorithm: splitter, isSorted, sort, map;
import std.range : drop, retro, take;
import std.conv : to, parse, ConvException;
import std.typecons : Flag, Yes, No;
import std.stdio;
import dhtslib;
import htslib.sam;
import htslib.hts_log;
import readstatus;
import stats;
import util;

void clipRead(SAMRecord* rec, ReadStatus* status)
{
    auto new_cigar = rec.cigar.dup;
    auto qual = rec.qscores();
    if (status.art_left)
    {
        //get artifact cigar
        // writeln((*rec)["am"].toString.splitter(";").front.splitter(",").drop(2).front);
        auto art_cigar = cigarFromString((*rec)["am"].toString.splitter(";")
                .front.splitter(",").drop(2).front);

        //assert left side is soft-clipped 
        if (art_cigar[0].op != Ops.SOFT_CLIP)
        {
            writeln((*rec)["am"]);
            debug assert(false);
        else return;
        }
        if (art_cigar.alignedLength > rec.length)
        {
            writeln((*rec)["am"]);
            debug assert(false);
        else return;
        }

        if (new_cigar[0].op == Ops.HARD_CLIP)
            new_cigar = new_cigar[1 .. $];
        assert(new_cigar[0].op == Ops.SOFT_CLIP);

        rec.b.core.pos += art_cigar.alignedLength - new_cigar[0].length;

        //trim sequence
        rec.sequence = rec.sequence[art_cigar.alignedLength .. $];
        rec.qscores(qual[art_cigar.alignedLength .. $].dup);

        auto len_to_clip = art_cigar.alignedLength;

        while (len_to_clip > 0)
        {
            if (new_cigar[0].op.isQueryConsuming)
            {
                if (len_to_clip < new_cigar[0].length)
                {
                    new_cigar[0].length = new_cigar[0].length - len_to_clip;
                    len_to_clip = 0;
                }
                else
                {
                    len_to_clip -= new_cigar[0].length;
                    new_cigar = new_cigar[1 .. $];
                }
            }
            else
            {
                new_cigar = new_cigar[1 .. $];
            }
        }
        new_cigar = CigarOp(art_cigar.alignedLength, Ops.HARD_CLIP) ~ new_cigar[];
    }
    if (status.art_right)
    {
        //get artifact cigar
        auto art_cigar = cigarFromString((*rec)["am"].toString.splitter(";")
                .drop(1).front.splitter(",").drop(2).front);

        //assert right side is soft-clipped 
        if (art_cigar[$ - 1].op != Ops.SOFT_CLIP)
        {
            writeln((*rec)["am"]);
            debug assert(false);
        else return;
        }
        if (art_cigar.alignedLength > rec.length)
        {
            writeln((*rec)["am"]);
            debug assert(false);
        else return;
        }

        //trim sequence
        rec.sequence = rec.sequence[0 .. $ - art_cigar.alignedLength];
        rec.qscores(qual[0 .. $ - art_cigar.alignedLength].dup);

        if (new_cigar[$ - 1].op == Ops.HARD_CLIP)
            new_cigar = new_cigar[0 .. $ - 1];
        assert(new_cigar[$ - 1].op == Ops.SOFT_CLIP);
        auto len_to_clip = art_cigar.alignedLength;

        while (len_to_clip > 0)
        {
            if (new_cigar[$ - 1].op.isQueryConsuming)
            {
                if (len_to_clip < new_cigar[$ - 1].length)
                {
                    new_cigar[$ - 1].length = new_cigar[$ - 1].length - len_to_clip;
                    len_to_clip = 0;
                }
                else
                {
                    len_to_clip -= new_cigar[$ - 1].length;
                    new_cigar = new_cigar[0 .. $ - 1];
                }
            }
            else
            {
                new_cigar = new_cigar[0 .. $ - 1];
            }
        }
        new_cigar = new_cigar[] ~ CigarOp(art_cigar.alignedLength, Ops.HARD_CLIP);
    }
    rec.cigar = new_cigar;
}

SAMRecord makeArtifactRecord(SAMRecord* original, bool left, bool mate)
{
    auto rec = SAMRecord(bam_dup1(original.b), original.h);
    rec.sequence = reverse_complement_sam_record(rec);
    if (left)
    {
        rec.cigar = cigarFromString(rec["am"].toString.splitter(";")
                .front.splitter(",").drop(2).front);
    }
    else
    {
        rec.cigar = cigarFromString(rec["am"].toString.splitter(";").drop(1)
                .front.splitter(",").drop(2).front);
    }
    //set unpaired, change strand, and supplementary
    rec.b.core.flag &= 0b1111_1111_0011_1100;
    rec.b.core.flag ^= 0b0000_0000_0001_0000;
    rec.b.core.flag |= 0b0000_1000_0000_0000;
    if (mate)
    {
        rec.tid = rec.mateTID;
    }
    if (left)
    {
        rec.pos = ZB(rec["am"].to!string.splitter(";").front.splitter(",").drop(1).front.to!long);
    }
    else
    {
        rec.pos = ZB(rec["am"].to!string.splitter(";").drop(1)
            .front.splitter(",").drop(1).front.to!long);
    }
    return rec;
}

int numericallyAwareStringComparison(string a, string b)
{
    while(a.length > 0 && b.length > 0)
    {
        if((a[0] > '9' || a[0] < '0') && (b[0] > '9' || b[0] < '0')) {
            if(a[0] == b[0]) {
                a = a[1..$];
                b = b[1..$];
                continue;
            }else if(a[0] < b[0])
                return -1;
            else return 1;
        } else {
            long a_int = -1;
            long b_int = -1;
            try {
                a_int = parse!long(a);
            } catch(ConvException)
            {

            }

            try {
                b_int = parse!long(b);
            } catch (ConvException) {
                
            }
            if(a_int == b_int) {
                continue;
            }else if(a_int < b_int) return -1;
            else return 1;
        }
    }
    if(a.length == b.length){
        return 0;
    } else if(a.length < b.length) {
        return -1;
    } else return 1;
} 

int filter(bool clip)(string cl, string[] args, ubyte con)
{
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
    
    Stats stats;
    static if (clip == true)
    {
        hts_log_warning("fade-out","Using the -c flag means the output SAM/BAM will not be sorted (regardless of prior sorting)");
        foreach (SAMRecord rec; bam.allRecords())
        {
            stats.read_count++;
            ReadStatus val;
            auto tag = rec["rs"];
            if (!tag.exists)
            {
                out_bam.write(rec);
                continue;
            }
            val.raw = tag.to!ubyte;
            stats.parse(val);
            if (!(val.art_left | val.art_right))
            {
                out_bam.write(rec);
            }
            else
            {
                clipRead(&rec, &val);
                out_bam.write(rec);
            }
        }
    }
    static if (clip == false)
    {
        import util : chunkBy;
        auto recordRange = bam.allRecords;
        auto origRange = recordRange.save;
        
        alias sort_pred = (a,b) => a.queryName.idup.numericallyAwareStringComparison(b.queryName.idup) < 0;
        
        if(recordRange.take(10).isSorted!sort_pred){
            hts_log_warning("fade-out","Output looks name-sorted, ejecting all reads with same readname if any have an artifact");
            foreach (recs; bam.allRecords.chunkBy!((a, b) => a.queryName == b.queryName))
            {
                auto group = recs.array;
                bool art_found = false;
                foreach (rec; group)
                {
                    stats.read_count++;
                    ReadStatus val;
                    auto tag = rec["rs"];
                    if (!tag.exists)
                    {
                        continue;
                    }
                    val.raw = tag.to!ubyte;
                    stats.parse(val);
                    if (val.art_left | val.art_right)
                    {
                        art_found = true;
                    }
                }
                if (!art_found)
                {
                    foreach (rec; group)
                    {
                        out_bam.write(rec);
                    }
                }
            }
        }else{
            hts_log_warning("fade-out","Output doesn't look name-sorted, ejecting by only reads with an artifact");
            foreach (rec; origRange)
            {
                stats.read_count++;
                ReadStatus val;
                auto tag = rec["rs"];
                if (!tag.exists)
                {
                    continue;
                }
                val.raw = tag.to!ubyte;
                stats.parse(val);
                if (!(val.art_left | val.art_right))
                {
                    out_bam.write(rec);
                }
            }
        }
    }
    stats.print;
    return 0;
}

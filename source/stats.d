module stats;
import std.stdio;
import std.algorithm:filter;
import std.array:array;
import std.range:drop;
import std.algorithm.iteration:splitter;
import std.conv:to;
import readstatus;
import dhtslib;
import dparasail;
import std.math: round;
import std.algorithm: sum;
import std.utf;

struct Stats {
    int read_count;
    int clipped;
    int sup;
    int art_sup;
    int art;
    int qual;
    int art_mate;
    int art_short;

    int aln_l;
    int aln_r;
    //0 Read is Softclipped
    // sc
    //1 Read has Supp Alignment
    // sup
    //2   Supp is on opposite strand from read
    // sup_opp_strand
    //3   Sc doesn't meet qual cutoff
    // qual
    //4   Read is artifact
    // art
    //5 Artifact aligns to mate region and not read
    // art_mate
    //6   Artifact is greater than 5 bp long but shorter than 15 (TODO: set empirically)
    // art_short
    //7 supp alignment not close to read or mate
    // far
    void parse(ReadStatus rs){
        clipped+=rs.sc;
        art+=(rs.art_left|rs.art_right);
        sup+=rs.sup;
        art_sup+=(rs.art_left|rs.art_right) & rs.sup;
        art_mate+=((rs.art_left & rs.mate_left)|(rs.art_right & rs.mate_right));
        aln_l+=rs.art_left;
        aln_r+=rs.art_right;
    }
    void print(){
        stderr.write("read count:\t");
        stderr.writeln(read_count);
        stderr.write("Clipped %:\t");
        stderr.writeln(clipped/float(read_count));
        stderr.write("% With Supplementary alns:\t");
        stderr.writeln(sup/float(read_count));
        stderr.write("Artifact rate:\t");
        stderr.writeln(art/float(read_count));
        stderr.write("% With Supplementary alns and artifacts:\t");
        stderr.writeln(art_sup/float(read_count));
        stderr.write("Artifact rate left only:\t");
        stderr.writeln(aln_l/float(read_count));
        stderr.write("Artifact rate right only:\t");
        stderr.writeln(aln_r/float(read_count));
    }
}

void statsfile(string[] args){
    import std.array:join,split;
    import std.format:format;
    auto bam = SAMReader(args[1]);
    auto outfile = File(args[2],"w");
    auto header = ["qname","rname","pos","cigar","art_start","art_end","aln_rname","aln_start","aln_end","art_cigar","stemloop","stemloop_rc","predicted_inverted_repeat","IR_identity","avgbq","art_avgbq","art_bq","IR_bq","flagbinary","flag","strand"];
    outfile.writeln(header.join('\t'));
    auto ps = Parasail("ACTGN",2,-1,10,2);
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
            auto am_fields = am_split[0].split(',');
            string[] towrite;
            towrite~=rec.queryName.idup;
            towrite~=bam.target_names[rec.tid];
            towrite~=rec.pos.to!string;
            towrite~=rec.cigar.toString;
            towrite~=(rec.pos-rec.cigar.ops.filter!(x=>x.op==Ops.SOFT_CLIP).front.length).to!string;
            towrite~=rec.pos.to!string;
            towrite~=am_fields[0];
            towrite~=am_fields[1];
            towrite~=(am_fields[1].to!int+cigarFromString(am_fields[2]).ref_bases_covered).to!string;
            towrite~=am_fields[2];
            towrite~=rec["as"].toString.splitter(";").front;
            towrite~=rec["ar"].toString.splitter(";").front;

            auto end = round(float(towrite[$-2].length) * 0.75).to!ulong;
            parasail_query p;
            p.seq1=toUTFz!(char *)(towrite[$-2][0..end]);
            p.seq2=toUTFz!(char *)(towrite[$-1]);
            p.seq1Len=cast(int)end;
            p.seq2Len=cast(int)towrite[$-1].length;
            p.result=ps.aligner!("sw","stats","striped","16")(p);
            towrite~=towrite[$-2][0 .. p.result.end_query + 1];
            towrite~=(float(parasail_result_get_matches(p.result)) / float(p.result.end_query + 1)).to!string;

            towrite~=(float(rec.qscores!false().sum) / float(rec.length)).to!string;
            auto bq = rec["ab"].toString[0..towrite[$-4].length].dup;
            towrite~=(float(bq.sum) / float(towrite[$-4].length)).to!string;
            bq[] = bq[]+33;
            towrite~=bq.idup;
            towrite~=towrite[$-1][0 .. p.result.end_query + 1];

            towrite~= format!"%08b"(rs.raw);
            towrite~= rs.raw.to!string;
            // this is reverse of the normal strand function as the artifact
            // is always on the opposite strand as the read
            towrite~= ["-","+"][rec.isReversed()];
            outfile.writeln(towrite.join("\t"));
        }
        if(rs.art_right){
            auto am_fields = am_split[1].split(',');
            string[] towrite;
            towrite~=rec.queryName.idup;
            towrite~=bam.target_names[rec.tid];
            towrite~=rec.pos.to!string;
            towrite~=rec.cigar.toString;
            towrite~=(rec.pos+rec.cigar.ref_bases_covered).to!string;
            towrite~=(rec.pos+rec.cigar.ref_bases_covered+rec.cigar.ops.filter!(x=>x.op==Ops.SOFT_CLIP).array[$-1].length).to!string;
            towrite~=am_fields[0];
            towrite~=am_fields[1];
            towrite~=(am_fields[1].to!int+cigarFromString(am_fields[2]).ref_bases_covered).to!string;
            towrite~=am_fields[2];
            towrite~=rec["as"].toString.splitter(";").drop(1).front;
            towrite~=rec["ar"].toString.splitter(";").drop(1).front;
           
            auto end = round(float(towrite[$-2].length) * 0.75).to!ulong;
            parasail_query p;
            p.seq1=toUTFz!(char *)(towrite[$-2][0..end]);
            p.seq2=toUTFz!(char *)(towrite[$-1]);
            p.seq1Len=cast(int)end;
            p.seq2Len=cast(int)towrite[$-1].length;
            p.result=ps.aligner!("sw","stats","striped","16")(p);
            towrite~=towrite[$-2][0 .. p.result.end_query + 1];
            towrite~=(float(parasail_result_get_matches(p.result)) / float(p.result.end_query + 1)).to!string;

            towrite~=(float(rec.qscores!false().sum) / float(rec.length)).to!string;
            auto bq = rec["ab"].toString[$-towrite[$-4].length .. $].dup;
            towrite~=(float(bq.sum) / float(towrite[$-4].length)).to!string;
            bq[]=bq[]+33;
            towrite~=bq.idup;
            towrite~=towrite[$-1][0 .. p.result.end_query + 1];
 
            towrite~= format!"%08b"(rs.raw);
            towrite~= rs.raw.to!string;
            // this is reverse of the normal strand function as the artifact
            // is always on the opposite strand as the read
            towrite~= ["-","+"][rec.isReversed()];
            outfile.writeln(towrite.join("\t"));
        }
    }
}


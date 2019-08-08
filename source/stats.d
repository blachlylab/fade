module stats;
import std.stdio;
import std.algorithm:filter;
import std.array:array;
import std.conv:to;
import readstatus;
import dhtslib;

struct Stats {
    int read_count;
    int clipped;
    int sup;
    int sup_opp;
    int art;
    int qual;
    int art_mate;
    int art_short;

    int aln_l;
    int aln_r;
    //0	Read is Softclipped
    // sc
    //1	Read has Supp Alignment
    // sup
    //2   Supp is on opposite strand from read
    // sup_opp_strand
    //3   Sc doesn't meet qual cutoff
    // qual
    //4   Read is artifact
    // art
    //5	Artifact aligns to mate region and not read
    // art_mate
    //6   Artifact is greater than 5 bp long but shorter than 15 (TODO: set empirically)
    // art_short
    //7 supp alignment not close to read or mate
    // far
    void parse(ReadStatus rs){
        clipped+=(rs.left.sc|rs.right.sc);
        sup+=(rs.left.sup|rs.right.sup);
        sup_opp+=(rs.left.sup_opp_strand|rs.right.sup_opp_strand);
        qual+=(rs.left.qual|rs.right.qual);
        art+=(rs.left.art|rs.right.art);
        art_mate+=(rs.left.art_mate|rs.right.art_mate);
        art_short+=(rs.left.art_short|rs.right.art_short);

        aln_l+=rs.left.art;
        aln_r+=rs.right.art;
    }
    void print(){
        stderr.write("read count:\t");
        stderr.writeln(read_count);
        stderr.write("Clipped %:\t");
        stderr.writeln(clipped/float(read_count));
        stderr.write("With Supplementary alns:\t");
        stderr.writeln(sup/float(read_count));
        stderr.write("With Supplementary alns on opposite strand:\t");
        stderr.writeln(sup_opp/float(read_count));
        stderr.write("Artifact rate:\t");
        stderr.writeln(art/float(read_count));
        stderr.write("Artifact rate left only:\t");
        stderr.writeln(aln_l/float(read_count));
        stderr.write("Artifact rate right only:\t");
        stderr.writeln(aln_r/float(read_count));
        stderr.write("Artifact rate short (<15bp):\t");
        stderr.writeln(art_short/float(read_count));
        stderr.write("Artifact rate mate:\t");
        stderr.writeln(art_mate/float(read_count));
    }
}

void statsfile(string[] args){
    import std.array:join,split;
    import std.format:format;
    auto bam = SAMReader(args[1]);
    auto outfile = File(args[2],"w");
    auto header = ["qname","rname","pos","cigar","art_start","art_end","aln_rname","aln_start","aln_end","flagbinary","flag"];
    outfile.writeln(header.join('\t'));
    foreach(SAMRecord rec;bam.all_records()){
        auto tag=rec["rs"];
        if(tag.data==null) continue;
        ReadStatus rs;
        rs.raw=tag.to!ushort;
        if(!(rs.left.art|rs.right.art)) continue;
        tag=rec["am"];
        if(tag.data==null) continue;
        auto am = tag.toString;
        auto am_split = am.split(';');
        if(rs.left.art){
            auto am_fields = am_split[0].split(',');
            string[] towrite;
            towrite~=rec.queryName.idup;
            towrite~=bam.target_names[rec.tid];
            towrite~=rec.pos.to!string;
            towrite~=rec.cigar.toString;
            towrite~=(rec.pos-rec.cigar.ops.filter!(x=>x.op==Ops.SOFT_CLIP).array[0].length).to!string;
            towrite~=rec.pos.to!string;
            towrite~=am_fields[0];
            towrite~=am_fields[1];
            towrite~=am_fields[2];
            towrite~= format!"%08b"(rs.left.raw);
            towrite~= rs.left.raw.to!string;
            outfile.writeln(towrite.join("\t"));
        }
        if(rs.right.art){
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
            towrite~=am_fields[2];
            towrite~= format!"%08b"(rs.right.raw);
            towrite~= rs.right.raw.to!string;
            outfile.writeln(towrite.join("\t"));
        }
    }
}


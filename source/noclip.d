module noclip;
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
import util;

void noclipfile(string[] args){
    import std.array:join,split;
    import std.format:format;
    auto bam = SAMReader(args[1]);
    auto outfile = File(args[2],"w");
    auto header = ["qname","sc_q_scores","sc_seq","sc_avg_bq","avg_bq","art_status"];
    outfile.writeln(header.join('\t'));
    auto ps = Parasail("ACTGN",3,-8,10,5);
    foreach(SAMRecord rec;bam.allRecords()){
        auto tag=rec["rs"];
        if(!tag.exists) continue;
        ReadStatus rs;
        rs.raw=tag.to!ubyte;
        if(!rs.sc) continue;
        auto clips = parse_clips(rec.cigar);
        if(clips[0].raw){
            string[] towrite;
            towrite~=rec.queryName.idup;
            auto q_scores = rec.qscores();
            auto clip_len = clips[0].length;
            auto sc_q_scores = q_scores[0..clip_len].dup;
            auto sc_avg_bq = (float(sc_q_scores.sum) / float(clip_len));
            auto avg_bq = (float(q_scores.sum) / float(rec.length));
            sc_q_scores[] = sc_q_scores[]+33;
            towrite~=(cast(char[])sc_q_scores).idup;
            towrite~=rec.sequence[0..clip_len].idup;
            towrite~=sc_avg_bq.to!string;
            towrite~=avg_bq.to!string;
            towrite~=rs.art_left.to!string;
            outfile.writeln(towrite.join("\t"));
        }
        if(clips[1].raw){
            string[] towrite;
            towrite~=rec.queryName.idup;
            auto q_scores = rec.qscores();
            auto clip_len = clips[1].length;
            auto sc_q_scores = q_scores[$-clip_len..$].dup;
            auto sc_avg_bq = (float(sc_q_scores.sum) / float(clip_len));
            auto avg_bq = (float(q_scores.sum) / float(rec.length));
            sc_q_scores[] = sc_q_scores[]+33;
            towrite~=(cast(char[])sc_q_scores).idup;
            towrite~=rec.sequence[$-clip_len..$].idup;
            towrite~=sc_avg_bq.to!string;
            towrite~=avg_bq.to!string;
            towrite~=rs.art_right.to!string;
            outfile.writeln(towrite.join("\t"));
        }
    }
}


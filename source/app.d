module main;
import std.stdio;
import std.getopt;
import std.parallelism:defaultPoolThreads;
import filter:filter;
import anno;
import stats;
import remap;
import noclip;



int artifact_floor_length=5;
int align_buffer_size=300;
int threads;
bool clip = false;
bool output_bam;
bool output_ubam;

string full_help = "Fragmentase Artifact Detection and Elimination\n"~
                "usage: ./fade [subcommand]\n"~
                "\tannotate: marks artifact reads in bam tags (must be done first)\n"~
                "\tout: eliminates artifact from reads(may require queryname sorted bam)\n"~
                "\tstats: reports extended information about artifact reads\n"~
                "\tstats-clip: reports extended information about all soft-clipped reads\n"~
                "\textract: extracts artifacts into a mapped bam";

void main(string[] args){
    if(args.length==1){
        auto res=getopt(args,config.bundling);
        defaultGetoptPrinter(
            full_help, res.options);
        stderr.writeln();
        return;
    }
    if(args[1]=="annotate"){
        auto res=getopt(args,config.bundling,
            "threads|t","extra threads for parsing the bam file",&threads,
            "min-length","Minimum number of bases for a soft-clip to be considered for artifact detection",&artifact_floor_length,
            "window-size|w","Number of bases considered outside of read or mate region for re-alignment",&align_buffer_size,
            "bam|b","output bam",&output_bam,
            "ubam|u","output uncompressed bam",&output_ubam);
        if (res.helpWanted | (args.length < 3)) {
            defaultGetoptPrinter(
                "Fragmentase Artifact Detection and Elimination\n"~
                "annotate: performs re-alignment of soft-clips and annotates bam records with bitflag (rs) and realignment tags (am)\n"~
                "usage: ./fade annotate [BAM/SAM input] [Indexed fasta reference]\n", res.options);
            stderr.writeln();
            return;
        }
        if(threads!=0){
            defaultPoolThreads(threads);
        }
        ubyte con = (((cast(ubyte) output_bam)<<1) | cast(ubyte) output_ubam);
        if(con>2) {
            stderr.writeln("please use only one of the b or u flags");
            return;
        }
        annotate(args[1..$],con);
    }else if(args[1]=="out"){
        auto res=getopt(args,config.bundling,
            "clip|c","clip reads instead of filtering them",&clip,
            "threads|t","extra threads for parsing the bam file",&threads,
            "bam|b","output bam",&output_bam,
            "ubam|u","output uncompressed bam",&output_ubam
        );

        if (res.helpWanted | (args.length < 3)) {
            defaultGetoptPrinter(
                "Fragmentase Artifact Detection and Elimination\n"~
                "out: removes all read and mates for reads contain the artifact (used after annotate and requires queryname sorted bam)"~
                " or, with the -c flag, hard clips out artifact sequence from reads\n"~
                "usage: ./fade out [BAM/SAM input]\n", res.options);
            stderr.writeln();
            return;
        }
        if(threads!=0){
            defaultPoolThreads(threads);
        }
        ubyte con = (((cast(ubyte) output_bam)<<1) | cast(ubyte) output_ubam);
        if(con>2) {
            stderr.writeln("please use only one of the b or u flags");
            return;
        }
        if(clip)
            filter!(true)(args[1..$],con);
        else
            filter!(false)(args[1..$],con);
    }else if(args[1]=="extract"){
        auto res=getopt(args,config.bundling,
            "threads|t","extra threads for parsing the bam file",&threads,
            "bam|b","output bam",&output_bam,
            "ubam|u","output uncompressed bam",&output_ubam
        );

        if (res.helpWanted | (args.length < 3)) {
            defaultGetoptPrinter(
                "Fragmentase Artifact Detection and Elimination\n"~
                "extract: extracts artifacts into a mapped bam\n"~
                "usage: ./fade extract [BAM/SAM input]\n", res.options);
            stderr.writeln();
            return;
        }
        if(threads!=0){
            defaultPoolThreads(threads);
        }
        ubyte con = (((cast(ubyte) output_bam)<<1) | cast(ubyte) output_ubam);
        if(con>2) {
            stderr.writeln("please use only one of the b or u flags");
            return;
        }
        remapArtifacts(args[1..$],con);
    }else if(args[1]=="stats"){
        auto res=getopt(args,config.bundling,
        "threads|t","threads for parsing the bam file",&threads);
        if (res.helpWanted | (args.length < 3)) {
            defaultGetoptPrinter(
                "Fragmentase Artifact Detection and Elimination\n"~
                "stats: reports extended information about artifact reads (used after annotate)", res.options);
            stderr.writeln();
            return;
        }
        if(threads!=0){
            defaultPoolThreads(threads);
        }
        statsfile(args[1..$]);
    }else if(args[1]=="stats-clip"){
        auto res=getopt(args,config.bundling,
        "threads|t","threads for parsing the bam file",&threads);
        if (res.helpWanted | (args.length < 3)) {
            defaultGetoptPrinter(
                "Fragmentase Artifact Detection and Elimination\n"~
                "stats-clip: reports extended information about all soft-clipped reads (used after annotate)", res.options);
            stderr.writeln();
            return;
        }
        if(threads!=0){
            defaultPoolThreads(threads);
        }
        noclipfile(args[1..$]);
    }else{
        auto res=getopt(args,config.bundling);
        if (res.helpWanted | (args.length < 2)) {
            defaultGetoptPrinter(
                full_help, res.options);
            stderr.writeln();
            return;
        }
    }
}




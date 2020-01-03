module main;
import std.stdio;
import std.getopt;
import std.parallelism:defaultPoolThreads;
import filter:filter;
import anno;
import stats;



int artifact_floor_length=5;
int artifact_short_cutoff=15;
int align_buffer_size=300;
int mate_size_est=151;
int qscore_cutoff=20;
int threads;
bool clip = false;
string artifact_bam="";

void main(string[] args){
    if(args.length==1){
        auto res=getopt(args,config.bundling);
        defaultGetoptPrinter(
            "Fragmentase Artifact Detection and Elimination\n"~
            "usage: ./fade [subcommand]\n"~
            "\tannotate: marks artifact reads in bam tags (must be done first)\n"~
            "\tout: eliminates artifact from reads(may require queryname sorted bam)", res.options);
        stderr.writeln();
        return;
    }
    if(args[1]=="annotate"){
        auto res=getopt(args,config.bundling,
            "threads|t","threads for parsing the bam file",&threads,
            "min-length","Minimum number of bases for a soft-clip to be considered for artifact detection",&artifact_floor_length,
            "short-length","Minimum number of bases for a soft-clip to not be considered short",&artifact_short_cutoff,
            "window-size|w","Number of bases considered outside of read or mate region for re-alignment",&align_buffer_size,
            "q-score|q","Minimum average base-quality score of a soft-clip to be considered an artifact",&qscore_cutoff,
            "mate-est|m","Read Mate size estimate in bases",&mate_size_est);
        if (res.helpWanted | (args.length < 3)) {
            defaultGetoptPrinter(
                "Fragmentase Artifact Detection and Elimination\n"~
                "annotate: performs re-alignment of soft-clips and annotates bam records with bitflag (rs) and realignment tags (am)", res.options);
            stderr.writeln();
            return;
        }
        if(threads!=0){
            defaultPoolThreads(threads);
        }
        annotate(args[1..$]);
    }else if(args[1]=="out"){
        auto res=getopt(args,config.bundling,
            "clip|c","clip reads instead of filtering them",&clip,
            "threads|t","threads for parsing the bam file",&threads,
            "artifact-bam|a","filename to extract artifact reads to (BAM/SAM)",&artifact_bam
        );

        if (res.helpWanted | (args.length < 3)) {
            defaultGetoptPrinter(
                "Fragmentase Artifact Detection and Elimination\n"~
                "out: removes all read and mates for reads contain the artifact (used after annotate and requires queryname sorted bam)"~
                " or, with the -c flag, hard clips out artifact sequence from reads", res.options);
            stderr.writeln();
            return;
        }
        if(threads!=0){
            defaultPoolThreads(threads);
        }
        if(clip && artifact_bam!="")
            filter!(true,true)(args[1..$],artifact_bam);
        else if(clip && artifact_bam=="")
            filter!(true,false)(args[1..$],artifact_bam);
        else if(!clip && artifact_bam!="")
            filter!(false,true)(args[1..$],artifact_bam);
        else
            filter!(false,false)(args[1..$],artifact_bam);
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
    }else{
        auto res=getopt(args,config.bundling);
        if (res.helpWanted | (args.length < 2)) {
            defaultGetoptPrinter(
                "Fragmentase Artifact Detection and Elimination\n"~
                "usage: ./fade [subcommand]\n"~
                "\tannotate: marks artifact reads in bam tags (must be done first)\n"~
                "\tout: eliminates artifact from reads(may require queryname sorted bam)", res.options);
            stderr.writeln();
            return;
        }
    }
}




module main;
import std.stdio;
import std.getopt;
import std.parallelism : defaultPoolThreads;
import filter : filter;
import std.array : join;
import std.format : format;
import dhtslib : SAMReader;
import htslib.hts_log;
import anno;
import stats;
import remap;
import noclip;
import _version;


int artifact_floor_length = 5;
int align_buffer_size = 300;
int threads;
bool clip = false;
bool output_bam;
bool output_ubam;

enum fade_header = "Fragmentase Artifact Detection and Elimination
version: %s
".format(VERSION);

string full_help = fade_header ~ "
usage: fade [subcommand]
    annotate: marks artifact reads in bam tags (must be done first)
    out: eliminates artifact from reads(may require queryname sorted bam)
    stats: reports extended information about artifact reads
    stats-clip: reports extended information about all soft-clipped reads
    extract: extracts artifacts into a mapped bam
";

string anno_help = fade_header ~ "
annotate: performs re-alignment of soft-clips and annotates bam records with bitflag (rs) and realignment tags (am)
usage: fade annotate [options] <input BAM/SAM> <Indexed fasta reference>
";

string out_help = fade_header ~ "
out: removes all read and mates for reads contain the artifact (used after annotate)
     or, with the -c flag, hard clips out artifact sequence from reads
     it is reccomended that the input SAM/BAM be queryname sorted
usage: fade out [options] <input BAM/SAM>
";

string extract_help = fade_header ~ "
extract: extracts artifacts into a mapped SAM/BAM (used after annotate)
usage: ./fade extract [options] <annotated BAM/SAM> 
";

string stats_help = fade_header ~ "
stats: reports extended information about all artifact reads (used after annotate)
usage: ./fade stats [options] <annotated BAM/SAM>  
";

string stats_clip_help = fade_header ~ "
stats-clip: reports extended information about all soft-clipped reads (used after annotate)
usage: ./fade stats-clip [options] <annotated BAM/SAM> 
";

int main(string[] args)
{
    auto cl = join(args," ");
    if (args.length == 1)
    {
        auto res = getopt(args, config.bundling);
        defaultGetoptPrinter(full_help, res.options);
        stderr.writeln();
        return 0;
    }
    if (args[1] == "annotate")
    {
        auto res = getopt(args, config.bundling, "threads|t",
                "extra threads for parsing the bam file", &threads, "min-length",
                "Minimum number of bases for a soft-clip to be considered for artifact detection",
                &artifact_floor_length,
                "window-size|w",
                "Number of bases considered outside of read or mate region for re-alignment",
                &align_buffer_size, "bam|b", "output bam", &output_bam,
                "ubam|u", "output uncompressed bam", &output_ubam);
        if (res.helpWanted | (args.length < 3))
        {
            defaultGetoptPrinter(anno_help,res.options);
            stderr.writeln();
            return 0;
        }
        if (threads != 0)
        {
            defaultPoolThreads(threads);
        }
        ubyte con = (((cast(ubyte) output_bam) << 1) | cast(ubyte) output_ubam);
        if (con > 2)
        {
            hts_log_error("fade-annotate", "Please use only one of the b or u flags");
            return 1;
        }
        return annotate(cl, args[1 .. $], con, artifact_floor_length, align_buffer_size);
    }
    else if (args[1] == "out")
    {
        auto res = getopt(args, config.bundling, "clip|c", "clip reads instead of filtering them",
                &clip, "threads|t", "extra threads for parsing the bam file",
                &threads, "bam|b", "output bam", &output_bam, "ubam|u",
                "output uncompressed bam", &output_ubam);

        if (res.helpWanted | (args.length < 3))
        {
            defaultGetoptPrinter(out_help,res.options);
            stderr.writeln();
            return 0;
        }
        if (threads != 0)
        {
            defaultPoolThreads(threads);
        }
        ubyte con = (((cast(ubyte) output_bam) << 1) | cast(ubyte) output_ubam);
        if (con > 2)
        {
            hts_log_error("fade-annotate", "Please use only one of the b or u flags");
            return 1;
        }
        if (clip)
            return filter!(true)(cl, args[1 .. $], con);
        else
            return filter!(false)(cl, args[1 .. $], con);
    }
    else if (args[1] == "extract")
    {
        auto res = getopt(args, config.bundling, "threads|t",
                "extra threads for parsing the bam file", &threads, "bam|b",
                "output bam", &output_bam, "ubam|u", "output uncompressed bam", &output_ubam);

        if (res.helpWanted | (args.length < 3))
        {
            defaultGetoptPrinter(extract_help, res.options);
            stderr.writeln();
            return 0;
        }
        if (threads != 0)
        {
            defaultPoolThreads(threads);
        }
        ubyte con = (((cast(ubyte) output_bam) << 1) | cast(ubyte) output_ubam);
        if (con > 2)
        {
            hts_log_error("fade-annotate", "Please use only one of the b or u flags");
            return 1;
        }
        return remapArtifacts(cl, args[1 .. $], con);
    }
    else if (args[1] == "stats")
    {
        auto res = getopt(args, config.bundling, "threads|t",
                "threads for parsing the bam file", &threads);
        if (res.helpWanted | (args.length < 3))
        {
            defaultGetoptPrinter(stats_help,res.options);
            return 0;
        }
        if (threads != 0)
        {
            defaultPoolThreads(threads);
        }

        File outfile;
        if(args.length == 3)
            outfile = File(args[2], "w");
        else if(args.length == 2)
            outfile = stdout;
        else{
            defaultGetoptPrinter(stats_help,res.options);
            return 1;
        }

        auto bam = SAMReader(args[1]);

        return statsfile(bam, outfile);
    }
    else if (args[1] == "stats-clip")
    {
        auto res = getopt(args, config.bundling, "threads|t",
                "threads for parsing the bam file", &threads);
        if (res.helpWanted | (args.length < 3))
        {
            defaultGetoptPrinter(stats_clip_help,res.options);
            return 0;
        }
        if (threads != 0)
        {
            defaultPoolThreads(threads);
        }
        File outfile;
        if(args.length == 3)
            outfile = File(args[2], "w");
        else if(args.length == 2)
            outfile = stdout;
        else{
            defaultGetoptPrinter(stats_clip_help,res.options);
            return 1;
        }

        auto bam = SAMReader(args[1]);

        return noclipfile(bam, outfile);
    }
    else
    {   
        auto res = getopt(args, config.bundling);
        if (res.helpWanted | (args.length < 2))
        {
            defaultGetoptPrinter(full_help, res.options);
            return 0;
        }else{
            hts_log_error("fade", args[1] ~ " is not a fade subcommand");
            defaultGetoptPrinter(full_help, res.options);
            return 1;
        }
    }
}

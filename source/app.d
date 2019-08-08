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
int sa_size_wiggle=5;
int threads;

void main(string[] args){
    auto res=getopt(args,config.bundling,
	"threads|t","threads for parsing the bam file",&threads);
	if (res.helpWanted) {
		defaultGetoptPrinter(
            "usage: ./fade [annotate] [bam] [reference fasta with fai] [out bam]\n"~
            "usage: ./fade [filter or clip] [bam] [out bam]\n"~
            "annotate: marks artifact reads in bam tags (must be done first)\n"~
            "filter: removes fragments (read and mate) with artifact (requires queryname sorted bam)\n"~
            "clip: removes artifact region only", res.options);
		stderr.writeln();
		return;
	}
	if(args.length<3){
		writeln("usage: ./fade [annotate] [bam] [reference fasta with fai] [out bam]\nusage: ./fade [filter or clip] [bam] [out bam]");
		return;
	}else{
		if(threads!=0){
			defaultPoolThreads(threads);
		}
    }
    if(args[1]=="annotate"){
        annotate(args[1..$]);
    }else if(args[1]=="filter"){
        filter!false(args[1..$]);
    }else if(args[1]=="clip"){
        filter!true(args[1..$]);
    }else if(args[1]=="stats"){
        statsfile(args[1..$]);
    }
}




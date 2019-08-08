module anno;
import std.algorithm:filter;
import std.conv:to;
import std.algorithm.searching:count;
import std.algorithm.iteration:splitter;
import std.array:array;
import dhtslib;
import dparasail;
import readstatus;
import analysis;
import util;

void annotate(string[] args){
	auto bam = SAMReader(args[1]);
	auto fai=IndexedFastaFile(args[2]);
	auto out_bam=SAMWriter(args[3],bam.header);
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
	auto p=Parasail("ACTGN",1,-1,1,3);
	foreach(SAMRecord rec;bam.all_records){
        ReadStatus status;
        if(rec.isSupplementary()||
            rec.isSecondary()||
            !rec.isMapped()||
            rec.cigar.ops.filter!(x=>x.op==Ops.SOFT_CLIP).count()==0
        ){
            rec["rs"]=status.raw;
            out_bam.write(&rec);
            continue;
        }
        CigarOp[2] clips=parse_clips(rec.cigar);
        if(clips[0].length!=0) status.left.sc=true;
        if(clips[1].length!=0) status.right.sc=true;
		if(rec["SA"].data!=null){
			string[] sup=rec["SA"].toString.splitter(",").array;
            if(sup[2][0]!=rec.strand){
                status.left.sup_opp_strand=true;
                status.right.sup_opp_strand=true;
                auto sa_cigar =cigarFromString(sup[3]);
                if(clips[0].length!=0){
                    if(sa_cigar.ops[0].op!=Ops.SOFT_CLIP){
                        status.left.sup=true;
                        status.left.art=true;
                        rec["rs"]=status.raw;
                        rec["am"]=sup[0]~","~sup[1]~","~(sup[1].to!int+sa_cigar.ops[0].length).to!string~","~sa_cigar.ops[0].length.to!string~";";
                        out_bam.write(&rec);
                        continue;
                    }
                }
                if(clips[1].length!=0){
                    if(sa_cigar.ops[$-1].op!=Ops.SOFT_CLIP){
                        status.right.sup=true;
                        status.right.art=true;
                        rec["rs"]=status.raw;
                        rec["am"]=";"~sup[0]~","~sup[1]~","~(sup[1].to!int+sa_cigar.ops[0].length).to!string~","~sa_cigar.ops[0].length.to!string;
                        out_bam.write(&rec);
                        continue;
                    }
                }
            }
		}
		//left soft-clip (left on reference not 5' neccesarily)
        string align_str;
		if(clips[0].length!=0){
            align_str~=align_clip(&bam,&fai,&p,&rec,&status,clips[0].length(),true);
		}
        align_str~=";";
		//right soft-clip
		if(clips[1].length()!=0){
            align_str~=align_clip(&bam,&fai,&p,&rec,&status,clips[1].length(),false);
		}
        // writeln(status.raw);
        rec["rs"]=status.raw;
        rec["am"]=align_str;
        assert(rec["rs"].check!ubyte);
        assert(rec["am"].check!string);
        out_bam.write(&rec);
	}
}

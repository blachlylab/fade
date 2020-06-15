module anno;
import std.algorithm:filter;
import std.conv:to;
import std.algorithm.searching:count;
import std.algorithm.iteration:splitter;
import std.array:array;
import core.sync.mutex:Mutex;
import std.parallelism:parallel;
import dhtslib;
import dparasail;
import readstatus;
import analysis;
import util;

void annotate(string[] args,ubyte con){
    auto bam = SAMReader(args[1]);
    // auto fai=IndexedFastaFile(args[2]);
    auto out_bam=getWriter(con,bam.header);
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
    auto p=Parasail("ACTGN",2,-3,10,2);
    auto m = new Mutex();
    foreach(SAMRecord rec;parallel(bam.all_records)){
        ReadStatus status;
        if(rec.isSupplementary()||
            rec.isSecondary()||
            !rec.isMapped()||
            rec.cigar.ops.filter!(x=>x.op==Ops.SOFT_CLIP).count()==0
        ){
            rec["rs"]=status.raw;
            m.lock;
            out_bam.write(&rec);
            m.unlock;
            continue;
        }
        CigarOp[2] clips=parse_clips(rec.cigar);
        if(clips[0].length!=0||clips[1].length!=0) status.sc=true;
        if(rec["SA"].exists) status.sup=true;
        //left soft-clip (left on reference not 5' neccesarily)
        Align_Result align_1,align_2;
        if(clips[0].length!=0){
            align_1=align_clip!true(&bam,args[2],&p,&rec,&status,clips[0].length());
        }
        //right soft-clip
        if(clips[1].length()!=0){
            align_2=align_clip!false(&bam,args[2],&p,&rec,&status,clips[1].length());
        }
        // writeln(status.raw);
        rec["rs"]=status.raw;
        rec["am"]=align_1.alignment~";"~align_2.alignment;
        rec["as"]=align_1.stem_loop~";"~align_2.stem_loop;
        rec["ar"]=align_1.stem_loop_rc~";"~align_2.stem_loop_rc;
        rec["ap"]=align_1.palindrome~";"~align_2.palindrome;
        assert(rec["rs"].check!ubyte);
        assert(rec["am"].check!string);
        m.lock;
        out_bam.write(&rec);
        m.unlock;
    }
    out_bam.close;
}

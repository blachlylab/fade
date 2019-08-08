module readstatus;
import std.bitmanip:bitfields;

/// ushort bitflag for clip status
union ReadStatus{
    ushort raw;
    struct{
        ClipStatus left;
        ClipStatus right;
    }
}
/// ubyte bitflag for indicating artifact status
union ClipStatus {
    /// raw ubyte
    ubyte raw;
    //ubyte read status encoding
    mixin(bitfields!(
        //0	Read is Softclipped
        bool,"sc",1,
        //1	Read has Supp Alignment
        bool,"sup",1,
        //2   Supp is on opposite strand from read
        bool,"sup_opp_strand",1,
        //3   Sc doesn't meet qual cutoff
        bool,"qual",1,
        //4   Read is artifact
        bool,"art",1,
        //5	Artifact aligns to mate region and not read
        bool,"art_mate",1,
        //6   Artifact is greater than 5 bp long but shorter than 15 (TODO: set empirically)
        bool,"art_short",1,
        //7 supp alignment not close to read or mate
        bool,"far",1,
    ));
}
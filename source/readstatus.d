module readstatus;
import std.bitmanip : bitfields;

/// ushort bitflag for clip status
union ReadStatus
{
    ubyte raw;
    mixin(bitfields!(
        /// 0	Read is Softclipped
        bool, "sc", 1,
        /// 1	Read's left sc is an artifact
        bool, "art_left", 1,
        /// 2   Read's right sc is an artifact
        bool, "art_right", 1,
        /// 3   Left Artifact aligns to mate region and not read
        bool, "mate_left", 1,
        /// 4   Right Artifact aligns to mate region and not read
        bool, "mate_right", 1,
        /// 5	read has supplementary alignment
        bool, "sup", 1,
        /// filler for ubyte size
        bool,"f1", 1,
        /// filler for ubyte size
        bool, "f2", 1,
    ));
}

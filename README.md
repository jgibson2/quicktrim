# quicktrim
Short Read trimming application written in C.

This application was created in an attempt to speed up the process of adapter trimming in Illumina and other types of short reads. Typically, the process of trimming is unnecessarily time-consuming, oftentimes taking up to or longer than the actual alignment step. Although there are many feature-rich tools for adapter trimming, they are generally inefficient due to language choice or other design choices.

## Design
quicktrim uses aligned memory to increase reading speed. In addition, quicktrim leverages OpenMP for parallelization, and includes other optimizations for fast reading and trimming of sequences. quicktrim utilizes a minimal number of dependencies, depending only the C standard library (note: for aligned memory to work, the C standard lib must be C11 or more recent or _POSIX_VERSION must be at least 200112) and OpenMP.

## Options
        o: Output base filename (default out)
        i: Single-ended input file
        1: Paired-end input file 1
        2: Paired-end input file 2
        a: 3' adapter
        A: 5' adapter
        q: Set quality cutoff (default 30)
        l: Set length cutoff (default 20)
        r: Set number of high-quality bases in a row for use with fast trimming (default 5)
        p: Phred score base (default 33)
        b: Change buffer size in bytes (default 100000) {NOTE: must be larger than record length / 10}
        v: Minimum 3' adapter overlap (default 10)
        s: Minimum 3' alignment score (default 8)
        V: Minimum 5' adapter overlap (default 10)
        S: Minimum 5' alignment score (default 8)
        m: Switch to faster, but less accurate mode of triming (HQ bases in a row)
        h: Print help menu

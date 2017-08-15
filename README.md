# quicktrim
Short Read trimming application written in C.

This application was created in an attempt to speed up the process of adapter trimming in Illumina and other types of short reads. Typically, the process of trimming is unnecessarily time-consuming, oftentimes taking up to or longer than the actual alignment step. Although there are many feature-rich tools for adapter trimming, they are generally inefficient due to language choice or other design choices.

## Design
quicktrim uses aligned memory to increase reading speed. In addition, quicktrim leverages OpenMP for parallelization, and includes other optimizations for fast reading and trimming of sequences.

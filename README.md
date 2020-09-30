The GDPM algorithm for computing closed itemsets w.r.t. closure levels. 

Copyright (c) 2020 Tatiana Makhalova



GDPM mines FIMI .dat files for closed itemsets. The output of the algorithm is the closure levels of closed itemsets.
The output folder "closure_structure/dat" is created in the source data directory. The folder contains the list of closed
itemsets with their supports. 


   -- How to compile --- 

Execute "make", then you can see "gdpm" in the same directory. 


   --- Command Line Options ---

GDPM requires at least file name. The parameters listed below are optional:

 E: run the GDPM-EXT version, otherwise the GDPM-INT version will be executed,
 Q: do not output the closure levels;
 -n[num]: the maximal closure level to be computed, by default all closure levels ar computed. 

   --- Example --- 

% gdmp test/example.dat

- to compute the whole closure structure by GDPM-INT and output the results to './closure_structure/dat/' by levels. The summary file 'summary_I.csv' will be written in the same folder.


% gdmp test/example.dat EQ -n 2 

- to compute the first 2 closure levels (of the closure structure) by GDPM-ENT and to do not output the closure levels to './closure_structure/dat/' by levels. The summary file 'summary_E.csv' will be written in the same folder.

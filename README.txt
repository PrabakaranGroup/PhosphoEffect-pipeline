RING-2.0 - Residue Interaction Network Generation 
--------------------------------------------------------------
Author: Damiano Piovesan
E-mail: damianopiovesan@gmail.com, silvio.tosatto@unipd.it
Date: 21st December 2015
Version: 2.0

-------------------------------------------------

0 - Remark


RING is a command line program written in C++. It is designed 
to run on POSIX systems (Linux/Unix). The C++ sources have been 
compiled with g++ 4.8.2. RING exists also as web server and it 
is freely available ar URL: 
	
	http://protein.bio.unipd.it/ring/ 


The RING method (first version) has been published in:

RING: networking interacting residues, evolutionary information and 
energetics in protein structures
Alberto J. M. Martin, Michele Vidotto, Filippo Boscariol, 
Tomàs Di Domenico, Ian Walsh and Silvio C. E. Tosatto
Bioinformatics (2011) 27 (14) 



-------------------------------------------------

1 - Installation

RING works out of the box as it is. However to run it outside the
installation path you have to set an environment variable, VICTOR_ROOT.
For example if you have RING installed in the folder /home/user/ring
you can just run:

	export VICTOR_ROOT=/home/user/ring



-------------------------------------------------

2 - Execution and examples

To execute RING just reach the installation path and type:

	./bin/Ring -i <pdb_file>


By default RING-2.0 calculate contacts of the first chain and providing
interaction types (HBONDS, VDW, IONIC, SSBOND, PIPISTACK, PICATION).
To customize your calculation and for a complete list of RING parameters
just read the program help by typing:

	./bin/Ring -h


For further details about the method, visit the RING web server ("about" 
"help" pages).

	http://protein.bio.unipd.it/ring/ 

	










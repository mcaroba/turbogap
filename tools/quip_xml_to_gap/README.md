Author(s): Mikhail S. Kuklin, Miguel A. Caro, Richard Jana

**Documentation**

_make_gap_files.py_ script transforms potential files from Quip to TurboGAP format for running MD.

make_gap_files.py [name_of_the_quip_xml] [name_of_the_output_file]

If Hirshfeld volumes GAP(s) are presented, then the command is:

make_gap_files.py [name_of_the_quip_xml] [name_of_the_output_file] [name_of_the_hirshfeld_xml]

**Example**

Without Hirshfeld volumes GAP(s):

make_gap_files.py pot.xml pot.gap

With Hirshfeld volumes GAP(s):

make_gap_files.py pot.xml pot.gap hirshfeld.xml

**Note**
- all necessary files should be in the folder called 'gap_files'
- compression files should be included to the folder 'gap_files' manually with the name compress_[number_of_many_body_descriptor].dat (for example: compress_1.dat)
- current implementation allows reading only ONE hirshfeld.xml file. However, user can create necessary files running code several times and manually changing gap file

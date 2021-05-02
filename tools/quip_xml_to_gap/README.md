Author(s): Mikhail S. Kuklin, Miguel A. Caro

Documentation

make_gap_files.py [name_of_the_quip_xml] [name_of_the_output_file]

If Hirshfeld volumes GAP(s) are presented, then the command is:

make_gap_files.py [name_of_the_quip_xml] [name_of_the_output_file] [name_of_the_hirshfeld_xml]

Example

Without Hirshfeld volumes GAP(s):

make_gap_files.py pot.xml pot.gap

With Hirshfeld volumes GAP(s):

make_gap_files.py pot.xml pot.gap hirshfeld.xml

Note that compression files should be included to the folder 'gap_files' manually with the name compress_[number_of_many_body_descriptor].dat (for example: compress_1.dat]
All necessary files should be in the folder called 'gap_files'.

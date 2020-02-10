This folder contains all of the code required to create the
GIP static file - ie, the file that handles all of the transformations
between the geographic and Apex-magnetic field frames. It does not include the data itself.

The process consists of creating and running 2 programs:

apex2000_prog 

apex_prog

The whole process is controlled by the script:

runscript.sh

Here are the contents of "runscript.sh":

ifort -o apex2000_prog apex2000.f         ! need to convert to your compiler

make                                      ! Makefile to create apex_prog 
                                          ! (again, need to convert to your compiler)

./apex2000_prog  < input_date > outfile   ! run apex_2000_prog using input_date as input

./apex_prog < input_date >> outfile       ! run apex_prog using input_date as input



The first program, apex2000_prog creates an intermediate file "Apex_grid_data" which is roughly 31MB
and produces a global grid based on the Apex coordinate system

The second program uses "Apex_grid_data" and creates the large file that is needed by GIP to
define all of the flux-tubes and the transforms between the geographic grid and these tubes
(in both directions).  This program also uses the data file "tiegcm_defined_apex_heights"
which is never edited.

The final output file is:

GIP_apex_coords_etc.2000.0.format

This is an ascii double-precision file of size 1.4GB

Note:  once you have calculated the final file "GIP_apex_coords_etc.2000.0.format" - this is the only output
file needed - ie, the intermediate file "Apex_grid_data" can be discarded..... it is not needed to run GIP.

Also: Within all of this there are several references to "2000" - this just reflects the fact that the
year 2000 was used extensively when this was being built.  The code will work for all dates from 1900 - 
2000 - all you need to do to is to edit the file "input_date" - and then rename the output file accordingly.
The whole system does not work for 2010 at this point - the underlying IGRF datasets need to be updated
for it work with 2010 and beyond.... it might work for 2005 - but I'm not sure....in any case, it needs to
be updated to cope with future dates.

So, to create files for 1975:

edit input_date
run runscript.sh
rename GIP_apex_coords_etc.2000.0.format -> GIP_apex_coords_etc.1975.0.format
Then define this 1975 file in the GIP input as the correct static file.

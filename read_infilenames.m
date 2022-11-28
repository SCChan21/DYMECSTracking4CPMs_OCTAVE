function [nfiles, filenames, start_date, timestep, ntimes_per_file]=read_infilenames(inpath, infiles)
%% reads the content of the infiles, ie the start date and list of filenames containing the data to use for tracking
%% Email: j.a.crook@leeds.ac.uk
%% Last update 14 June 2016.

% Changed 20190522
% textscan does not work the same in Octave, so it must be changed to work in a more ugly way
% file is now formatted as yyyy,mm,dd,hh,mm,timestep,ntimes_per_file,file1 file2 file3... (SPACE SEPERATED, NO COMMA BETWEEN FILENAMES NOR AT THE END!!!!)

fidn=fopen([inpath infiles],'r');
uglycommaseperatedstuff = textscan(fidn,'%n %n %n %n %n %n %n %s', 'delimiter', ',');
start_date.year   = uglycommaseperatedstuff{1};
start_date.month  = uglycommaseperatedstuff{2};
start_date.day    = uglycommaseperatedstuff{3};
start_date.hour   = uglycommaseperatedstuff{4};
start_date.minute = uglycommaseperatedstuff{5};
timestep          = uglycommaseperatedstuff{6};
ntimes_per_file   = uglycommaseperatedstuff{7};
filenames         = strsplit(uglycommaseperatedstuff{8}{1});
[n, nfiles]       = size(filenames);
fclose(fidn);
clear fidn;
return
%end

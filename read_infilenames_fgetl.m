function [nfiles, filenames, start_date, timestep, ntimes_per_file]=read_infilenames_fgetl(inpath, infiles)
%% reads the content of the infiles, ie the start date and list of filenames containing the data to use for tracking
%% Email: j.a.crook@leeds.ac.uk
%% Last update 14 June 2016.

% Changed 20190522
% textscan does not work the same in Octave, so it must be changed to work in a more ugly way
% file is now formatted as yyyy,mm,dd,hh,mm,timestep,ntimes_per_file,file1 file2 file3... (SPACE SEPERATED, NO COMMA BETWEEN FILENAMES NOR AT THE END!!!!)

fidn=fopen([inpath infiles],'r');
start_date.year   = str2num(fgetl(fidn));
start_date.month  = str2num(fgetl(fidn));
start_date.day    = str2num(fgetl(fidn));
start_date.hour   = str2num(fgetl(fidn));
start_date.minute = str2num(fgetl(fidn));
timestep          = str2num(fgetl(fidn));
ntimes_per_file   = str2num(fgetl(fidn));
filenames = {}
while ~feof(fidn)
    tline = fgetl(fidn);
    filenames(1, end + 1) = tline;
    disp(tline)
end
[n, nfiles] = size(filenames);
fclose(fidn);
clear fidn;
return
%end

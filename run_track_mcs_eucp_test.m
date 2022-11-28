#!/usr/bin/octave -qf
%% Code to set up inputs to track_mcs and run it

inpath=['./'];
testpath='eucp_test'
trackpath='.'

outpath=[testpath '/' trackpath '/']
config_file=[testpath '/'  trackpath '/rpole12_precip.cfg']
infiles=[testpath '/rpole_12_MultiLine.txt']
rifn=@read_infilenames_fgetl

ok=track_mcs(inpath, config_file, infiles, outpath, rifn);

%end

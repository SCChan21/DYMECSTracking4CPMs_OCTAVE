%% Code to set up inputs to track_mcs and run it

inpath=['./default_test/'];
outpath=['./default_test/'];
config_file='vera_precip.cfg';
infiles='precip_filenames_30062014-31072014.txt';
rifn=@read_infilenames

ok=track_mcs(inpath, config_file, infiles, outpath);
%end

function config=read_config_file(inpath, config_file)
%% reads the content of the config file required for tracking
%% Email: j.a.crook@leeds.ac.uk
%% Last update 14 June 2016.

% Changed 20190522
% textscan does not work the same in Octave, so it must be changed to work in a more ugly way
% file is now formatted as minarea,blockradius,squaresize,halo,olapthresh,threshold,less_than,varname (COMMA SEPERATED, DO NOT PUT COMMA AFTER varname!!!!)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read config file to get minarea, blockradius, squaresize, halo, olapthresh, threshold, less_than, varname
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fflush(stdout);
fflush(stderr);
%
configfile = [inpath config_file];
fidn = fopen(configfile,'r');
uglycommaseperatedstuff = textscan(fidn,'%n %n %n %n %n %n %n %s', 'delimiter', ',');
minarea     = uglycommaseperatedstuff{1};
blockradius = uglycommaseperatedstuff{2};
squaresize  = uglycommaseperatedstuff{3};
halo        = uglycommaseperatedstuff{4};
olapthresh  = uglycommaseperatedstuff{5};
threshold   = uglycommaseperatedstuff{6};
less_than   = uglycommaseperatedstuff{7};
varname     = uglycommaseperatedstuff{8}{1};
%
config.minarea     = minarea;
config.blockradius = blockradius;
config.squaresize  = squaresize;
config.halo        = halo;
config.olapthresh  = olapthresh;
config.threshold   = threshold;
config.less_than   = less_than;
config.varname     = varname;
fclose(fidn);
clear fidn;
%
disp(['minarea     = ' num2str(minarea)]);
disp(['blockradius = ' num2str(blockradius)]);
disp(['squaresize  = ' num2str(squaresize)]);
disp(['halo        = ' num2str(halo)]);
disp(['olapthresh  = ' num2str(olapthresh)]);
disp(['threshold   = ' num2str(threshold)]);
disp(['less_than   = ' num2str(less_than)]);
disp(['varname     = ' varname]);
%
return
%end

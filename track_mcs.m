function ok=track_mcs(inpath, config_file, infiles, outpath, rifn, GaussianSmooth = 0, CThreeSixty = 0)
%% Tracking code written by Thorwald Stein and Julia Crook
%% Email: t.h.m.stein@reading.ac.uk j.a.crook@leeds.ac.uk
%% Acquired for use with SEVIRI brightness temperatures and modified to work with Cascade data
%% Simplified for external users.
%% Last update 14 June 2016.

pkg ("load", "auto");

pkg load netcdf
import_netcdf
pkg load statistics
[exitcode, whichout] = system('which gnuplot-wx'); % For some bizzare reason, it kept that funny \n at the end...
disp(['gnuplot found at ' whichout(1:end-1)]);
gnuplot_binary(whichout(1:end-1)); 
%gnuplot_binary(whichout(1:end-1) ' -geometry 600x400+600+400'); 
setenv("GNUTERM","x11")

% Forced flush stdout/stderr: Similar to Python, stdout and stderr are sometimes held in buffer
fflush(stdout);
fflush(stderr);

ok=1;

config=read_config_file(inpath, config_file);
[nfiles, ncfilenames, start_date, timestep, ntimes_per_file]=rifn(inpath, infiles);
%[nfiles, ncfilenames, start_date, timestep, ntimes_per_file]=read_infilenames(inpath, infiles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Q is the struct array from the previous image
%% M will be struct array for current image
%% Initialised here as an empty array
%% Struct arrays have fields:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% label                : array with labels indicating individual regions
%% S(i)                 : struct array containing properties of region i
%% S(i).area            : number of pixels in region i
%% S(i).centroid        : [x,y] centre of region i
%% S(i).box             : [x(ul),y(ul),width_x,width_y] bounding rectangle
%% S(i).was             : original label (for storm)
%% S(i).life            : lifetime of storm (with consistent label)
%% S(i).u               : horizontal displacement from T-1 to T
%% S(i).v               : vertical displacement from T-1 to T
%% S(i).parent          : major cell from break-up (name children)
%% S(i).child           : minor cell from break-up (name parent)
%% S(i).wasdist         : distance of centroids of current storm and original
%% S(i).track.xpos(j)   : x-position of centroid of region i at time (in lifetime) j
%% S(i).track.ypos(j)   : y-position of centroid of region i at time (in lifetime) j
%% S(i).accreted        : previously more than 1 cell (merged)
%% S(i).cell            : list of cores C associated with storm i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% C(j).storm           : storm associated with core j
%% C(j).index           : index location of core j
%% C(j).maxrain         : value associated with core j
%% C(j).was             : original label (for cell)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIRST
%% LOOP THROUGH OR LOAD IMAGE TO IDENTIFY CELLS
%% AND COMPARE WITH PREVIOUS IMAGE(S) FOR TRACKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

firstfile=1;
newwas=[];
newcell=[];

t=1; % counter for which timestep within a file
tix=1; % index into the data if more than one timestep in a file 
n=1;
NT=ntimes_per_file;
month=start_date.month;
day=start_date.day;
hour=start_date.hour;
minute=start_date.minute;

while n <= nfiles

  warning off 
  tic

  if month<10
    mstr=['0' num2str(month)];
  else
    mstr=[num2str(month)];
  end;
  
  if day<10
    dstr=['0' num2str(day)];
  else
    dstr=[num2str(day)];
  end;
  
  if hour<10
    hstr=['0' num2str(hour)];
  else
    hstr=[num2str(hour)];
  end;
  
  if minute<10
    minstr=['0' num2str(minute)];
  else
    minstr=[num2str(minute)];
  end;

  if (t == 1)
 
    % we need to read the next file
    disp(['Current ncfile = ' ncfilenames{n}]);
    this_ncfname=ncfilenames{n};
 
    % Check file existance (skip if needed)
    if isempty(dir(this_ncfname))
        disp([this_ncfname ' is missing']);
        while t <= NT
          if CThreeSixty
            [month,day,hour,minute]=go_to_next_time_360d(month,day,hour,minute,timestep);
          else
            [month,day,hour,minute]=go_to_next_time(month,day,hour,minute,timestep);
          end; 
          t=t+1;
        end;
        n=n+1;
        t=1;
        tix=1;
        continue;
    end;

    % Read file (this reads the entire file? -- could be a source of memory problem)
    ncid = netcdf.open(this_ncfname, 'NOWRITE');
    indata = ncread(this_ncfname, config.varname); % Broken pipe without ; ?????
    netcdf.close(ncid);
 
    % squeeze: UM data has 4 dimensions - 3rd one for altitude and is not needed so remove it
    indata=squeeze(indata);
    ncfname_split=strsplit(this_ncfname,'/'); % split off the pathname and filename
    lastpart=ncfname_split{end};
    outfile_bname=strsplit(lastpart,'\.nc'); % remove the nc file extension
    disp([outfile_bname{1} ' nc File loaded ' num2str(toc)]); % toc = clock timer
    this_NT = size(indata,3);
    if this_NT < NT
        % when there has been a missing timestep in a cascade file it has always been the 1st timestep
        nmissing = NT - this_NT
        while t <= nmissing
          if CThreeSixty
            [month,day,hour,minute]=go_to_next_time_360d(month,day,hour,minute,timestep);
          else
            [month,day,hour,minute]=go_to_next_time(month,day,hour,minute,timestep);
          end;
          t = t+1;
        end;
        continue;
    end 
  end
  this_image=indata(:,:,tix);
  if firstfile
      prevhour=[];
      prevmin=[];
      % NLINES and NPIXELS must be the same in all files
      NLINES=size(indata,1);
      NPIXELS=size(indata,2);
      [X,Y]=meshgrid(1:NPIXELS,1:NLINES);
      old_image=[];
      firstfile=0;
  end;

  % base the output filename on the ncfilename and the index t
  outfilename=[outfile_bname{1} '_' num2str(t)];
  disp(['outfilename = ' outfilename]);

  [M, newwas, newcell] = do_tracking(this_image, old_image, NLINES, NPIXELS, X, Y, month, day, hour, minute, prevhour, prevmin, timestep, outpath, outfilename, config, newwas, newcell, Q, GaussianSmooth);
  prevhour = hour;
  prevmin = minute;
  if CThreeSixty
    [month,day,hour,minute] = go_to_next_time_360d(month,day,hour,minute, timestep);
  else
    [month,day,hour,minute] = go_to_next_time(month,day,hour,minute, timestep);
  end;
  disp(['Writing and plotting done in ' num2str(toc)]);
  old_image=this_image;
  if max(M.label(:))>0 %% Set Q (T-1) to be M (T) for next time step
    clear Q;
    Q=M;
  else  %%  reset Q
    Q=[];
  end;
  clear M;
  % set up t and NT for next time round loop
  t=t+1;
  tix=tix+1;
  if t > NT
    % we need to read next file
    n=n+1;
    t=1;
    tix=1;
  end
end %% WHILE n <= nfiles

%clear all
return
%end

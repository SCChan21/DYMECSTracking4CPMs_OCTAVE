function [M,newwas, newcell]=do_tracking(this_image, old_image, NLINES, NPIXELS, X, Y, month, day, hour, minute, prevhour, prevmin, timestep, outpath, outfilename, config, newwas, newcell, Q, GaussianSmooth = 0)
%% Tracking code written by Thorwald Stein and updates by Julia Crook
%% Email: t.h.m.stein@reading.ac.uk j.a.crook@leeds.ac.uk
%% Acquired for use with SEVIRI brightness temperatures and modified to work with Cascade data
%% Simplified for external users.
%% Last update 14 June 2016.
%% 25/10/16: added meanval, minval and maxval to storm structure


misval=-999;    %% missing value
halosq=config.halo^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Q is the struct array from the previous image
%% M will be struct array for current image
%% Initialised here as an empty array
%% Struct arrays have fields:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% label                : array with labels indicating individual regions
%% S(i)                 : struct array containing properties of region i
%% S(i).area            : number of pixels in region i
%% S(i).meanval         : mean value within this storm
%% S(i).minval          : minimum value within this storm
%% S(i).maxval          : maximum value within this storm
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

matlab_file=[outpath outfilename '.mat'];
li = dir(matlab_file);
if isempty(li)
  do_track = 1;
  try
    M = define_storms(this_image,config.minarea,config.threshold,config.blockradius,misval,config.less_than,1,GaussianSmooth);
  catch
    disp(['Ah Oh... something messed up in define_storms']);
  	continue
  end;
  tic
else
  load(matlab_file);
  disp([outfilename ' matlab File loaded ' num2str(toc)]);
  do_track = 0;
  if isempty(M.S)
      maxwas=1;
  else
    tempwas=[];
    for ms=1:length(M.S)
      tempwas(end+1)=M.S(ms).was;
    end
    maxwas=max(tempwas)+1;
  end
  if isempty(newwas)
    newwas = maxwas;
  else
    if newwas < maxwas
      newwas=maxwas;
    end
  end
  if isempty(M.C)
    maxcell=1;
  else
    tempwas=[];
    for ms=1:length(M.C)
      tempwas(end+1)=M.C(ms).was;
    end
    maxcell=max(tempwas)+1;
  end
  if isempty(newcell)
    newcell=maxcell;
  else
    if newcell < maxcell
      newcell=maxcell;
    end
  end
  tic
end;

close all

if do_track

if isempty(prevhour)
  mistime=1;
else
  if prevhour==hour
    mistime=(minute-prevmin)/timestep;
  elseif prevhour<hour
    mistime=((hour-prevhour)*60+minute-prevmin)/timestep;
  else
    mistime=((hour+24-prevhour)*60+minute-prevmin)/timestep;
  end
end

disp([hstr minstr ' - ' num2str(mistime)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIRST TIME ONLY
%% (A) FILL "WAS" AND "LIFE" VALUES 
%% IN STORM ARRAY TO SET INITIAL LABELS
%% AND STORM LIFE TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['Start tracking at ' num2str(toc)]);
tic
pack
disp(['Finished packing at ' num2str(toc)]);
tic %% Check time spent on relabelling and tracking

if isempty(Q) 
waslabels=[];
tempmax=max(M.label(:));
ii=1;
cellnum=1;
for labelnum=1:tempmax
  newmax=max(M.label(:));
  if labelnum>newmax;
    break;
  end;
  centrind=find(M.label==labelnum);
  if isempty(centrind)
    continue;
  end;
  if length(centrind)<config.minarea
    M.label(centrind)=0;
    M.locmax(centrind)=0;
  end;
  M.label(centrind)=ii;
  M.S(ii).area=length(centrind);
  M.S(ii).meanval=mean(this_image(centrind));
  M.S(ii).minval=min(this_image(centrind));
  M.S(ii).maxval=max(this_image(centrind));
  M.S(ii).centroid = [mean(X(centrind)) mean(Y(centrind))];
  M.S(ii).box = [min(X(centrind)) min(Y(centrind)) max(X(centrind))-min(X(centrind)) max(Y(centrind))-min(Y(centrind))];
  M.S(ii).was = ii;
  M.S(ii).life = 1;
  M.S(ii).track.xpos = M.S(ii).centroid(1);
  M.S(ii).track.ypos = M.S(ii).centroid(2);
  M.S(ii).u = 0;
  M.S(ii).v = 0;
  M.S(ii).parent = [misval];
  M.S(ii).child = [misval];
  M.S(ii).wasdist = [misval];
  M.S(ii).accreted = [misval];
  nhind=find(M.locmax==1 & M.label==ii);
  M.S(ii).cell=misval;
  if 1-isempty(nhind)
    for kk=1:length(nhind)
      if M.S(ii).cell(end)==misval
        M.S(ii).cell(end)=cellnum;
      else
        M.S(ii).cell(end+1)=cellnum;
      end;
      M.C(cellnum).was = cellnum;
      M.C(cellnum).storm = ii;
      M.C(cellnum).index = nhind(kk);
      M.C(cellnum).maxrain=this_image(nhind(kk));
      cellnum=cellnum+1;
    end;
  end;
  waslabels(end+1)=M.S(ii).was;
  %% FIND ELLIPSE VALUES - crystal5 function not available
  if M.S(ii).area<0
    [xa,ya]=find(X>=M.S(ii).box(1) & X<=M.S(ii).box(1)+M.S(ii).box(3) & Y>=M.S(ii).box(2) & Y<=M.S(ii).box(2)+M.S(ii).box(4));
    A=flipud(M.label(min(xa):max(xa),min(ya):max(ya)));
    A(find(A~=ii))=0;
    A(find(A>0))=1;
    [Dm,Area,alpha,theta,Dmax,Dshort,Dlong,leng,Dx,Dz] ...
      =crystal5(A, 0, 0, 1);
     M.S(ii).ellipse.area=pi*Dlong*Dshort/4;
     M.S(ii).ellipse.axialratio=alpha;
     M.S(ii).ellipse.orientation=theta;
%   pause
%   close
  else
     M.S(ii).ellipse.area=misval;
     M.S(ii).ellipse.axialratio=misval;
     M.S(ii).ellipse.orientation=misval;  
  end;
  ii=ii+1;
end;
newwas=max(waslabels)+1;
newcell=cellnum;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANY OTHER TIME:
%% FILL EACH STORM WITH INFORMATION AND DO
%% TRACKING 
%% FOR EACH STORM, FIND OVERLAP WITH STORMS AT T-1 (Q)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp(['Figure initiated'])
%toc

elseif 1-isempty(M)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Q and M are not empty, so use fft to get velocities
%% and update uvlabel in Q accordingly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  if config.less_than 
      oldbt=old_image;
      newbt=this_image;
      threshold2=config.threshold+10;
      oldbt(find(Q.label<1))=threshold2;
	  oldbt=oldbt-threshold2;
      newbt(find(M.label<1))=threshold2;
	  newbt=newbt-threshold2;
  else
      oldbt=log(old_image);
      newbt=log(this_image);
      oldbt(find(Q.label<1))=NaN;
      newbt(find(Q.label<1))=NaN;
  end;
  % if NLINES or NPIXELS is not divisible by squaresize we will not be able
  % to look at whole image - miss out pixels at edge
  max_corx=floor(2*NLINES/config.squaresize)-1;
  max_cory=floor(2*NPIXELS/config.squaresize)-1;
  start_x=1+floor((NLINES-(config.squaresize+(max_corx-1)*config.squaresize/2))/2);
  start_y=1+floor((NPIXELS-(config.squaresize+(max_cory-1)*config.squaresize/2))/2);
  maxline=config.squaresize*(max_corx+1)/4;
  maxpixel=config.squaresize*(max_cory+1)/4;
  [xmat,ymat]=meshgrid(-1*maxpixel:(config.squaresize/2):maxpixel,-1*maxline:(config.squaresize/2):maxline);
  buu=NaN*ones(size(xmat));
  bvv=buu;
  bww=buu;
  bss=buu;

  for corx=1:max_corx
    x1=(config.squaresize/2)*(corx-1)+start_x;
    x2=x1+config.squaresize-1;
    for cory=1:max_cory
      y1=(config.squaresize/2)*(cory-1)+start_y;
      y2=y1+config.squaresize-1;

      oldsquare=(oldbt(x1:x2,y1:y2));
      newsquare=(newbt(x1:x2,y1:y2));
      if config.less_than
        templen=length(find(oldsquare<=config.threshold-threshold2));
      else
        templen=length(isnan(oldsquare));
      end;
      bss(corx+1,cory+1)=templen;
      disp(newsquare)    
      sumold=nansum(oldsquare(:)); % sum(x, 'omitnan') becomes nansum(x), requires package statistics
      sumnew=nansum(newsquare(:)); % sum(x, 'omitnan') becomes nansum(x), requires package statistics
      if sumold==0 || sumnew==0 || isnan(sumold) || isnan(sumnew)
        buu(corx+1,cory+1)=NaN;
        bvv(corx+1,cory+1)=NaN;
        bww(corx+1,cory+1)=NaN;
      else
        %% 1=TUKEY WINDOW
        %% 0=NO WINDOW (RECTANGULAR)
        [dy,dx,amplitude]=fft_track(oldsquare,newsquare,1);
        buu(corx+1,cory+1)=dx/mistime;
        bvv(corx+1,cory+1)=-dy/mistime;%% indices are upside down so need minus to get real-world v-velocity
        bww(corx+1,cory+1)=amplitude;
      end
    end;
  end;
  %% ACTUAL VELOCITY
  baa=sqrt(buu.^2+bvv.^2);
  display('buu:');
  display(buu);
  display('bvv:');
  display(bvv);
  ind=find(isnan(baa)==0);
  buu(find(baa>mean(baa(ind))+2*std(baa(ind))))=NaN;
  bvv(find(baa>mean(baa(ind))+2*std(baa(ind))))=NaN;
  %% ADHOC SOLUTION FOR HCLIMcom-HCLIM38-AROME
  #{
  if length(ind) > 0
    buu(find(baa>mean(baa(ind))+2*std(baa(ind))))=NaN;
    bvv(find(baa>mean(baa(ind))+2*std(baa(ind))))=NaN;
  else
    buu(:)=NaN;
    bvv(:)=NaN;
  end; 
  #}
  %% END END END
  uumat=inpaint_nans(buu,5);
  vvmat=inpaint_nans(bvv,5);
  X_zero_centred=X-ceil(NPIXELS/2);
  Y_zero_centred=Y-ceil(NLINES/2);
  newumat=interp2(xmat,ymat,flipud(uumat),X_zero_centred, Y_zero_centred,'cubic');
  newvmat=interp2(xmat,ymat,flipud(vvmat),X_zero_centred, Y_zero_centred,'cubic');
  uv_filename=[outpath outfilename '_UV'];
  save(uv_filename,'newumat','newvmat')
  
newlabel=zeros(size(Q.label));
for jj=1:length(Q.S)
  labelind=find(Q.label==jj);
  u=nanmean(newumat(labelind));  % mean(x, 'omitnan') becomes nanmean(x), requires package statistics
  v=-nanmean(newvmat(labelind)); %% Y-DIR FLIPPED UD ; % mean(x, 'omitnan') becomes nanmean(x), requires package statistics
  % have to assume zero velocity if u or v are NaN
  if isnan([u])
      u=0;
  end
  if isnan([v])
      v=0;
  end
  if (u==0 && v==0) 
    newlabel(labelind)=jj;
  else
    [yind,xind]=find(Q.label==jj);
    for ii=1:length(labelind)
      newxind=xind(ii)+round(u);
      newyind=yind(ii)+round(v);
      if newxind>size(newlabel,2) || newyind>size(newlabel,1) || newxind<1 || newyind<1
        continue;
      elseif newlabel(newyind,newxind)>0
        olddist=(X(newyind,newxind)-Q.S(newlabel(newyind,newxind)).centroid(1))^2+(Y(newyind,newxind)-Q.S(newlabel(newyind,newxind)).centroid(2))^2;
        newdist=(X(newyind,newxind)-Q.S(jj).centroid(1))^2+(Y(newyind,newxind)-Q.S(jj).centroid(2))^2;
        if newdist<olddist
          newlabel(newyind,newxind)=jj;
        end
      else
        newlabel(newyind,newxind)=jj;
      end;
    end;
  end;
end;
Q.uvlabel=newlabel;
for jj=1:length(Q.S)   
  centrind=find(Q.uvlabel==jj);
  Q.S(jj).uvcentroid = [mean(X(centrind)) mean(Y(centrind))];
  Q.S(jj).uvarea = length(centrind);
end;
clear newlabel
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now loop through M and check for
%% overlap with advected Q storms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  tempmax=max(M.label(:));
  jj=1;
  wasnum=zeros(1,length(M.S));
  if 1-isempty(M.C)
    wascell=zeros(1,length(M.C));
  else
    wascell=[];
  end;
  cellnum=1;
  qbins=[0:length(Q.S)];
  qarea=1;
  qlife=1;
  for qq=1:length(Q.S)
    qarea(qq+1)=Q.S(qq).uvarea;
    qlife(qq+1)=Q.S(qq).life;
  end;
  for labelnum=1:tempmax
    newmax=max(M.label(:));
    if labelnum>newmax
      break;
    end;
    centrind=find(M.label==labelnum);
    if isempty(centrind)
      continue;
    end;
    if length(centrind)<config.minarea
      M.label(centrind)=0;
      M.locmax(centrind)=0;
      continue;
    end;
%    disp(num2str(jj))
    M.label(centrind)=jj;
    M.S(jj).area=length(centrind);
    M.S(jj).meanval=mean(this_image(centrind));
    M.S(jj).minval=min(this_image(centrind));
    M.S(jj).maxval=max(this_image(centrind));
    M.S(jj).centroid = [mean(X(centrind)) mean(Y(centrind))];
    M.S(jj).box = [min(X(centrind)) min(Y(centrind)) max(X(centrind))-min(X(centrind)) max(Y(centrind))-min(Y(centrind))];
    M.S(jj).was = jj;
    M.S(jj).life = 1;
    M.S(jj).track.xpos = M.S(jj).centroid(1);
    M.S(jj).track.ypos = M.S(jj).centroid(2);
    u=nanmean(newumat(centrind)); % mean(x, 'omitnan') becomes nanmean(x), requires package statistics
    % have to assume zero velocity if u or v are NaN
    if isnan([u]) 
        u=0;
    end;
    M.S(jj).u = u;
    v=nanmean(newvmat(centrind)); % mean(x, 'omitnan') becomes nanmean(x), requires package statistics
    if isnan([v])
       v=0;
    end
    M.S(jj).v = v;
    M.S(jj).parent = [misval];
    M.S(jj).child = [misval];
    M.S(jj).wasdist = [misval];
    M.S(jj).accreted = [misval];
    M.S(jj).cell=misval;
    nhind=find(M.locmax==1 & M.label==jj);
    if 1-isempty(nhind)
    for kk=1:length(nhind)
      if M.S(jj).cell(end)==misval
        M.S(jj).cell(end)=cellnum;
      else
        M.S(jj).cell(end+1)=cellnum;
      end;
      M.C(cellnum).was = cellnum;
      M.C(cellnum).storm = jj;
      M.C(cellnum).index = nhind(kk);
      M.C(cellnum).maxrain=this_image(nhind(kk));
      cellnum=cellnum+1;
    end;
    end;
  %% FIND ELLIPSE VALUES - crystal5 function not available
  if M.S(jj).area<0
    [xa,ya]=find(X>=M.S(jj).box(1) & X<=M.S(jj).box(1)+M.S(jj).box(3) & Y>=M.S(jj).box(2) & Y<=M.S(jj).box(2)+M.S(jj).box(4));
    A=flipud(M.label(min(xa):max(xa),min(ya):max(ya)));
    A(find(A~=jj))=0;
    A(find(A>0))=1;
    [Dm,Area,alpha,theta,Dmax,Dshort,Dlong,leng,Dx,Dz] ...
      =crystal5(A, 0, 0, 1);
     M.S(jj).ellipse.area=pi*Dlong*Dshort/4;
     M.S(jj).ellipse.axialratio=alpha;
     M.S(jj).ellipse.orientation=theta;
   else
     M.S(jj).ellipse.area=misval;
     M.S(jj).ellipse.axialratio=misval;
     M.S(jj).ellipse.orientation=misval;
   end
    %disp(['Storm ' num2str(jj)])
    %toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHECK OVERLAP WITH QHIST
%% IF NO OVERLAP, THEN
%% GENERATE (halo)km RADIUS AROUND CENTROID
%% CHECK FOR OVERLAP WITHIN (halo)km OF CENTROID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    qhist=hist(Q.uvlabel(find(M.label==jj)),qbins)/M.S(jj).area+hist(Q.uvlabel(find(M.label==jj)),qbins)./qarea;
    newblob=[];
    if max(qhist(2:end))<config.olapthresh
      newblob=0*X;
      blobind=find((X-M.S(jj).centroid(1)).^2 + (Y-M.S(jj).centroid(2)).^2<halosq);
      newblob(blobind)=newblob(blobind)+1;
      qhist=hist(Q.uvlabel(find(newblob==1)),qbins)/M.S(jj).area+hist(Q.uvlabel(find(newblob==1)),qbins)./qarea;
    end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IF OVERLAP, THEN 
%% - INHERIT "WAS"
%% - UPDATE "LIFE" AND "TRACK" AND "WASDIST"
%% - INHERIT "U" AND "V" (ONLY UPDATE IF SINGLE OVERLAP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if max(qhist(2:end))>=config.olapthresh
      numlaps=find(qhist(2:end)>=config.olapthresh);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IF MORE THAN ONE BEST OVERLAP
%% KEEP PROPERTIES OF NEAREST IN CENTROID       
%% IF MORE THAN ONE LARGEST, KEEP NEAREST IN CENTROID     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
      if length(numlaps)>1
        lapdist=[];
        sectlap=[];
        for kkind=1:length(numlaps)
          lapdist(kkind)=sqrt((M.S(jj).centroid(1)-Q.S(numlaps(kkind)).uvcentroid(1))^2 + (M.S(jj).centroid(2)-Q.S(numlaps(kkind)).uvcentroid(2))^2);
          sectlap(kkind)=length(find(Q.uvlabel==numlaps(kkind) & M.label==jj));
        end;
        kkmax=find(sectlap==max(sectlap));
        if length(kkmax)>1
          kkmax=kkmax(find(lapdist(kkmax)==min(lapdist(kkmax))));
        end;
        kkmax=kkmax(1); %% ALL ELSE EQUAL, JUST PICK FIRST
        M.S(jj).was = Q.S(numlaps(kkmax)).was;
        qnum(jj) = numlaps(kkmax);
        M.S(jj).life = Q.S(numlaps(kkmax)).life + mistime;
        M.S(jj).track.xpos(1:Q.S(numlaps(kkmax)).life) = Q.S(numlaps(kkmax)).track.xpos(1:Q.S(numlaps(kkmax)).life);
        M.S(jj).track.ypos(1:Q.S(numlaps(kkmax)).life) = Q.S(numlaps(kkmax)).track.ypos(1:Q.S(numlaps(kkmax)).life);
        M.S(jj).track.xpos(M.S(jj).life) = M.S(jj).centroid(1);
        M.S(jj).track.ypos(M.S(jj).life) = M.S(jj).centroid(2);
        lapdist=sqrt((M.S(jj).centroid(1)-Q.S(numlaps(kkmax)).centroid(1))^2 + (M.S(jj).centroid(2)-Q.S(numlaps(kkmax)).centroid(2))^2);
        M.S(jj).wasdist = length(find(Q.uvlabel==numlaps(kkmax) & M.label==jj));
        qind=numlaps(kkmax);
        alllaps=find(qhist(2:end)>=config.olapthresh);
        for kkind=1:length(alllaps)
          if alllaps(kkind)==numlaps(kkmax)
            continue;
          end;
          if M.S(jj).accreted(end)==misval
            M.S(jj).accreted(end)=Q.S(alllaps(kkind)).was;
          else
            M.S(jj).accreted(end+1)=Q.S(alllaps(kkind)).was;
          end;
        end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SINGLE OVERLAP     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

      else
        lapdist=sqrt((M.S(jj).centroid(1)-Q.S(numlaps).centroid(1))^2 + (M.S(jj).centroid(2)-Q.S(numlaps).centroid(2))^2);
        M.S(jj).was = Q.S(numlaps).was;
        qnum(jj) = numlaps;
        M.S(jj).life = Q.S(numlaps).life + mistime;
        M.S(jj).track.xpos(1:Q.S(numlaps).life) = Q.S(numlaps).track.xpos(1:Q.S(numlaps).life);
        M.S(jj).track.ypos(1:Q.S(numlaps).life) = Q.S(numlaps).track.ypos(1:Q.S(numlaps).life);
        M.S(jj).track.xpos(M.S(jj).life) = M.S(jj).centroid(1);
        M.S(jj).track.ypos(M.S(jj).life) = M.S(jj).centroid(2);
        M.S(jj).wasdist = length(find(Q.uvlabel==numlaps & M.label==jj));
        qind=numlaps;
      end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ATTACH CELLS TO NEWLY LABELLED STORMS      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

      if M.S(jj).cell(end)>misval
        if qind==0
          for kk=1:length(M.S(jj).cell)
%%            M.S(jj).cell(kk) = newcell;
            M.C(M.S(jj).cell(kk)).storm = jj;
            M.C(M.S(jj).cell(kk)).was = newcell;
            wascell(M.S(jj).cell(kk)) = newcell;
            newcell=newcell+1;
%            if newcell>9999
%              newcell=1;
%            end;
            %
            % Below was commented out
            %
%            disp('plot 1 is here!');
%            texstring=[num2str(M.C(M.S(jj).cell(kk)).maxrain)];
%            plot(X(M.C(M.S(jj).cell(kk)).index),Y(M.C(M.S(jj).cell(kk)).index),'r.');
%            text(X(M.C(M.S(jj).cell(kk)).index),Y(M.C(M.S(jj).cell(kk)).index),texstring,'fontsize',10,'color','blue');
            %
            % end of commented out block
            %
          end;
        else %% INHERIT CELL LABELS FROM PREVIOUS STORM......
          if M.S(jj).u == misval || M.S(jj).v ==misval %% CANNOT TRACK
%%            disp('AAA')
            for kk=1:length(M.S(jj).cell)
%%              M.S(jj).cell(kk) = newcell;
              M.C(M.S(jj).cell(kk)).storm = jj;
              M.C(M.S(jj).cell(kk)).was = newcell;
              wascell(M.S(jj).cell(kk)) = newcell;
              newcell=newcell+1;
     %         if newcell>9999
      %          newcell=1;
      %        end;
            end;
%%            disp('AAA111')
          else %% TRACK BACK AND CHECK WITH CELLS IN Q
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if M.S(jj).accreted(end)~=misval %% TRICKY!!!
%            disp('BBB')
            kkdist=[];
            kknum=[];
            %%  numlaps contains indices of Q storms merged into M storm
              for nn=1:length(numlaps)
                if Q.S(numlaps(nn)).cell(end)~=misval
                 for qq=1:length(Q.S(numlaps(nn)).cell)
                  for kk=1:length(M.S(jj).cell)
                    backtrack(kk,1)=X(M.C(M.S(jj).cell(kk)).index)-M.S(jj).u-X(Q.C(Q.S(numlaps(nn)).cell(qq)).index);
                    backtrack(kk,2)=Y(M.C(M.S(jj).cell(kk)).index)-M.S(jj).v-Y(Q.C(Q.S(numlaps(nn)).cell(qq)).index);
                  end;
                  backdist=sqrt(backtrack(:,1).^2+backtrack(:,2).^2);
                  kknum(nn,qq)=min(find(backdist==min(backdist(:))));
                  kkdist(nn,qq)=backdist(kknum(nn,qq));
                  clear backtrack backdist
                 end;
                end;
              end;
%            disp('BBB111')
%              kkdist
%              kknum
%              pause
              if 1-isempty(kknum)
              for kk=1:length(M.S(jj).cell)
                [nnind,qqind]=find(kknum==kk & kkdist<config.halo);
                if 1-isempty(qqind)
                  if length(qqind)>1
                    qmax=0;
                    for qq=1:length(qqind)
                      if Q.C(Q.S(numlaps(nnind(qq))).cell(qqind(qq))).maxrain>qmax;
                        qmax=Q.C(Q.S(numlaps(nnind(qq))).cell(qqind(qq))).maxrain;
                        M.C(M.S(jj).cell(kk)).was=Q.C(Q.S(numlaps(nnind(qq))).cell(qqind(qq))).was;
                        wascell(M.S(jj).cell(kk)) = Q.C(Q.S(numlaps(nnind(qq))).cell(qqind(qq))).was;
                      end;
                    end;
                  else
                    M.C(M.S(jj).cell(kk)).was=Q.C(Q.S(numlaps(nnind)).cell(qqind)).was;
                    wascell(M.S(jj).cell(kk))=Q.C(Q.S(numlaps(nnind)).cell(qqind)).was;
                  end
                  clear nnind qqind            
                else %% NEW CELL
                  for kk=1:length(M.S(jj).cell)
%                    M.S(jj).cell(kk) = newcell;
                    M.C(M.S(jj).cell(kk)).storm = jj;
                    M.C(M.S(jj).cell(kk)).was = newcell;
                    wascell(M.S(jj).cell(kk)) = newcell;
                    newcell=newcell+1;
       %             if newcell>9999
        %              newcell=1;
        %            end;
                  end;
                end; 
              end; %% LOOP THROUGH M CELLS
              else %% NEW CELL (if 1-isempty(kknum))
                  for kk=1:length(M.S(jj).cell)
%                    M.S(jj).cell(kk) = newcell;
                    M.C(M.S(jj).cell(kk)).storm = jj;
                    M.C(M.S(jj).cell(kk)).was = newcell;
                    wascell(M.S(jj).cell(kk)) = newcell;
                    newcell=newcell+1;
        %            if newcell>9999
        %              newcell=1;
        %            end;
                  end;
              end;
%              disp('BBB222')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
            elseif Q.S(qind).cell(end)==misval %% CANNOT TRACK
%              disp('CCC')
              for kk=1:length(M.S(jj).cell)
%                M.S(jj).cell(kk) = newcell;
                M.C(M.S(jj).cell(kk)).storm = jj;
                M.C(M.S(jj).cell(kk)).was = newcell;
                wascell(M.S(jj).cell(kk)) = newcell;
                newcell=newcell+1;
        %        if newcell>9999
        %          newcell=1;
        %        end;
              end;
%              disp('CCC111')
            else %% SINGLE OVERLAP
%              disp('DDD')
              kknum=[];
              kkdist=[];
              for qq=1:length(Q.S(qind).cell)
                for kk=1:length(M.S(jj).cell)
                  backtrack(kk,1)=X(M.C(M.S(jj).cell(kk)).index)-M.S(jj).u-X(Q.C(Q.S(qind).cell(qq)).index);
                  backtrack(kk,2)=Y(M.C(M.S(jj).cell(kk)).index)-M.S(jj).v-Y(Q.C(Q.S(qind).cell(qq)).index);
                end;
                backdist=sqrt(backtrack(:,1).^2+backtrack(:,2).^2);
                kknum(qq)=min(find(backdist==min(backdist(:))));
                kkdist(qq)=backdist(kknum(qq));
                clear backtrack backdist
              end;
%              disp('DDD111')
%              kknum
%              kkdist
%              pause
              for kk=1:length(M.S(jj).cell)
                qqind=find(kknum==kk & kkdist<config.halo);
                if 1-isempty(qqind)
                  if length(qqind)>1
                    qmax=0;
                    for qq=1:length(qqind)
                      if Q.C(Q.S(qind).cell(qqind(qq))).maxrain>qmax;
                        qmax=Q.C(Q.S(qind).cell(qqind(qq))).maxrain;
                        M.C(M.S(jj).cell(kk)).was=Q.C(Q.S(qind).cell(qqind(qq))).was;
                        wascell(M.S(jj).cell(kk))=Q.C(Q.S(qind).cell(qqind(qq))).was;
                      end;
                    end;
                  else
                    M.C(M.S(jj).cell(kk)).was=Q.C(Q.S(qind).cell(qqind)).was;
                    wascell(M.S(jj).cell(kk))=Q.C(Q.S(qind).cell(qqind)).was;
                  end
                  clear qqind
                else %% NEW CELL
                  for kk=1:length(M.S(jj).cell)
%                    M.S(jj).cell(kk) = newcell;
                    M.C(M.S(jj).cell(kk)).storm = jj;
                    wascell(M.S(jj).cell(kk)) = newcell;
                    M.C(M.S(jj).cell(kk)).was = newcell;
                    newcell=newcell+1;
        %            if newcell>9999
        %              newcell=1;
        %            end;
                  end;
                end;            
              end; %% LOOP THROUGH M CELLS
%              disp('DDD222')
            end;   %% IF/ELSE MERGED OR SINGLE OVERLAP
          end;     %% IF/ELSE NEW OR OLD STORM
        end;       %% IF THERE ARE M CELLS
      end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IF NO OVERLAP, THEN (NEW STORM)
%% - "WAS" SET TO CURRENT MAX LABEL + 1
%% - UPDATE "LIFE" AND "TRACK" AND "WASDIST" FOR A NEW STORM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
      M.S(jj).was = newwas;
      qnum(jj) = 0;
      M.S(jj).life = 1;
      M.S(jj).track.xpos = M.S(jj).centroid(1);
      M.S(jj).track.ypos = M.S(jj).centroid(2);
      newwas=newwas+1;
  %    if newwas>9999
  %      newwas=1;
   %   end;
      if M.S(jj).cell(end)>misval
          for kk=1:length(M.S(jj).cell)
%            M.S(jj).cell(kk) = newcell;
            M.C(M.S(jj).cell(kk)).storm = jj;
            M.C(M.S(jj).cell(kk)).was = newcell;
            wascell(M.S(jj).cell(kk))=newcell;
            newcell=newcell+1;
%            if newcell>9999
%              newcell=1;
%            end;
            %
            % below was commented out
            %
%            disp('plot 2 is here!');
%            texstring=[num2str(M.C(M.S(jj).cell(kk)).maxrain)];
%            plot(X(M.C(M.S(jj).cell(kk)).index),Y(M.C(M.S(jj).cell(kk)).index),'r.');
%            text(X(M.C(M.S(jj).cell(kk)).index),Y(M.C(M.S(jj).cell(kk)).index),texstring,'fontsize',10,'color','blue');
            %
            % end of commented out block
            %
          end;
      end;
    end;
    wasnum(jj)=M.S(jj).was;
    jj=jj+1;
  end;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QUICK SANITY CHECK
%% ACCRETED SHOULD NEVER BE A VALUE
%% SIMILAR TO EXISTING STORM ID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for jj=1:length(M.S)
    if M.S(jj).accreted(end)==misval
      continue;
    else
      for acnum=1:length(M.S(jj).accreted)
        acind=find(wasnum-M.S(jj).accreted(acnum)==0);
        if 1-isempty(acind)
          M.S(jj).accreted(acnum)=misval;
        end
      end
      acnew=find(M.S(jj).accreted>misval);
      if 1-isempty(acnew)
        M.S(jj).accreted=M.S(jj).accreted(acnew);
      else
        M.S(jj).accreted=misval;
      end
    end
  end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TRACKING MERGING BREAKING
%% MULTIPLE STORMS AT T (M) MAY HAVE SAME LABEL "WAS"
%% FIND STORM WITH CENTROID AT T NEAREST TO CENTROID AT T-1
%% THIS IS THE "PARENT" STORM, 
%% "PARENT" VECTOR WITH INDICES OF NEW LABELS FOR "CHILD" STORMS
%% STORMS WITH SAME WAS BUT FURTHER FROM CENTROID ARE "CHILD", VALUE "PARENT"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %disp(['Start checking degenerate labels'])
  %toc
  
  for jj=1:length(M.S)
    if M.S(jj).wasdist==misval
      continue;
    end;
    wasind=find(wasnum==wasnum(jj));
    wassep=zeros(1,length(wasind));
    for kkind=1:length(wasind)
      wassep(kkind)=M.S(wasind(kkind)).wasdist;
    end;
    kkmax=find(wassep==max(wassep));
    clear wassep
    kkmax=min(kkmax);
    children=[];
    for kkind=1:length(wasind)
      if kkind~=kkmax
        M.S(wasind(kkind)).child = M.S(wasind(kkind)).was;
        M.S(wasind(kkind)).was = newwas;
        newwas=newwas+1;
  %      if newwas>9999
  %        newwas=1;
   %     end;
        wasnum(wasind(kkind)) = M.S(wasind(kkind)).was;
        children(end+1) = M.S(wasind(kkind)).was;
        M.S(wasind(kkind)).wasdist = [misval];
        xx=M.S(wasind(kkind)).track.xpos(:);
        yy=M.S(wasind(kkind)).track.ypos(:);
        if M.S(wasind(kkind)).cell(end)>misval
          for tt=1:length(M.S(wasind(kkind)).cell)
%%            if M.C(M.S(wasind(kkind)).cell(tt)).storm<0
              M.C(M.S(wasind(kkind)).cell(tt)).storm = wasind(kkind);
%%            end;
          end;
        end;
      else
        kkthis = kkind;
      end;
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% UPDATE PARENT STORM WITH CHILDREN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for kkind=kkthis
        if 1-isempty(children)
          M.S(wasind(kkind)).parent = children;
        end;
        xx=M.S(wasind(kkind)).track.xpos(:);
        yy=M.S(wasind(kkind)).track.ypos(:);
        if M.S(wasind(kkind)).cell(end)>misval
          for tt=1:length(M.S(wasind(kkind)).cell)
%%            if M.C(M.S(wasind(kkind)).cell(tt)).storm<0
              M.C(M.S(wasind(kkind)).cell(tt)).storm = wasind(kkind);
%%            end;
          end;
        end;
    end;      
    
  end;

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHECK CELL LABEL REDUNDANCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  if 1-isempty(wascell)
    for icell=1:length(wascell)
      cellind=find(wascell==wascell(icell));
      if length(cellind)==1
        continue;
      elseif M.S(M.C(icell).storm).parent(end)>misval
        continue;
      else %% NEW LABEL!!!
        M.C(icell).was = newcell;
        wascell(icell) = newcell;
        newcell = newcell+1;
   %     if newcell>9999
   %       newcell=1;
   %     end
      end
    end
  end
  
  clear wascell

else
  newwas=1;
  newcell=1;
end;

disp(['Tracking done in ' num2str(toc)])  
tic
pack
disp(['Finished packing at ' num2str(toc)])
tic

if isempty(li)
  disp('Saving')
  save(matlab_file,'M');    
end;

%pause

%if length(M.S)>50
%  do_ascii=2;
%else
  do_ascii=1;
%end

if do_ascii
fileascii=[outpath outfilename '.tmp'];
fid=fopen(fileascii,'w');
fprintf(fid,'missing_value=%-.0f\n',misval);
fprintf(fid,'time=%-.0f:%-.0f\n',hour,minute);
for jj=1:length(M.S)

  fprintf(fid,'storm %-.0f area=%-.0f centroid=%-.2f,%-.2f box=%-.1f,%-.1f,%-.0f,%-.0f life=%-.0f u=%-.2f v=%-.2f ',[M.S(jj).was M.S(jj).area M.S(jj).centroid(1) M.S(jj).centroid(2) M.S(jj).box(1) M.S(jj).box(2) M.S(jj).box(3) M.S(jj).box(4) M.S(jj).life-1 M.S(jj).u M.S(jj).v]);
%  if jj==1
%    fid=fopen(fileascii,'a');
%  end
  meanval=M.S(jj).meanval;
  minval=M.S(jj).minval;
  maxval=M.S(jj).maxval;
  fprintf(fid,'mean=%f min=%f max=%f accreted=',[meanval minval maxval]);
  if length(M.S(jj).accreted)>1
    for acind=1:length(M.S(jj).accreted)-1
      fprintf(fid,'%-.0f,',M.S(jj).accreted(acind));
    end;
    fprintf(fid,'%-.0f parent=',M.S(jj).accreted(end));
  else
    fprintf(fid,'%-.0f parent=',M.S(jj).accreted(end));
  end;
  fprintf(fid,'%-.0f child=',M.S(jj).child);
  if length(M.S(jj).parent)>1
    for acind=1:length(M.S(jj).parent)-1
      fprintf(fid,'%-.0f,',M.S(jj).parent(acind));
    end;
    fprintf(fid,'%-.0f cell=',M.S(jj).parent(end));
  else
    fprintf(fid,'%-.0f cell=',M.S(jj).parent(end));
  end;
  if M.S(jj).cell(end)>misval
    if length(M.S(jj).cell)>1
      for acind=1:length(M.S(jj).cell)-1
        fprintf(fid,'%-.0f,',M.C(M.S(jj).cell(acind)).was);
      end;
      fprintf(fid,'%-.0f\n',M.C(M.S(jj).cell(end)).was);
    else
      fprintf(fid,'%-.0f\n',M.C(M.S(jj).cell(end)).was);
    end;
    for ii=1:length(M.S(jj).cell)
      fprintf(fid,'cell %-.0f stormid=%-.0f centroid=%-.2f,%-.2f maxr=%f\n',[M.C(M.S(jj).cell(ii)).was M.S(M.C(M.S(jj).cell(ii)).storm).was X(M.C(M.S(jj).cell(ii)).index) Y(M.C(M.S(jj).cell(ii)).index) M.C(M.S(jj).cell(ii)).maxrain]);
    end;
  else
    fprintf(fid,'%-.0f\n',M.S(jj).cell(end));
  end;
end;
fclose(fid);
clear fid
unix(['mv ' fileascii ' ' fileascii(1:end-3) 'txt']);
end;

end; %% IF do_track
clear wasnum
disp(['Writing and plotting done in ' num2str(toc)]);

return
%end

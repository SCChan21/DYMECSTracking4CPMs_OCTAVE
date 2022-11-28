function M=define_storms(image_data,minarea,threshold,blockradius,misval,less_than,method,GaussianSmooth=0,GaussianNoise=0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tracking code written by Thorwald Stein
%% Email: t.h.m.stein@reading.ac.uk
%% and modified by Julia Crook to handle tests for > threshold or < threshold
%% 
%% Simplified for external users.
%% Last update 14 May 2015.
%% 25/10/16: added meanval, minval and maxval to storm structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THIS CODE MAKES A BINARY IMAGE OUT OF GIVEN image_data
%% USING threshold. 
%% THE BINARY IMAGE IS THEN USED AS INPUT
%% FOR LABELLING STORMS.
%% ALGORITHM ON PAGE 37-39 HARALICK AND SHAPIRO.
%% RESULT IS A STRUCT "M" CONTAINING STORMS "S" 
%% AND CELLS (LOCAL MAX/MIN) "C"
%% AND A MATRIX "labelbt" LABELLING EACH REGION S
%% AND A MATRIX "localmax" LABELLING EACH LOCAL MIN/MAX C
%% TRY pcolor(M.labelbt) TO SEE REGIONS
%% THEN LOOP THROUGH M.C TO ADD LOCATIONS OF LOCAL MAX/MIN
%% LOOP THROUGH M.S TO REPLACE misval IN STORM PROPERTIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% image_data   :: 2D vector to test for threshold
%% minarea      :: minimum size (number of pixels) of feature
%% threshold    :: threshold to distinguish features from background
%% blockradius  :: to define square region around a value to check if it's a local maximum/minimum
%% misval       :: value assigned to missing data
%% less_than    :: if 0 look for image data>threshold, else look for image_data < threshold
%% method       :: if 0 the don't include diagonals in neighbours, 1 do include diagonals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off
tic

pkg load image

image_data0 = image_data;

%%%% New: Add Gaussian Noise (this deals with some observations where data level become discrete and this messes up cell alog
%%%% Variance is set to 0.0001 (i.e. STD = 0.01)
if GaussianNoise == 1
  disp(['Applying Gaussian noise with SDEV = 0.01']);
  image_data = imnoise(image_data, 'gaussian', 0.0, 0.01*0.01);
elseif GaussianNoise == 2
  disp(['Applying Speckle (Multiplicative) Gaussian noise with SDEV = 0.025']);
  image_data = imnoise(image_data,  'speckle', 0.025*0.025);
else
  % do nothing
end;

%%%% New: Add Gaussian Blur (this deals with inputs that are higher-res than normal
if GaussianSmooth
  kernel = fspecial('Gaussian', [5 5], 0.5);
  disp(['Applying ' num2str(rows(kernel)) 'x' num2str(columns(kernel)) ' Gaussian blur kernel, weights between = ' num2str(min(kernel(:))) ' ' num2str(max(kernel(:)))]);
  image_data = imfilter(image_data, kernel, "same");
end;

%%%% Label image pixels of interest
th_data = image_data; %% Initialize th_data array
th_data(find(image_data >  threshold))=1;
th_data(find(image_data <= threshold))=0;

if less_than > 0
  disp('less than flag (reverse threshold) activated');
  th_data=1-th_data;
end;

disp(['Pre-adjusted largest value in image   = ' num2str(max(image_data0(:)))]);
disp(['Adjusted largest value in image       = ' num2str(max(image_data(:)))]);
disp(['Pre-adjusted smallest value in image  = ' num2str(min(image_data0(:)))]);
disp(['Adjusted smallest value in image      = ' num2str(min(image_data(:)))]);
disp(['Number of valid elements              = ' num2str(sum(th_data(:)))]);

labels=0*th_data;
localmax=labels;

NLINES=size(th_data,1);
NPIXELS=size(th_data,2);
newlabel=1;

disp(['Pre loop stuff done in ' num2str(toc)]);
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIRST LOOP 
%% TOP-DOWN (LEFT-RIGHT) SCAN THROUGH BINARY ARRAY
%% FIRST SET OF LABELS (labelbt) AND EQUIVALENCES ASSIGNED (eqtable)
%%
%% LOCALMAX ASSIGNED IF LOCAL MAXIMUM WITHIN BLOCKRADIUS 
%% LOCALMAX USED TO INDICATE INDIVIDUAL CELLS WITHIN A STORM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for L=1:NLINES
%  disp(['Initialize local equivalence table for line ' num2str(L)])
  eqtable=[];
%  disp('Process the line')
  for P=1:NPIXELS
    if th_data(L,P)==1
      % is this a min or max pixel?
      NHD=neighbourhood(L,P,NPIXELS,NLINES,blockradius);
      if less_than > 0
        if image_data(L,P)==min(image_data(NHD))
          localmax(L,P)=1;
        end;
      else
        if image_data(L,P)==max(image_data(NHD))
          localmax(L,P)=1;
        end;
      end;
          
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      A=neighbours_new(L,P,NPIXELS,NLINES,method);
      A=A(find(labels(A)>0));
      if isempty(A)
        NL=newlabel;
        newlabel=newlabel+1;
      else
        NL=min(labels(A));
        for X=1:length(A)
          if labels(A(X))~=NL
            if isempty(eqtable)
              eqtable=NL;
              eqtable(1,2)=labels(A(X));
            else
              [temp_A,temp_B]=find(eqtable==labels(A(X)));
              if isempty(temp_A)
                [temp_M,temp_Z]=find(eqtable==NL);
                if isempty(temp_M)
                  eqtable(size(eqtable,1)+1,1)=labels(A(X));
                  eqtable(size(eqtable,1),2)=NL;
                else
                  temp_O=find(eqtable(temp_M,:)>0);
                  eqtable(temp_M,max(temp_O)+1)=labels(A(X));
                end;
              else
                [temp_M,temp_Z]=find(eqtable==NL);
                if isempty(temp_M)
                  temp_O=find(eqtable(temp_A,:)>0);
                  eqtable(temp_A,max(temp_O)+1)=NL;
                end;
              end;
            end;
          end;
        end;
      end;
      labels(L,P)=NL;
    else
      localmax(L,P)=0;
    end;
  end;

%  disp(['Find equivalence classes detected on this line'])
  
  if isempty(eqtable)
    continue;
  else
    eqtable(find(eqtable==0))=max(labels(:))+1;
    for E=1:size(eqtable,1)
      eqlabel(E)=min(eqtable(E,:));
    end;  
%    disp(['Relabel the parts of line L with their equivalence class labels'])
    for P=1:NPIXELS
      if th_data(L,P)==1
        [temp_E, temp_F]=find(eqtable==labels(L,P));
        if 1-isempty(temp_E)
          labels(L,P)=eqlabel(temp_E);
        end;
      end;
    end;
    clear eqlabel
  end;
end;

disp(['First top-down loop done in ' num2str(toc)]);
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SECOND LOOP 
%% DOWN-UP (RIGHT-LEFT) SCAN THROUGH LABEL ARRAY
%% TO FIND EQUIVALENT LABELS AND SET UNIFORM LABEL PER REGION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for L=NLINES:-1:1
%  disp(['Initialize local equivalence table for line ' num2str(L)])
  eqtable=[];
%  disp(['Initialize all labels on line ' num2str(L) ' to zero'])
  for P=1:NPIXELS
    if labels(L,P)~=0
      A=neighboursflip_new(L,P,NPIXELS,NLINES,method);
      A=A(find(labels(A)>0)); 
      NL=labels(L,P);   
      for X=1:length(A)
        if labels(A(X))~=NL
          if isempty(eqtable)
            eqtable=NL;
            eqtable(1,2)=labels(A(X));
          else
            [temp_A,temp_B]=find(eqtable==labels(A(X)));
            if isempty(temp_A)
              [temp_M,temp_Z]=find(eqtable==NL);
              if isempty(temp_M)
                eqtable(size(eqtable,1)+1,1)=labels(A(X));
                eqtable(size(eqtable,1),2)=NL;
              else
                temp_O=find(eqtable(temp_M,:)>0);
                eqtable(temp_M,max(temp_O)+1)=labels(A(X));
              end;
            else
              [temp_M,temp_Z]=find(eqtable==NL);
              if isempty(temp_M)
                temp_O=find(eqtable(temp_A,:)>0);
                eqtable(temp_A,max(temp_O)+1)=NL;
              end;
            end;
          end;
        end;
      end;
    end;
  end;

%  disp(['Find equivalence classes detected on this line'])
  
  if isempty(eqtable)
    continue;
  else
    eqtable(find(eqtable==0))=max(labels(:))+1;
    for E=1:size(eqtable,1)
      eqlabel(E)=min(eqtable(E,:));
    end;  
%    disp(['Relabel the parts of line L with their equivalence class labels'])
    for P=1:NPIXELS
      if labels(L,P)~=0
        [temp_E, temp_F]=find(eqtable==labels(L,P));
        if 1-isempty(temp_E)
          labels(L,P)=eqlabel(temp_E);
        end;
      end;
    end;
    clear eqlabel
  end;
end;

disp(['Second down-up loop done in ' num2str(toc)]);
tic

maxnum=0;
maxcell=0;
for ii=1:max(labels(:))
  ind=find(labels==ii);
  if isempty(ind)
    continue;
  elseif length(ind)<minarea
    labels(ind)=0;
    localmax(ind)=0;
  else
    maxnum=maxnum+1;
  end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THIRD LOOP 
%% THROUGH LABEL NUMBERS
%% GENERATE STORM STRUCTURE WITH CHARACTERISTICS
%% FOR EACH LABEL INDEX
%% GENERATE CELL STRUCTURE IF localmax=1 WITHIN STORM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if maxnum>0

  M.S(maxnum).area=[misval];
  M.S(maxnum).meanval=[misval];
  M.S(maxnum).minval=[misval];
  M.S(maxnum).maxval=[misval];
  M.S(maxnum).centroid=[misval misval];
  M.S(maxnum).box=[misval misval misval misval];
  M.S(maxnum).was=[misval];
  M.S(maxnum).life=[misval];
  M.S(maxnum).u = [misval];
  M.S(maxnum).v = [misval];
  M.S(maxnum).parent = [misval];
  M.S(maxnum).child = [misval];
  M.S(maxnum).wasdist = [misval];
  M.S(maxnum).track.xpos = [misval];
  M.S(maxnum).track.ypos = [misval];
  M.S(maxnum).accreted = [misval];
  M.S(maxnum).cell=[misval];
  M.S(maxnum).ellipse.area=[misval];
  M.S(maxnum).ellipse.axialratio=[misval];
  M.S(maxnum).ellipse.orientation=[misval];
  M.S(maxnum).cell=misval;

  maxcell=sum(localmax(:));
  if maxcell>0
    M.C(maxcell).storm=[misval];
    M.C(maxcell).index=[misval];
    M.C(maxcell).maxrain=[misval];
    M.C(maxcell).was=[misval];
  else
    M.C=[];
  end;

else
  M.S=[];
  M.C=[];
end;

disp(['Total number of storms at threshold ' num2str(threshold) ': ' num2str(length(M.S))]);
disp(['Total number of cells: ' num2str(length(M.C))]);
M.label=labels;
M.locmax=localmax;

disp(['Labelling done in ' num2str(toc)]);

return

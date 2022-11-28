function A=neighbourhood(L,P,NPIXELS,NLINES,blockradius)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RETURNS NEIGHBOURHOOD OF PIXEL (L,P)
%% WITHIN BLOCKRADIUS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  A=[];

  if L<=blockradius
    if P<=blockradius
      for P_ind=1:blockradius+P
        A(length(A)+1:length(A)+blockradius+L)=(P_ind-1)*NLINES+[1:blockradius+L];
      end;
    elseif P>NPIXELS-blockradius
      for P_ind=P-blockradius:NPIXELS
        A(length(A)+1:length(A)+blockradius+L)=(P_ind-1)*NLINES+[1:blockradius+L];
      end;
    else
      for P_ind=P-blockradius:P+blockradius
        A(length(A)+1:length(A)+blockradius+L)=(P_ind-1)*NLINES+[1:blockradius+L];
      end;
    end;
  elseif L>NLINES-blockradius
    if P<=blockradius
      for P_ind=1:blockradius+P
        A(length(A)+1:length(A)+NLINES-L+blockradius+1)=(P_ind-1)*NLINES+[L-blockradius:NLINES];
      end;
    elseif P>NPIXELS-blockradius
      for P_ind=P-blockradius:NPIXELS
        A(length(A)+1:length(A)+NLINES-L+blockradius+1)=(P_ind-1)*NLINES+[L-blockradius:NLINES];
      end;
    else
      for P_ind=P-blockradius:P+blockradius
        A(length(A)+1:length(A)+NLINES-L+blockradius+1)=(P_ind-1)*NLINES+[L-blockradius:NLINES];
      end;
    end;
  elseif P<=blockradius
    for P_ind=1:blockradius+P
      A(length(A)+1:length(A)+2*blockradius+1)=(P_ind-1)*NLINES+[L-blockradius:L+blockradius];
    end;
  elseif P>NPIXELS-blockradius
    for P_ind=P-blockradius:NPIXELS
      A(length(A)+1:length(A)+2*blockradius+1)=(P_ind-1)*NLINES+[L-blockradius:L+blockradius];
    end;
  else
    for P_ind=P-blockradius:P+blockradius
      A(length(A)+1:length(A)+2*blockradius+1)=(P_ind-1)*NLINES+[L-blockradius:L+blockradius];
    end;
  end;  

return

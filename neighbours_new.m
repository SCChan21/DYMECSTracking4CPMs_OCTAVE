function A=neighbours_new(L,P,NPIXELS,NLINES,method)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RETURNS NEIGHBOURS OF PIXEL (L,P)
%% ON LINE L AND L-1 ONLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if method==0 %% only connected by sides (4 neighbours)
  if L==1
    if P==1
      A=[(P-1)*NLINES+L];
    else
      A=[(P-1)*NLINES+L; (P-2)*NLINES+L];
    end;
  elseif P==1
    A=[(P-1)*NLINES+L; (P-1)*NLINES+L-1];
  else
    A=[(P-1)*NLINES+L; (P-1)*NLINES+L-1; (P-2)*NLINES+L];
  end;
else % method==1 %% connected through diagonals (8 neighbours)
  if L==1
    if P==1
      A=[(P-1)*NLINES+L];
    else
      A=[(P-1)*NLINES+L; (P-2)*NLINES+L];
    end;
  elseif P==1
    A=[(P-1)*NLINES+L; (P-1)*NLINES+L-1; P*NLINES+L-1];
  elseif P==NPIXELS
    A=[(P-1)*NLINES+L; (P-1)*NLINES+L-1; (P-2)*NLINES+L; (P-2)*NLINES+L-1];
  else
    A=[(P-1)*NLINES+L; (P-1)*NLINES+L-1; (P-2)*NLINES+L; (P-2)*NLINES+L-1; P*NLINES+L-1];
  end;
end

return

function [dx,dy,amp]=fft_track(s1,s2,method)

len=length(s1);
%padlen=2.^ceil(log2(2*len));

if method

%% APPLY TUKEY WINDOW (TAPERED COSINE)

alpha=max([0.1 10/len]);
xhan=[0.5:len-.5];
hann1=ones(size(xhan));
hann1(find(xhan<alpha*len/2))=0.5*(1+cos(pi*(2*xhan(find(xhan<alpha*len/2))/(alpha*len)-1)));
hann1(find(xhan>len*(1-alpha/2)))=0.5*(1+cos(pi*(2*xhan(find(xhan>len*(1-alpha/2)))/(alpha*len)-2/alpha+1)));
hann2=hann1'*hann1;

else

xhan=[0.5:len-.5];
hann1=ones(size(xhan));
hann2=hann1'*hann1;

end

%% FIND CONVOLUTION S1,S2 USING FFT 

b1=s1.*hann2;
b2=s2.*hann2;

m1=b1-mean(b1(:));
m2=b2-mean(b2(:));
 
normval=sqrt(sum(m1(:).^2)*sum(m2(:).^2));
ffv=ifft2(fft2(m2).*conj(fft2(m1)));

[val ind]=max(ffv(:));
[ii,jj]=ind2sub(size(ffv),ind);
dx=ii-1;dy=jj-1; %% ii is change in row so change in y
if dx>len/2
dx=dx-len;
end
if dy>len/2
dy=dy-len;
end
amp=val/normval;

return

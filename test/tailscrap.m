function diff=tailscrap()
% test code for circular convolution with the scrapped portion 
% at the tail of the buffer, rather than the front
%
% The idea is to rotate the zero-padded h (impulse response) buffer
% to the left nh-1 samples, rotating the junk samples as well.
% This could be very handy in avoiding buffer copies during fast filtering.
nh=10;
nfft=256;

#h=ones(1,nh);
h=rand(1,nh);

#x=[1 zeros(1,nfft-1)];
x=rand(1,nfft);

hpad=[ h(nh) zeros(1,nfft-nh) h(1:nh-1) ]; 

y = ifft( fft(hpad) .* fft(x) );
yfilt = filter(h,1,x);
yfilt_no_trans = yfilt(nh:nfft);

#y2=y(nh:nfft);
y2=y(1:nfft-nh+1);
diff=y2 - yfilt_no_trans;
end

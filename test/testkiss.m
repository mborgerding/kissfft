function snr= testkiss( nfft , prec ,scale_t2f ,scale_f2t )
if nargin<1, nfft=1024
endif
if nargin<2, prec='float'
endif
if nargin<3, scale_t2f=1
endif
if nargin<4, scale_f2t=1
endif

if strcmp(prec,'short'),
    kfft='./kffts';
elseif strcmp(prec,'double')
    kfft='./kfftd';
else
    kfft='./kfft';
endif

siglen = nfft;
sig=floor(32767*rand(1,siglen)) + j*floor(32767*rand(1,siglen));

filesave('time.dat',prec,sig);
cmd = sprintf('%s -n %d time.dat freq.dat',kfft,nfft);
system(cmd);

Fsigcomp=fileload('freq.dat',prec,1) * scale_t2f;
Fsig=fft(sig);

%x=linspace(0,2*pi*(nfft-1)/nfft,nfft);
%plot(x,abs(Fsig),'r',x,abs(Fsigcomp),'g')
    diff = Fsig - Fsigcomp;
    noise_pow =  sum( conj(diff).*diff );
    sig_pow = sum( conj(Fsig).*Fsig );
    snr_t2f = 10*log10( sig_pow / noise_pow )
    avg_scale = mean(  abs(Fsig) ./ abs(Fsigcomp) );
    var_scale = var(  abs(Fsig) ./ abs(Fsigcomp) );

cmd = sprintf('%s -i -n %d freq.dat time2.dat',kfft,nfft);
system(cmd);
sigcomp=fileload('time2.dat',prec,1) * scale_f2t;
%sig=ifft(Fsigcomp);

    diff = sig - sigcomp;
    noise_pow =  sum( conj(diff).*diff );
    sig_pow = sum( conj(sig).*sig );
    snr_f2t = 10*log10( sig_pow / noise_pow )
    avg_scale = mean(  abs(sig) ./ abs(sigcomp) );
    var_scale = var(  abs(sig) ./ abs(sigcomp) );

snr=[snr_t2f snr_f2t];    
endfunction


function snr= testkiss( nfft , nrows, prec )
printf('### testing SNR for %d x %d point %s FFTs\n' , nfft,nrows,prec);

if strcmp( prec ,'short')
    scale_t2f=nfft*nrows;
    scale_f2t=nfft*nrows;
else
    scale_t2f=1;
    scale_f2t=1/(nfft*nrows);
endif

kfft= sprintf('./fftutil_%s',prec);

sig=floor(32767*rand(nrows,nfft)) + j*floor(32767*rand(nrows,nfft));

filesave('time.dat',prec,sig);
if nrows > 1
    cmd = sprintf('%s -r %d -n %d time.dat freq.dat',kfft,nrows,nfft);
else    
    cmd = sprintf('%s -n %d time.dat freq.dat',kfft,nfft);
endif 

system(cmd);

Fsigcomp=fileload('freq.dat',prec,1,nrows ) * scale_t2f;

if nrows == 1,
    Fsig=fft(sig);
else    
    Fsig=fft2(sig);
endif

    diff = Fsig - Fsigcomp;
    noise_pow =  sum( conj(diff).*diff );
    sig_pow = sum( conj(Fsig).*Fsig );
    snr_t2f = 10*log10( sig_pow / noise_pow )
    avg_scale = mean(  abs(Fsig) ./ abs(Fsigcomp) );
    var_scale = var(  abs(Fsig) ./ abs(Fsigcomp) );

if nrows > 1
    cmd = sprintf('%s -r %d -i -n %d freq.dat time2.dat',kfft,nrows,nfft);
else    
    cmd = sprintf('%s -i -n %d freq.dat time2.dat',kfft,nfft);
endif

system(cmd);
sigcomp=fileload('time2.dat',prec,1,nrows) * scale_f2t;

    diff = sig - sigcomp;
    noise_pow =  sum( conj(diff).*diff );
    sig_pow = sum( conj(sig).*sig );
    snr_f2t = 10*log10( sig_pow / noise_pow )
    avg_scale = mean(  abs(sig) ./ abs(sigcomp) );
    var_scale = var(  abs(sig) ./ abs(sigcomp) );

snr=[snr_t2f snr_f2t];
endfunction


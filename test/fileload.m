function data = fileload( fname , prec , iscomplex )

f = fopen(fname,"r", "native");
data = fread(f,Inf,prec)';
len=length(data);
fclose(f);

if iscomplex,
    data = (data(1:2:len) + j*data(2:2:len) );
endif    

endfunction

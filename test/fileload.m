function data = fileload( fname , prec , iscomplex ,nrows )

f = fopen(fname,"r", "native");

data = fread(f,Inf,prec);

fclose(f);

len = length(data);
if iscomplex,
    data = (data(1:2:len) + j*data(2:2:len) );
    len = len/2;
endif

tdata = zeros(len/nrows,nrows);
tdata(:) = data;
data = tdata .';

endfunction

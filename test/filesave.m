function filesave( fname , prec , data )

f = fopen(fname,"w", "native");
len=length(data);

if is_complex(data),
    flat=zeros(1,2*len);
    flat(1:2:2*len) = real(data);
    flat(2:2:2*len) = imag(data);
    data = flat;
endif

fwrite(f,data,prec);
fclose(f);

endfunction

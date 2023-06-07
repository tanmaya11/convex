function indices =  supersetIndices(n)
  indices = mat2cell(double(dec2bin(0:2^n-1,n))==48,ones(2^n,1));
end
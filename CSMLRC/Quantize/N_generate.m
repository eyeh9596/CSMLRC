function a = N_generate(n, dim)
a = 1:n;
for i=1:(dim-1)
    a = repmat(a,1,n);
    b = repmat(1:n,n^i,1);
    b = b(:)';
    a = [a;b];
    
end

a = a-1;

c = 0;
for i =1:dim
    c = c+a(i,:)*n^(i-1);
end

a = a+1;

end
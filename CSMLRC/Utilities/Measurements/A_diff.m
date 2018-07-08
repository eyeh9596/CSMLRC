function b = A_diff(x,OMEGA,signvec,P_block,Phi)

N = length(P_block);
block_size = sqrt(N);
n = length(P_image)/(block_size^2);
fx = zeros(N,n);
x = bsxfun(@times,signvec,x(:));
x=reshape(x,[256 256]);
%x=reshape(x,[512 512]);
B=im2col(x,[block_size block_size],'distinct');
B=B(P_block,:);
for i=1:size(B,2)
    B_temp=reshape(B(:,i),[block_size,block_size]);
    B_temp=Phi*B_temp*Phi';
    fx(:,i)=B_temp(:);
end
%fx=Phi*B;
fx=fx';
if iscell(OMEGA)
    for i=1:length(OMEGA)
        b{i}=fx(OMEGA{i});
    end
else
    b=fx(OMEGA);
end

end

% rand('state',1)
% randn('state',1)
% signvec = exp(1i*2*pi*rand(n,1));
% inds=[1:m];
% I=speye(n);
% SubsampM=I(inds(1:m),:);
% k = 3277;
% M=@(x) real(SubsampM(1:k,:)*reshape(fft2(reshape(bsxfun(@times,signvec,x(:)),[height,width])),[n,1]).*(1/sqrt(n)).*sqrt(n/m));
% Mt=@(x) bsxfun(@times,conj(signvec),reshape(ifft2(reshape(SubsampM(1:k,:)'*x(:),[height,width])),[n,1]))*sqrt(n)*sqrt(n/k);
% U=@(x) x(:);
% Ut= @(x) x(:);
% d=ones(m,1)*n/m;
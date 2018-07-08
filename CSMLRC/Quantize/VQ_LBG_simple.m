%%take data columns as the samples
function codebk=VQ_LBG_simple(data, init)

[dim, nn] = size(data);
codebk = init;
flag = 0;
dist2 = 10e8;
mtemp = size(init,2);
while (flag ==0)
    bin = zeros(1,mtemp);
    idx = zeros(1,nn);
    yy = zeros(dim, mtemp);
    for i=1:nn
        dist = sum((codebk(:,1:mtemp)-repmat(data(:,i),[1, mtemp])).^2,1);
        [~,idx(i)] = min(dist);
        bin(idx(i)) = bin(idx(i))+1;
        yy(:,idx(i)) = yy(:,idx(i))+data(:,i);
    end
    codebk(:,1:mtemp) = yy(:,1:mtemp)./repmat(bin,dim,1);
    dist1 = 0;
    for i = 1:nn
        dist1 = dist1+sum((data(:,i)-codebk(:,idx(i))).^2,1);
    end
    dist1 = dist1/(nn*dim);
    if (abs(dist2-dist1)/dist1>0.001)
        flag = 0;
        dist2 = dist1;
    else
        flag = 1;
        
    end

end


return
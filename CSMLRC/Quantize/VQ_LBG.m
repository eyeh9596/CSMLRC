%%take data columns as the samples
function codebk=VQ_LBG(data, cb_size)

[dim, nn] = size(data);
yy = zeros(dim,cb_size);
codebk = mean(data, 2);
mtemp = 1;
while (mtemp < cb_size)
    
    mtemp = mtemp*2;
    for i = 1:mtemp
        j = ceil(i/2);
        if (i == fix(i/2)*2)
            del = -0.001;
        else
            del = 0.001;
        end
        yy(:,i) = codebk(:,j).*(1+del);
    end
    flag = 0;
    dist2 = 10e8;
    while (flag ==0)
        bin = zeros(1,mtemp);
        idx = zeros(1,nn);
        codebk = zeros(dim, mtemp);
        for i=1:nn
            dist = sum((yy(:,1:mtemp)-repmat(data(:,i),[1, mtemp])).^2,1);
            [~,idx(i)] = min(dist);
            bin(idx(i)) = bin(idx(i))+1;
            codebk(:,idx(i)) = codebk(:,idx(i))+data(:,i);
        end
        codebk(:,1:mtemp) = codebk(:,1:mtemp)./repmat(bin,dim,1);
        dist1 = 0;
        for i = 1:nn
            dist1 = dist1+sum((data(:,i)-codebk(:,idx(i))).^2,1);
        end
        dist1 = dist1/(nn*dim);
        if (abs(dist2-dist1)/dist1>0.0005)
            flag = 0;
            dist2 = dist1;
        else
            flag = 1;
            
        end
        
    end
    
end

return
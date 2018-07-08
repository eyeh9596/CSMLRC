function [data_output]=Vector_Quantize(data_input, codebk, dim, flag)

if flag ==1
    nn = fix(length(data_input)/dim);
    data_input = reshape(data_input, [dim, nn]);
    data_output = zeros(1,nn);
    cb_size = size(codebk,2);
    for i=1:nn
        [~,data_output(i)]=min(sum((repmat(data_input(:,i),1,cb_size)-codebk).^2,1));
    end
    
elseif flag ==0
    nn = size(data_input,2);
    data_output = zeros(dim, nn);
    for i=1:nn
        data_output(:,i) = codebk(:,data_input(i));
    end
    data_output = (data_output(:))';
end

end
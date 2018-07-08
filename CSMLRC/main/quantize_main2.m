for i =1:11
    T{i} = (par.y{i}-quantize.Mu(i))./quantize.Sigm(i);
end

T = cell2mat(T);
T = (T(:))';
%T = reshape(T, 2, length(T)/2);


distr = 'Gaussian';
dim = 2;
bit = 5;
cb_size = 2^(bit*dim);

indx = N_generate(32, dim);
for i =1:size(indx,2)
    for j =1:dim
        init(j,i) = y(indx(j,i));
    end
end
%init = init.*0.75;

filename=['..\Quantize\VQ_dim',num2str(dim),'_size',num2str(cb_size),'_',distr,'.mat'];
if exist(filename,'file')
    load(filename);
else
    data_input = randn(1,500000);
    nn = fix(length(data_input)/dim);
    data_input = reshape(data_input, [dim, nn]);
    [IDX,codebk] = kmeans(data_input',cb_size,'start',init');
    codebk = codebk';
    save(filename,'codebk');
end


[data_output]=Vector_Quantize(T, codebk, dim, 1);
[data_input_in]=Vector_Quantize(data_output, codebk, dim, 0);
sum(sum((T-data_input_in).^2))
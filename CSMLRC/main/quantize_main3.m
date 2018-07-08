clear
randn('state',2)
addpath('..\Quantize')
dim = 2;
%% 构建初始码书
bit = 5;
q_l = 2^bit;
distr='Uniform';
filename=['baq_',num2str(q_l),'_',distr,'.mat'];
load(filename);
%y = y./max(y);

indx = N_generate(q_l, dim);
for i =1:size(indx,2)
    for j =1:dim
        init(j,i) = y(indx(j,i));
    end
end
%%

Mu = 1;
init_Norm = zeros(1,size(init,2));
for i =1 :size(init,2)
    init_Norm(:,i) = norm(init(:,i));
    init_Norm2(:,i) = 1*init_Norm(:,i)./(max(abs(init(:,i))));
end
%data_input_Mu = repmat(log(1+Mu.*data_input_Norm)./log(1+Mu),2,1).*data_input./repmat(data_input_Norm2,2,1);

init = init./repmat(init_Norm2,2,1);
figure,plot(init(1,:),init(2,:),'.')
for i =1 :size(init,2)
    init_Norm(:,i) = norm(init(:,i));
end

init_Mu = (repmat(exp(init_Norm).*(1+Mu)-1,dim,1))./Mu.*init;
figure,plot(init_Mu(1,:),init_Mu(2,:),'.')
%%
data_input = randn(1,20000);
cb_size = 2^(bit*dim);
nn = fix(length(data_input)/dim);
data_input = reshape(data_input, [dim, nn]);
%%

for i =1 :200
randn('state',i)
    data_input = randn(1,20000);
    nn = fix(length(data_input)/dim);
    data_input = reshape(data_input, [dim, nn]);
    [IDX,codebk] = kmeans(data_input',cb_size,'start',init');
    codebk = codebk';
end

[data_output]=Vector_Quantize(data_input(:), init, dim, 1);
[data_input_in]=Vector_Quantize(data_output, init, dim, 0);
sum(sum(((data_input(:))'-data_input_in).^2))





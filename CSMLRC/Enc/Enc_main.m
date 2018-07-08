function [PSNR]=Enc_main(Image_name,bpp,quantize)

delete('..\channel\transmit_data.mat')
%%
addpath('..\Utilities');
addpath('..\JPEG');
addpath('..\CCSDS');
addpath('..\JPEG2000');
addpath('..\HEVC');
addpath('..\Utilities\Measurements');
addpath('..\Quantize')
addpath(genpath('..\BCS\BCS-SPL-1.5-1'));
addpath(genpath('..\BCS\BCS-SPL-DPCM-1.0-2'));
addpath('..\BCS\WaveletSoftware');
addpath('..\entropy_code');
addpath(genpath('..\dlwt'));
%% load original images
Test_image_dir = 'D:\cz\Test_Images\256\1\';
ori_im = imread( fullfile(Test_image_dir, Image_name)) ;
ori_im=double(ori_im);

%% 低秩预处理
measure.lowrank_preprocess=1;
par.win        =    6;    % Patch size
par.nblk       =    45;
par.step       =    min(6, par.win-1);
par.beta       =    0.01;
par.nSig       =   4.66;           %  4.45
switch bpp
    case 0.25
        par.c1 =  2;
    case 0.5
        par.c1 =  2;            %  1.96    Threshold for weighted SVT
    case 0.75
        par.c1  =  1;            %  1.96    Threshold for weighted SVT
    case 1
        par.c1 = 0;
    otherwise
        par.c1 = 0;
end

par.c0         =   0;


if measure.lowrank_preprocess==1
    predict_image_dir = '../lowrank_image/1/';
    filename = [predict_image_dir,'low_',num2str(par.c1),'_',Image_name(1:end-4),'.mat'];
    if exist(filename,'file')
        load(filename);
    else
        it = 1;
        flag = 1;
        f = ori_im;
        beta = par.beta;
        for i=1:4
            blk_arr      =     Block_matching( f, par);
            U_arr        =     zeros(par.win^4, size(blk_arr,2), 'single');
            [rim, wei, U_arr]      =   Low_rank_appro(f, par, blk_arr, U_arr, it, flag );
            rim     =    (rim+beta*f)./(wei+beta);
            f = rim;
            
        end
        lowrank_im = f;
        save(filename, 'lowrank_im');
    end
else
    lowrank_im=ori_im;
end
%fprintf('low preprocess PSNR: %f',psnr(lowrank_im, ori_im, 255));
%% Quantization parameter setting
quantize.Mu=[];
quantize.Sigm=[];
quantize.distr='Gaussian';%'Gaussian'or'Uniform','Student','Train'
quantize.method='BAQ';%'BAQ'
quantize.packsize=round(90./quantize.bit);


%% measurement parameter setting
measure.Test_image_dir = fullfile(Test_image_dir, Image_name);
measure.Image_name=Image_name;
measure.dim=2;
measure.model='Hadamard_matrix';%'Hadamard_matrix'or'Orthogonal','Dct_matrix','FFT'
measure.block_size = 16;
measure.N_dc = length(ori_im(:))/(measure.block_size^2);
[measure.image_width,measure.image_height]=size(ori_im);

Rate_proportion=quantize.Rate_proportion;
measure.rate_allocation(1)=measure.N_dc;
measure.rate_allocation(2)=ceil((round(measure.image_width*measure.image_height*bpp/quantize.bit(2)*quantize.dim(2)*Rate_proportion(1))-measure.N_dc)/quantize.dim(2))*quantize.dim(2);
Level=quantize.Level;
for i=3:Level+1
    measure.rate_allocation(i)=ceil((round(measure.image_width*measure.image_height*bpp/quantize.bit(i)*quantize.dim(i)*Rate_proportion(i-1)))/quantize.dim(i))*quantize.dim(i);
end
[measure.prograssive_step]=OMEGA_slip(measure.rate_allocation);
measure.c1 = par.c1;
q=1:length(ori_im(:));
step=measure.prograssive_step;
measure.OMEGA=cell([1,length(step)]);
for i=1:length(measure.prograssive_step)
    measure.OMEGA{i}= q(step(i,1):step(i,2));
end


%% 测量
[par,measure]=CS_encode(lowrank_im,measure);
[par,quantize]=quantize_cell(par,quantize,measure.OMEGA,1);
[par,quantize]=quantize_cell(par,quantize,measure.OMEGA,0);
bin=par.bin(1:2);
Mu=quantize.Mu(1:2);
Sigm=quantize.Sigm(1:2);
% huff_enc = mat2huff(cell2mat(bin(1:2)));
% Trans.data{1}=huff_enc.code;
% Trans.hist{1}=huff_enc.hist;
% quantize.size{1}=huff_enc.size;

%% 预测
y=cell2mat(par.dec(1:2));
%S = length(huff_enc.code)*16+length(huff_enc.hist)*8;
for L=1:Level-1
    
    OMEGA=measure.OMEGA(1:L+1);
    [par]=CS_predict(y, OMEGA, par, measure, quantize.distr, quantize.bit(1:L+1));
    predict_im{L}=par.predict_im;
    PSNR{L} = par.PSNR;
    
    
    %% 残差
    
    OMEGA=cell2mat(measure.OMEGA(L+2));
    A=@(z)A_bp2(z,OMEGA,measure.P_image,measure.P_block,measure.Phi);
    par.y{L+2}=A(lowrank_im)-A(predict_im{L});
    
    %% 量化
    
    
    [par,quantize]=quantize_cell(par,quantize,measure.OMEGA,1);
    [par,quantize]=quantize_cell(par,quantize,measure.OMEGA,0);
    bin(L+2)=par.bin(L+2);
    Mu(L+2)=quantize.Mu(L+2);
    Sigm(L+2)=quantize.Sigm(L+2);
    
    
    %     huff_enc = mat2huff(cell2mat(bin(L+2)));
    %     Trans.data{L+1}=huff_enc.code;
    %     Trans.hist{L+1}=huff_enc.hist;
    %     quantize.size{L+1}=huff_enc.size;
    %     quantize.min=huff_enc.min;
    
    y=[y,A(predict_im{L})+par.dec{L+2}];
    %fprintf('L = %d, Psnr = %f\n',L,psnr(par.predict_im, ori_im,255));
    %S = S+length(huff_enc.code)*16;
    
end
L=L+1;
OMEGA=measure.OMEGA(1:L+1);
[par]=CS_predict(y, OMEGA, par, measure, quantize.distr, quantize.bit(1:L+1));
predict_im{L}=par.predict_im;
PSNR{L} = par.PSNR;
%fprintf('L = %d, Psnr = %f\n',L,psnr(par.predict_im, ori_im,255));
%%
re_code = 'None';
if strcmp(re_code,'None')
    
else
    res_im = uint8(ori_im-round(par.predict_im)+128);
    
    fp = fopen(['../',re_code,'/lynx.raw'],'w');
    fwrite(fp,res_im,'uint8');
    fclose(fp);
    switch re_code
        case 'DLWT'
            res_im = double(res_im);
            coeff = dlwt_mex(res_im,4,0);
            figure(1),imshow(coeff,[])
            
            I_re = dlwt_mex(coeff,4,1);
            figure(2),imshow(I_re,[])
            
            
        case 'JPEG'
            
            y = im2jpeg(res_im,0.6);
            re_dec = jpeg2im(y);
            (size(y.huffman.code,2)*16+size(y.huffman.hist,2)*8+16+8)/measure.image_width/measure.image_height
            fp = fopen('../JPEG/dec.raw');
            fwrite(fp,re_dec,'uint8');
            fclose(fp);
            
        case 'CCSDS'
            
            c_param.in_name = '../CCSDS/lynx.raw';
            c_param.out_name = '../CCSDS/str.bin';
            c_param.height = measure.image_height;
            c_param.width = measure.image_width;
            c_param.Bitrate = 0.25;
            c_param.BitDetpth=8;
            c_param.segnum=512;
            
            ccsds(c_param,1); % CCSDS信源编码
            c_param.in_name = '../CCSDS/str.bin';
            c_param.out_name = '../CCSDS/dec.raw';
            ccsds(c_param,0); % CCSDS信源解码
            
        case 'JPEG2000'
            
            in_name = '../JPEG2000/lynx.raw';
            out_name = '../JPEG2000/str.bin';
            high = measure.image_height;
            width = measure.image_width;
            BitDetpth = 8;
            Bitrate = 0.25;
            
            jpeg_2000(in_name,out_name,1,BitDetpth,Bitrate,width,high) %%JPEG2000信源编码
            in_name='../JPEG2000/str.bin';
            out_name='../JPEG2000/dec.raw';
            jpeg_2000(in_name,out_name,0,BitDetpth,Bitrate,width,high); %%JPEG2000信源解码
            
        case 'HEVC'
            
            in_name = '../HEVC/lynx.raw';
            out_name = '../HEVC/str.bin';
            high = measure.image_height;
            width = measure.image_width;
            %Bitrate = HEVC_Bitrate_parameter(Image_name,1,1);
            Bitrate = round(60000*0.3);
            BitDetpth = 8;
            
            HEVC(in_name,out_name,Bitrate,width,high,BitDetpth,'I',1);%% HEVC信源编码
            D=dir(out_name);
            raptor_compr_size=D.bytes;
            raptor_compr_size*8/high/width
            in_name='../HEVC/str.bin';
            out_name='../HEVC/dec.raw';
            HEVC(in_name,out_name,Bitrate,width,high,BitDetpth,'I',0);%% HEVC信源解码
            
        otherwise
            
    end
    %%
    fp=fopen(['../',re_code,'/dec.raw'],'r');
    re_dec = fread(fp,[measure.image_width, measure.image_height],'uint8');
    fclose(fp);
    
    %psnr(double(res_im),double(re_dec),255)
    mse_value = mean((double(res_im(:)) - double(re_dec(:))).^2 );
    %PSNR = 10*log10( 255*255 / mse_value )
end
%%
% fp = fopen('../Enc/enc.bin','w');
% for L =1:Level+1
%     N = quantize.bit(L);
%     map = cell(2^N,1);
%     for i = 1:2^N
%         map{i} = dec2bin(i-1,N);
%     end
%     hx{L} = cell2mat((map(par.bin{L}+1))');
%     
% end
% hx = cell2mat(hx);
% for i = 1:8:length(hx)
%     if i+7>length(hx)
%         l = length(hx)-i;
%         t = bin2dec([hx(i:i+l),dec2bin(0,8-(l+1))]);
%         fwrite(fp,t,'uint8');
%     else
%         t = bin2dec(hx(i:i+7));
%         fwrite(fp,t,'uint8');
%     end
% end
% fclose(fp);
%%
% fp = fopen('../Enc/enc.bin','r');
% hx_dec = fread(fp,'uint8');
%
%
% fclose(fp)

%%
Trans.Mu=Mu;
Trans.Sigm=Sigm;
Trans.data = par.bin;
quantize=rmfield(quantize,'Mu');
quantize=rmfield(quantize,'Sigm');

save('..\channel\transmit_data.mat','Trans','quantize','measure')


B_all = 0;
for i =1:length(par.y)
    B_all = B_all+length(par.bin{i})*quantize.bit(i);
end
B_all/(measure.image_width*measure.image_height);

end
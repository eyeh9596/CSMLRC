function  [Re]    =  CS_decode(Trans,quantize,measure)
randn('state',0);
rand('state',0);
ori_im=imread( measure.Test_image_dir);
ori_im=double(ori_im);% For computing PSNR only
Level=quantize.Level;
%%
% bin=Trans.data;
% huff_enc.min=quantize.min;
% huff_enc.code=bin{1};
% huff_enc.hist=Trans.hist{1};
% huff_enc.size=quantize.size{1};
% data_temp=huff2mat(huff_enc);
% Trans.data{1}=data_temp(measure.OMEGA{1});
% Trans.data{2}=data_temp(measure.OMEGA{2});
% for L=2:Level
%
%     huff_enc.code=bin{L};
%     huff_enc.hist=Trans.hist{L};
%     huff_enc.size=quantize.size{L};
%     Trans.data{L+1}=huff2mat(huff_enc);
% end

%% ≤‚ ‘
par.bin=Trans.data;
quantize.Mu=Trans.Mu;
quantize.Sigm=Trans.Sigm;
[par,quantize]=quantize_cell(par,quantize,measure.OMEGA,0);
y=cell2mat(par.dec(1:2));
for L=1:Level-1
    
    OMEGA=measure.OMEGA(1:L+1);
    [par]=CS_predict(y, OMEGA, par, measure, quantize.distr, quantize.bit(1:L+1));
    predict_im{L}=par.predict_im;
    
    
    %% ≤–≤Ó
    
    OMEGA=cell2mat(measure.OMEGA(L+2));
    A=@(z)A_bp2(z,OMEGA,measure.P_image,measure.P_block,measure.Phi);
    par.dec{L+2}=par.dec{L+2}+A(predict_im{L});
    
    %% ¡øªØ
    
    y=[y,A(predict_im{L})+par.dec{L+2}];
    %fprintf('L = %d, Psnr = %f\n',L,psnr(par.predict_im, ori_im,255));
    
    
end
L=L+1;
OMEGA=measure.OMEGA(1:L+1);
[par]=CS_predict(y, OMEGA, par, measure, quantize.distr, quantize.bit(1:L+1));
predict_im{L}=par.predict_im;
%%
rec_im=predict_im{L};
PSNR_nlr     =   csnr( rec_im, ori_im, 0, 0 );
%SSIM_nlr = cal_ssim(rec_im, ori_im,0,0);
SSIM_nlr = msssim(rec_im, ori_im);
%%


%%
Re.Rec_im=rec_im;
Re.PSNR_nlr=PSNR_nlr;
Re.SSIM_nlr = SSIM_nlr;


return;





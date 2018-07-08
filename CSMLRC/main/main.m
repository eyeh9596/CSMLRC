addpath('..\Enc');
addpath('..\Dec');
addpath('..\channel');
%Im_name_all={'barbara.tif','lena.tif','Monarch.tif','peppers.tif','boats.tif','parrots.tif','foreman.tif','house.tif'};
%P=1:9;
%quantize.bit = [5,5,3,3,2,2,2,1,1,1,1,1];
% quantize.bit = [5,5,5,4,4,3,3,2,2,1,1,2,3,3,4,4,4,3,3,3,2,2,2];
% quantize.dim = [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2];
quantize.bit = [5,5,4,4,3,3,2,2,1,1,1];%case bpp = 1.00
quantize.bit = [5,5,4,4,3,3,2,2,1,1,1];%case bpp = 0.75
quantize.bit = [5,5,3,3,3,2,2,2,1,1,1];%case bpp = 0.50
quantize.bit = [5,5,3,3,3,2,2,2,1,1,1];%case bpp = 0.25
quantize.dim = [1,1,1,1,1,1,1,1,1,1,1];
%  quantize.bit=[5,5,5,3,3,3,3,2,2,1,1,1];


for k=1:8
    for i=1:4
        bpp=0.25*i;
        if (bpp>0.60)
            quantize.bit = [5,5,4,4,3,3,2,2,1,1,1];%case bpp = 1.00
        elseif bpp <= 0.60
            quantize.bit = [5,5,3,3,3,2,2,2,1,1,1];%case bpp = 0.50
        end

        for m=1
            for n=1
                Image_name = [num2str(k),'.tif'];
                %quantize.Rate_proportion=[0.05+(m-1)*0.05,0.15+0.05*n,1-0.05-(m-1)*0.05-0.15-0.05*n]; %渐进分配bit
                quantize.Level=10;
                quantize.Rate_proportion=ones(1,quantize.Level)*(1/quantize.Level);
                PSNR=Enc_main(Image_name,bpp,quantize);
                
                %packet_lost((n-1)*0.05);
                [Re]=Dec_main();
               % PSNR_sum{k}(i,:)=cell2mat(PSNR);
                load('..\channel\transmit_data.mat');
                %%记录结果
                bits_num=length(cell2mat(Trans.data))*16/measure.image_width/measure.image_height;
                fp=fopen(['..\results\1\',Image_name(1:end-4),'_NLR_CS','.txt'],'a');
                fprintf(fp,'image:%s; bpp:%f; NLR-PSNR:%f; NLR-SSIM:%f; Proportion:%s; bits:%s\n',Image_name,bpp,Re.PSNR_nlr,Re.SSIM_nlr,strrep(num2str(quantize.Rate_proportion),'        ','-'),strrep(num2str(quantize.bit),'  ','_'));
                fprintf('image:%s; bpp:%f; NLR-PSNR:%f; NLR-SSIM:%f; Proportion:%s; bits:%s\n',Image_name,bpp,Re.PSNR_nlr,Re.SSIM_nlr,strrep(num2str(quantize.Rate_proportion),'        ','_'),strrep(num2str(quantize.bit),'  ','_'));
                fclose(fp);
                Test_image_dir = 'D:\cz\Test_Images\256\1\';
                ori_im = imread( fullfile(Test_image_dir, Image_name)) ;
                %pip_imshow(ori_im,48,68);
                %pip_imshow(Re.Rec_im,48,68);%202,146
            end
        end
    end
end


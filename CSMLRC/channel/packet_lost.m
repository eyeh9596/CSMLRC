function []=packet_lost(PLR)

	load('..\channel\transmit_data.mat');
    l1=0;
    l=length(Trans.data);
    if isempty(Trans.data{end})
        l=length(Trans.data)-1;
    end
    for i=1:l
        l1=l1+ceil(length(Trans.data{i})/quantize.packsize(i));
    end
    
	for i=2:l
		v=Trans.data{i};
		[~,~,v_lossmark]=setpackloss(v,length(v),quantize.packsize(i),PLR); %N_short��ʾ��ȥ���Ĳ���ֵʣ�µĲ�������,%v_lossmarkΪ���յ��Ĳ���ֵ�ı�־��1Ϊ������0Ϊ���յ���		
		[n,~]=find(v_lossmark>0);
		Trans.data{i}(n)=[];
		measure.OMEGA{i}([quantize.dim(i)*n;quantize.dim(i)*n-1])=[];%quantize.dim Ϊ2

    end
    l2=0;
    for i=1:l
        l2=l2+ceil(length(Trans.data{i})/quantize.packsize(i));
    end
    
    quantize.PLR=1-l2/l1;
    save('..\channel\transmit_data.mat','Trans','quantize','measure')

end
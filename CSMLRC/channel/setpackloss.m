


function [v_short,N_short,v_lossmark]=setpackloss(v_bin,N,packsize,PLR)


n=ceil(N/packsize);                       %���Ĵ�С������
num=randperm(n);
xx=ceil(n*PLR);
lossmark=zeros(n,1);
v_lossmark=zeros(N,1); 
 for i=1:xx
    lossmark(num(i))=1;            %������λ����1��ʾ��δ������0��ʾ
 end

  for i=1:N
     if(lossmark(ceil(i/packsize))==1)
         v_lossmark(i)=1;
     end 
 end
 v_short=v_bin;
 N_short=N-xx*packsize;
 if(lossmark(n)==1)                    %��Ϊ���һ������ֻ��4������ֵ��������һ�����Ķ�����־Ϊ�棬ֻ�ܼ�ȥ4���������Ǽ�ȥ20��
     N_short=N_short+(n*packsize-N);
 end

 
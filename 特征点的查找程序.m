wavename='sym8';
m=val(1,1:k);
m1=(m-1024)./200;
E=m1';
[C L]=wavedec(E,3,wavename);
cA3=appcoef(C,L,wavename,3);
cD1=detcoef(C,L,1);
cD2=detcoef(C,L,2);
cD3=detcoef(C,L,3);
[thr,sorh,keepapp]=ddencmp('den','wv',E);
thr1=thr./log(2);
thr2=thr./log(3);
thr3=thr./log(4);
TR=[thr1,thr2,thr3];
[XC,CXC,LXC,PERF0,PERF2]=wdencmp('lvd',E,wavename,3,TR,sorh);
ecgdata=XC;
level=4;    sr=360; 
swa=zeros(4,k);%�洢��ò��Ϣ
swd=zeros(4,k);%�洢ϸ����Ϣ
signal=ecgdata(0*k+1:1*k); %ȡ���ź�
for i=1:k-3
   swa(1,i+3)=1/4*signal(i+3-2^0*0)+3/4*signal(i+3-2^0*1)+3/4*signal(i+3-2^0*2)+1/4*signal(i+3-2^0*3);
   swd(1,i+3)=-1/4*signal(i+3-2^0*0)-3/4*signal(i+3-2^0*1)+3/4*signal(i+3-2^0*2)+1/4*signal(i+3-2^0*3);
end
j=2;
while j<=level
   for i=1:k-24
     swa(j,i+24)=1/4*swa(j-1,i+24-2^(j-1)*0)+3/4*swa(j-1,i+24-2^(j-1)*1)+3/4*swa(j-1,i+24-2^(j-1)*2)+1/4*swa(j-1,i+24-2^(j-1)*3);
     swd(j,i+24)=-1/4*swa(j-1,i+24-2^(j-1)*0)-3/4*swa(j-1,i+24-2^(j-1)*1)+3/4*swa(j-1,i+24-2^(j-1)*2)+1/4*swa(j-1,i+24-2^(j-1)*3);
   end
   j=j+1;
end
%����ԭ�źźͳ߶�ϵ����С��ϵ��
%figure(10);
%subplot(level+1,1,1);plot(ecgdata(1:k),'b');grid on ;axis tight;
%title('ECG�ź���j=1,2,3,4�߶��µĳ߶�ϵ���Ա�');
%for i=1:level
%    subplot(level+1,1,i+1);
%    plot(swa(i,:),'b');axis tight;grid on; xlabel('time');ylabel(strcat('a  ',num2str(i)));
%end
%figure(11);
%subplot(level,1,1); plot(ecgdata(1:k),'b'); grid on;axis tight;
%title('ECG�źż�����j=1,2,3,4�߶��µ�С��ϵ��');
%for i=1:level
%    subplot(level+1,1,i+1);
%    plot(swd(i,:),'b'); axis tight;grid on;
%    ylabel(strcat('d   ',num2str(i)));
%end

%**************************************����������ֵ��**********************%
ddw=zeros(size(swd));
pddw=ddw;
nddw=ddw;
%С��ϵ���Ĵ���0�ĵ�
posw=swd.*(swd>0);
%б�ʴ���0
pdw=((posw(:,1:k-1)-posw(:,2:k))<0);%�ҳ�б�ʴ���0��ֵ������ֵΪ1������Ϊ0
%������ֵ��
pddw(:,2:k-1)=((pdw(:,1:k-2)-pdw(:,2:k-1))>0);%�ҵ�����ֵ����Ϊ1������Ϊ0
%С��ϵ��С��0�ĵ�
negw=swd.*(swd<0);
ndw=((negw(:,1:k-1)-negw(:,2:k))>0);%�ҳ�б��С��0��ֵ������ֵΪ1������Ϊ0
%������ֵ��
nddw(:,2:k-1)=((ndw(:,1:k-2)-ndw(:,2:k-1))>0);%�ҵ���Сֵ����Ϊ1������Ϊ0
%������
ddw=pddw|nddw;
ddw(:,1)=1;
ddw(:,k)=1;
%�����ֵ���ֵ,��������0
wpeak=ddw.*swd;
wpeak(:,1)=wpeak(:,1)+1e-10;%???
wpeak(:,k)=wpeak(:,k)+1e-10;%???
%�������߶��¼�ֵ��
%figure(12);
%for i=1:level
%   subplot(level,1,i);
%    plot(wpeak(i,:),'b'); axis tight;grid on;
%ylabel(strcat('j=   ',num2str(i)));
%end
%subplot(4,1,1);
%title('ECG�ź���j=1,2,3,4�߶��µ�С��ϵ����ģ����ֵ��');
interva2=zeros(1,k);
intervaqs=zeros(1,k);
Mj1=wpeak(1,:);
Mj3=wpeak(3,:);
Mj4=wpeak(4,:);
%�����߶�3��ֵ��
%figure(13);
%plot (Mj3,'b');
%title('�߶�3��С��ϵ����ģ����ֵ��');
%��������ֵ��ƽ��
posi=Mj3.*(Mj3>0);
thposi=(max(posi(1:round(k/4)))+max(posi(round(k/4):2*round(k/4)))+max(posi(2*round(k/4):3*round(k/4)))+max(posi(3*round(k/4):4*round(k/4))))/4;
posi=(posi>thposi/3);
%�󸺼���ֵ��ƽ��
nega=Mj3.*(Mj3<0);
thnega=(min(nega(1:round(k/4)))+min(nega(round(k/4):2*round(k/4)))+min(nega(2*round(k/4):3*round(k/4)))+min(nega(3*round(k/4):4*round(k/4))))/4;
nega=-1*(nega<thnega/4);
%�ҳ���0��
interva=posi+nega;
loca=find(interva);
for i=1:length(loca)-1%loca����
    if abs(loca(i)-loca(i+1))<80%����Ǿ���ֵ��
       diff(i)=interva(loca(i))-interva(loca(i+1));%С��80����ĵ㸳ֵΪ-2
    else
       diff(i)=0;
    end
end
%�ҳ���ֵ��
loca2=find(diff==-2);
%������ֵ��
interva2(loca(loca2(1:length(loca2))))=interva(loca(loca2(1:length(loca2))));
%������ֵ��
interva2(loca(loca2(1:length(loca2))+1))=interva(loca(loca2(1:length(loca2))+1));%ʣ�µĶ�������ֵ
intervaqs(1:k-10)=interva2(11:k);
countR=zeros(1,1);
countQ=zeros(1,1);
countS=zeros(1,1);
mark1=0;
mark2=0;
mark3=0;
i=1;
j=1;
Rnum=0;
%*************************��������ֵ�Թ��㣬��R����ֵ��������QRS����㼰�յ�*******************%
while i<k
    if interva2(i)==-1
       mark1=i;
       i=i+1;
       while(i<k&interva2(i)==0)
          i=i+1;
       end
       mark2=i;
%�󼫴�ֵ�ԵĹ����
       mark3= round((abs(Mj3(mark2))*mark1+mark2*abs(Mj3(mark1)))/(abs(Mj3(mark2))+abs(Mj3(mark1))));
%R������ֵ��
       R_result(j)=mark3-10;%Ϊ��-10������ֵ��
       countR(mark3-10)=1;
%���QRS�����
       kqs=mark3-10;
       markq=0;
     while (kqs>1)&&( markq< 3)
         if Mj1(kqs)~=0
            markq=markq+1;
         end
         kqs= kqs -1;
     end
  countQ(kqs)=-1;
  
%���QRS���յ�  
  kqs=mark3-10;
  marks=0;
  while (kqs<k)&&( marks<3)
      if Mj1(kqs)~=0
         marks=marks+1;
      end
      kqs= kqs+1;
  end
  countS(kqs)=-1;
  i=i+60;
  j=j+1;
  Rnum=Rnum+1;
 end
i=i+1;
end
%************************ɾ�����㣬����©���**************************%
num2=1;
while(num2~=0)
   num2=0;
%j=3,�����
   R=find(countR);
   Q=find(countQ);
   S=find(countS);
%�������
   R_R=R(2:length(R))-R(1:length(R)-1);
   RRmean=mean(R_R);
%������R�����С��0.4RRmeanʱ,ȥ��ֵС��R��
for i=2:length(R)
    if (R(i)-R(i-1))<=0.4*RRmean
        num2=num2+1;
        if signal(R(i))>signal(R(i-1))
            countR(R(i-1))=0;
			countQ(Q(i-1))=0;
			countS(S(i-1))=0;
        else
            countR(R(i))=0;
			countQ(Q(i))=0;
			countS(S(i))=0;
        end
    end
end
end

num1=2;
while(num1>0)
   num1=num1-1;
   R=find(countR);
   R_R=R(2:length(R))-R(1:length(R)-1);
   RRmean=mean(R_R);
%������R���������1.6RRmeanʱ,��С��ֵ,����һ�μ��R��
for i=2:length(R)
    if (R(i)-R(i-1))>1.6*RRmean
        Mjadjust=wpeak(4,R(i-1)+80:R(i)-80);
        points2=(R(i)-80)-(R(i-1)+80)+1;
%��������ֵ��
        adjustposi=Mjadjust.*(Mjadjust>0);
        adjustposi=(adjustposi>thposi/4);
%�󸺼���ֵ��
        adjustnega=Mjadjust.*(Mjadjust<0);
        adjustnega=-1*(adjustnega<thnega/5);
%������
        interva4=adjustposi+adjustnega;
%�ҳ���0��
        loca3=find(interva4);
        diff2=interva4(loca3(1:length(loca3)-1))-interva4(loca3(2:length(loca3)));
%����м���ֵ��,�ҳ�����ֵ��
        loca4=find(diff2==-2);
        interva3=zeros(points2,1)';
        for j=1:length(loca4)
           interva3(loca3(loca4(j)))=interva4(loca3(loca4(j)));
           interva3(loca3(loca4(j)+1))=interva4(loca3(loca4(j)+1));
        end
        mark4=0;
        mark5=0;
        mark6=0;
    while j<points2
         if interva3(j)==-1;
            mark4=j;
            j=j+1;
            while(j<points2&interva3(j)==0)
                 j=j+1;
            end
            mark5=j;
%������
            mark6= round((abs(Mjadjust(mark5))*mark4+mark5*abs(Mjadjust(mark4)))/(abs(Mjadjust(mark5))+abs(Mjadjust(mark4))));
            countR(R(i-1)+80+mark6-10)=1;
            j=j+60;
         end
         j=j+1;
     end
    end
 end
end
%����ԭͼ����������
Mj4posi=Mj4.*(Mj4>0);
%��������ֵ��ƽ��
Mj4thposi=(max(Mj4posi(1:round(k/4)))+max(Mj4posi(round(k/4):2*round(k/4)))+max(Mj4posi(2*round(k/4):3*round(k/4)))+max(Mj4posi(3*round(k/4):4*round(k/4))))/4;
Mj4posi=(Mj4posi>Mj4thposi/3);
Mj4nega=Mj4.*(Mj4<0);
%�󸺼���ֵ��ƽ��
Mj4thnega=(min(Mj4nega(1:round(k/4)))+min(Mj4nega(round(k/4):2*round(k/4)))+min(Mj4nega(2*round(k/4):3*round(k/4)))+min(Mj4nega(3*round(k/4):4*round(k/4))))/4;
Mj4nega=-1*(Mj4nega<Mj4thnega/4);
Mj4interval=Mj4posi+Mj4nega;
Mj4local=find(Mj4interval);
Mj4interva2=zeros(1,k);
for i=1:length(Mj4local)-1
    if abs(Mj4local(i)-Mj4local(i+1))<80
       Mj4diff(i)=Mj4interval(Mj4local(i))-Mj4interval(Mj4local(i+1));
    else
       Mj4diff(i)=0;
    end
end
%�ҳ���ֵ��
Mj4local2=find(Mj4diff==-2);
%������ֵ��
Mj4interva2(Mj4local(Mj4local2(1:length(Mj4local2))))=Mj4interval(Mj4local(Mj4local2(1:length(Mj4local2))));
%������ֵ��
Mj4interva2(Mj4local(Mj4local2(1:length(Mj4local2))+1))=Mj4interval(Mj4local(Mj4local2(1:length(Mj4local2))+1));
mark1=0;
mark2=0;
mark3=0;
Mj4countR=zeros(1,1);
Mj4countQ=zeros(1,1);
Mj4countS=zeros(1,1);
flag=0;
while i<k
    if Mj4interva2(i)==-1
       mark1=i;
       i=i+1;
       while(i<k&Mj4interva2(i)==0)
          i=i+1;
       end
       mark2=i;
%�󼫴�ֵ�ԵĹ����,��R4�м�ֵ֮���������R�㡣
       mark3= round((abs(Mj4(mark2))*mark1+mark2*abs(Mj4(mark1)))/(abs(Mj4(mark2))+abs(Mj4(mark1))));
       Mj4countR(mark3)=1;
       Mj4countQ(mark1)=-1;
       Mj4countS(mark2)=-1;
       flag=1;
    end
    if flag==1
        i=i+200;
        flag=0;
    else
        i=i+1;
    end
end
%%%%%�Գ߶�4��R���ⲻ���ã���Ҫ�Ľ��ĵط�
%figure(20);
%plot(Mj4,'b');
%title('j=4');
%hold on;
%plot(Mj4countR,'r');
%plot(Mj4countQ,'g');
%plot(Mj4countS,'g');
%%%%%%%%%%%%%%%%%%%%%%%%%%Mj4������ҵ�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rlocated=find(Mj4countR);
Qlocated=find(Mj4countQ);
Slocated=find(Mj4countS);
countPMj4=zeros(1,1);
countTMj4=zeros(1,1);
countP=zeros(1,1);
countPbegin=zeros(1,1);
countPend=zeros(1,1);
countT=zeros(1,1);
countTbegin = zeros(1,1);
countTend = zeros(1,1);
windowSize=100;
%%%%%%%%%%%%%%%%%%%%%%%P�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rlocated Qlocated ���ڳ߶�4�µ�����
for i=2:length(Rlocated)
    flag=0;
    mark4=0;
    RRinteral=Rlocated(i)-Rlocated(i-1);
    for j=1:5:(RRinteral*2/3)
       % windowEnd=Rlocated(i)-30-j;
       windowEnd=Qlocated(i)-j;
        windowBegin=windowEnd-windowSize;
        if windowBegin<Rlocated(i-1)+RRinteral/3
            break;
        end
        %���ڵļ���Сֵ
        %windowposi=Mj4.*(Mj4>0);
        %windowthposi=(max(Mj4(windowBegin:windowBegin+windowSize/2))+max(Mj4(windowBegin+windowSize/2+1:windowEnd)))/2;
        [windowMax,maxindex]=max(Mj4(windowBegin:windowEnd));
        [windowMin,minindex]=min(Mj4(windowBegin:windowEnd));
        if minindex < maxindex &&((maxindex-minindex)<windowSize*2/3)&&windowMax>0.01&&windowMin<-0.1
            flag=1;
            mark4=round((maxindex+minindex)/2+windowBegin);
            countPMj4(mark4)=1;
            countP(mark4-20)=1;
            countPbegin(windowBegin+minindex-20)=-1;
            countPend(windowBegin+maxindex-20)=-1;
        end
        if flag==1
            break; 
        end
    end
    if mark4==0&&flag==0 %���û��P������R������1/3����ֵ-1
        mark4=round(Rlocated(i)-RRinteral/3);
        countP(mark4-20)=-1;
    end
end
%plot(countPMj4,'g');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%T�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear windowBegin windowEnd maxindex minindex windowMax windowMin mark4 RRinteral;

windowSizeQ=100;
for i=1:length(Rlocated)-1;
    flag=0;
    mark5=0;
    RRinteral=Rlocated(i+1)-Rlocated(i);
    for j=1:5:(RRinteral*2/3)
       % windowBegin=Rlocated(i)+30+j;
       windowBegin=Slocated(i)+j;
        windowEnd  =windowBegin+windowSizeQ;
        if windowEnd >Rlocated(i+1)-RRinteral/4
            break;
        end
        %%%%%�󴰿��ڵ������Сֵ
        [windowMax,maxindex]=max(Mj4(windowBegin:windowEnd));
        [windowMin,minindex]=min(Mj4(windowBegin:windowEnd));
        if minindex < maxindex &&((maxindex-minindex)<windowSizeQ)&&windowMax>0.1&&windowMin<-0.1
            flag=1;
            mark5=round((maxindex+minindex)/2+windowBegin);
            countTMj4(mark5)=1;
            countT(mark5-20)=1;%�ҵ�T����ֵ��
            %%%%%ȷ��T����ʼ����յ�
            countTbegin(windowBegin+minindex-20)=-1;
            countTend(windowBegin+maxindex-20)=-1;
        end
        if flag==1
            break;
        end
    end
    if mark5==0 %���û��T������R���� ���1/3����ֵ-2
        mark5=round(Rlocated(i)+ RRinteral/3);
        countT(mark5)=-2;
    end
end
%plot(countTMj4,'g');
%hold off;        
figure(21);
plot(ecgdata(0*k+1:1*k),'b'),grid on,axis tight,axis([1,k,-2,5]);
title('ECG�źŵĸ������μ��');
hold on
%plot(countR,'r');
%plot(countQ,'k');
%plot(countS,'k');
R_Q=R-Q;
S_R=S-R;
for i=1:length(R)
windowBegin=R(i)-15;
windowEnd=R(i)+15;
[windowMax,maxindex(i)]=max(XC(windowBegin:windowEnd));
maxindex(i)=R(i)+maxindex(i)-16;
end
R_result=maxindex;
windowBegin=1;
for i=1:length(Q)
windowBegin=Q(i)-R_Q(i);
if windowBegin<=1
windowBegin=1 ; 
windowEnd=Q(i)+R_Q(i);
[windowMin,minindex(i)]=min(XC(windowBegin:windowEnd));
minindex(i)=minindex(i);
else
windowBegin=Q(i)-R_Q(i);
windowEnd=Q(i)+R_Q(i);
[windowMin,minindex(i)]=min(XC(windowBegin:windowEnd));
minindex(i)=Q(i)+minindex(i)-R_Q(i)-1;
end
end
Q_result=minindex;

for i=1:length(S)
windowBegin=S(i)-S_R(i);
windowEnd=ceil(S(i)+5/4*S_R(i));
if windowEnd >= k
windowEnd=k;
end
[windowMin,minindex(i)]=min(XC(windowBegin:windowEnd));
minindex(i)=S(i)+minindex(i)-S_R(i)-1;
end
S_result=minindex;
for i=1:Rnum
    if R_result(i)==0;
        break
    end
    plot(R_result(i),ecgdata(R_result(i)),'bo','MarkerSize',10,'MarkerEdgeColor','g');
end
for i=1:Rnum
    if Q_result(i)==0;
        break
    end
    plot(Q_result(i),ecgdata(Q_result(i)),'bo','MarkerSize',10,'MarkerEdgeColor','r');
end
for i=1:Rnum
    if S_result(i)==0;
        break
    end
    plot(S_result(i),ecgdata(S_result(i)),'bo','MarkerSize',10,'MarkerEdgeColor','k');
end
plot(countP,'r');
plot(countT,'k');
plot(countPbegin,'g');
plot(countPend,'g');
plot(countTbegin,'k');
plot(countTend,'k
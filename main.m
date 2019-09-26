%�˳���Ϊ����������ͼ��λģ���㷨ģ�����
%���Ǽ����Ǳ��е����ݶ�����С�ı�׼ֵ
%�������f=1����������Ϊ�뾶Ϊ1������
%Ȼ��ȡ��С�̶����񳡣������������ȡ����һ�ǵ㣬��ΪĿ���ǵ�
%�ò����м������в�������Ϊ���������Ҫ����Ŀ��Ҫ������в������ݼ���һ��������Ŷ�
clear all
load('����2  �����Ǳ�.mat')

star=zeros(4908,3);%���Ǳ��е����ݴ�������ϵת��Ϊֱ������ϵ
star(:,1)=cos(star_data(:,3)*pi/180).*cos(star_data(:,2)*pi/180);
star(:,2)=cos(star_data(:,3)*pi/180).*sin(star_data(:,2)*pi/180);
star(:,3)=sin(star_data(:,3)*pi/180);
D=randsample(star_data(:,1),1);%���Ǳ������ѡȡһ���ǵ���ΪĿ���
[A,I]=sort(star_data(:,3));%�����Ǳ��γ��˳���Ǵ��ҵģ�Ҫ�������Ա�ȶԺ�ɸѡ
vi=12;%����Ŀ������
V=24;%��������������Χ����λΪ��
N=20;%�趨ѡȡ�ǵ����Ŀ
%% ��ͬ������Ŀ�����ɸѡ
if (star_data(end,2)-star_data(D,2))<vi %������ѡ����ǵ�λ���Ǳ�ྭһ�е�ĩ�ˣ������Ǳ����һ���ǳྭ����vi
    D1=vi-(star_data(end,2)-star_data(D,2));
    i=1;
    j=1;
    while star_data(i,2)<=D1
        i=i+1;
    end
    while (star_data(D,2)-star_data(D-j,2))<=vi
        j=j+1;
    end
    if (A(end,1)-star_data(D,3))<vi %������ѡ����ǵ�λ���Ǳ��γһ�е�ĩ�ˣ������Ǳ����һ���ǳ�γ����vi
        D2=V-(A(end,1)-star_data(D,3));
        m=1;
        [~,Da]=intersect(I,D);%�ҵ��Ǳ��γ��Ϣ����֮��Ŀ���ǵ��λ��
        n=4908-Da;
        while (A(Da,1)-A(Da-m,1))<=D2
            m=m+1;
        end
        G1=[D+j:4908,1:i];
        G=intersect(G1,I(Da-m:end,1));%Ѱ��ͬʱ�����γ��Χ�ͳྭ��Χ���ǵ㣬��Ѱ�Ҷ��߽���
        if size(G,1)>N %�����ȡһ�������ǵ�
          G=randsample(G,N);
        end
        [p1,p2]=projection(D,G,star); %����projection
        [pre,epsion,perm]=callo(G,D,star);%����callo
        figure
        plot3(p1(:,1),p1(:,2),p1(:,3),'w.','MarkerSize',0.5)%����ѡȡ�ǵ���������ʾ��ͼ
        set(gca,'color',[0,0,0])
        title('��������ǵ���άʾ��ͼ')
        figure
        plot(p2(:,1),p2(:,2),'w.','MarkerSize',0.5)%����ѡȡ�ǵ�ͶӰ�ڶ�άƽ���ʾ��ͼ
        hold on
        plot(0,0,'r*') %����Ŀ���
        set(gca,'color',[0,0,0])
        title('��������ǵ��άͶӰʾ��ͼ')
    end
    if (star_data(D,3)-A(1,1))<vi %������ѡ����ǵ�λ���Ǳ��γһ�еĿ��ˣ������Ǳ��һ���ǳ�γ����vi
       D3=V-(star_data(D,3)-A(1,1));
       m1=1;
       [~,Da]=intersect(I,D); 
       n1=Da;
       while (A(Da+m1,1)-star_data(D,3))<=D3
          m1=m1+1;
       end
       G1=[D+j:4908,1:i];
       G=intersect(G1,I(1:(Da+m1),1));
       if size(G,1)>N
          G=randsample(G,N);
       end
       [p1,p2]=projection(D,G,star);
       [pre,epsion,perm]=callo(G,D,star);
       figure
       plot3(p1(:,1),p1(:,2),p1(:,3),'w.','MarkerSize',0.5)
       set(gca,'color',[0,0,0])
       title('��������ǵ���άʾ��ͼ') 
       figure
       plot(p2(:,1),p2(:,2),'w.','MarkerSize',0.5)
       hold on
       plot(0,0,'r*') %����Ŀ���
       set(gca,'color',[0,0,0])
       title('��������ǵ��άͶӰʾ��ͼ')
    end
    % ����ͬ���沽������������ȡ������жϣ���ִ����ͬ�Ĳ���
    [~,Da]=intersect(I,D); 
    m=1;
    n=1;
    while (A(Da+m,1)-A(Da,1))<=vi
        m=m+1;
    end
    while (A(Da,1)-A(Da-n,1))<=vi
        n=n+1;
    end
    G1=[D+j:4908,1:i];
    G=intersect(G1,I(Da-n:Da+m,1));
    if size(G,1)>N
      G=randsample(G,N);
    end
    [p1,p2]=projection(D,G,star);
    [pre,epsion]=callo(G,D,star);
    figure
    plot3(p1(:,1),p1(:,2),p1(:,3),'w.','MarkerSize',0.5)
    set(gca,'color',[0,0,0])
    title('��������ǵ���άʾ��ͼ') 
    figure
    plot(p2(:,1),p2(:,2),'w.','MarkerSize',0.5)
    hold on
    plot(0,0,'r*') %����Ŀ���
    set(gca,'color',[0,0,0])
    title('��������ǵ��άͶӰʾ��ͼ')
end

if (star_data(D,2)-star_data(1,2))<vi
    D1=vi-(star_data(D,2)-star_data(1,2));
    i=1;
    j=1;
    while (star_data(end,2)-star_data(end-i,2))<=D1
        i=i+1;
    end
    while star_data(D+j,2)-(star_data(D,2))<=vi
        j=j+1;
    end
    if (A(end,1)-star_data(D,3))<vi
        D2=V-(A(end,1)-star_data(D,3));
        m=1;
        [~,Da]=intersect(I,D);
        n=4908-Da;
        while (A(Da+m,1)-A(Da,1))<=D2
            m=m+1;
        end
        G1=[1:D+j,4908-i:4908];
        G=intersect(G1,I(Da-m:end,1));
        if size(G,1)>N
          G=randsample(G,N);
        end
        [p1,p2]=projection(D,G,star);
        [pre,epsion,perm]=callo(G,D,star);
        figure
        plot3(p1(:,1),p1(:,2),p1(:,3),'w.','MarkerSize',0.5)
        set(gca,'color',[0,0,0])
        title('��������ǵ���άʾ��ͼ') 
        figure
        plot(p2(:,1),p2(:,2),'w.','MarkerSize',0.5)
        hold on
        plot(0,0,'r*') %����Ŀ���
        set(gca,'color',[0,0,0])
        title('��������ǵ��άͶӰʾ��ͼ')
    end
    if (star_data(D,3)-A(1,1))<vi
       D3=40-(star_data(D,3)-A(1,1));
       m1=1;
       [~,Da]=intersect(I,D); 
       n1=Da;
       while (A(Da+m1,1)-star_data(D,3))<=D3
           m1=m1+1;
       end
        G1=[D+j:4908,1:i];
        G=intersect(G1,I(1:Da+m,1));
        if size(G,1)>N
          G=randsample(G',N);
        end
        [p1,p2]=projection(D,G,star);
        [pre,epsion,perm]=callo(G,D,star);
        figure
        plot3(p1(:,1),p1(:,2),p1(:,3),'w.','MarkerSize',0.5)
        set(gca,'color',[0,0,0])
        title('��������ǵ���άʾ��ͼ') 
        figure
        plot(p2(:,1),p2(:,2),'w.','MarkerSize',0.5)
        hold on
        plot(0,0,'r*') %����Ŀ���
        set(gca,'color',[0,0,0])
        title('��������ǵ��άͶӰʾ��ͼ')
    end
    m=1;
    n=1;
    [~,Da]=intersect(I,D);
    while (A(Da+m,1)-A(Da,1))<=vi
        m=m+1;
    end
    while (A(Da,1)-A(Da-n,1))<=vi
        n=n+1;
    end
%     [~,Da]=intersect(I,D);
    G1=[D+j:4908,1:i];
    G=intersect(G1,I(Da-m:Da+n,1));
    if size(G,1)>N
      G=randsample(G,N);
    end
    [p1,p2]=projection(D,G,star);
    [pre,epsion]=callo(G,D,star);
    figure
    plot3(p1(:,1),p1(:,2),p1(:,3),'w.','MarkerSize',0.5)
    set(gca,'color',[0,0,0])
    title('��������ǵ���άʾ��ͼ') 
    figure
    plot(p2(:,1),p2(:,2),'w.','MarkerSize',0.5)
    hold on 
    plot(0,0,'r*') %����Ŀ���
    set(gca,'color',[0,0,0])
    title('��������ǵ��άͶӰʾ��ͼ')
end
m=1;
n=1;
while (star_data(D+m,2)-star_data(D,2))<=vi
    m=m+1;
end
while (star_data(D,2)-star_data(D-n,2))<=vi
    n=n+1;
end
if (A(end,1)-star_data(D,3))<vi
    D2=V-(A(end,1)-star_data(D,3));
    m1=1;
    [~,Da]=intersect(I,D);
%     n1=4908-Da;
    while (A(Da-m1,1)-A(Da,1))<D2
         m1=m1+1;
    end
    G1=D-n:D+m;
    G=intersect(G1,I(Da-m1:end,1));
    if size(G,1)>N
       G=randsample(G,N);
    end
    [p1,p2]=projection(D,G,star);
    [pre,epsion,perm]=callo(G,D,star);
    figure
    plot3(p1(:,1),p1(:,2),p1(:,3),'w.','MarkerSize',0.5)
    set(gca,'color',[0,0,0])
    title('��������ǵ���άʾ��ͼ')
    figure
    plot(p2(:,1),p2(:,2),'w.','MarkerSize',0.5)
    hold on
    plot(0,0,'r*') %����Ŀ���
    set(gca,'color',[0,0,0])
    title('��������ǵ��άͶӰʾ��ͼ')
end
if (star_data(D,3)-A(1,1))<vi
   D3=V-(star_data(D,3)-A(1,1));
   m1=1;
   [~,Da]=intersect(I,D); 
   while (A(Da+m1,1)-star_data(D,3))<=D3
      m1=m1+1;
   end
   G1=D-n:D+m;
   G=intersect(G1,I(1:Da+m,1));
   if size(G,1)>N
      G=randsample(G,N);
   end
   [p1,p2]=projection(D,G,star);
   lo=zeros(100,1);
   [pre,epsion,perm]=callo(G,D,star);
   figure
   plot3(p1(:,1),p1(:,2),p1(:,3),'w.','MarkerSize',0.5)
   set(gca,'color',[0,0,0])
   title('��������ǵ���άʾ��ͼ')
   figure
   plot(p2(:,1),p2(:,2),'w.','MarkerSize',0.5)
   hold on 
   plot(0,0,'r*') %����Ŀ���
   set(gca,'color',[0,0,0])
   title('��������ǵ��άͶӰʾ��ͼ')
end
[~,Da]=intersect(I,D);
m1=1;
n1=1;
while (A(Da+m1,1)-star_data(D,3))<=vi
    m1=m1+1;
end
while (star_data(D,2)-A(D-n1))<=vi
    n1=n1+1;
end
G1=D-n:D+m;
G=intersect(G1,I((Da-n1:Da+m1),1));
if size(G,1)>N
  G=randsample(G,N);
end
[p1,p2]=projection(D,G,star);
[pre,epsion,perm,pre1,pre2]=callo(G,D,star);
figure
plot3(p1(:,1),p1(:,2),p1(:,3),'w.','MarkerSize',0.5)
set(gca,'color',[0,0,0])
title('��������ǵ���άʾ��ͼ')
figure
plot(p2(:,1),p2(:,2),'w.','MarkerSize',0.5)
hold on
plot(0,0,'r*') %����Ŀ���
set(gca,'color',[0,0,0])
title('��������ǵ��άͶӰʾ��ͼ')



    
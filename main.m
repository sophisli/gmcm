%此程序为星敏感器星图定位模型算法模拟测试
%我们假设星表中的数据都是误差极小的标准值
%首先令焦距f=1，即视天球为半径为1的球体
%然后取大小固定的像场，随机在天球上取任意一星点，作为目标星点
%该测试中假设所有测量误差均为随机误差，即需要给题目中要求的所有测试数据加上一定的随机扰动
clear all
load('附件2  简易星表.mat')

star=zeros(4908,3);%将星表中的数据从球坐标系转换为直角坐标系
star(:,1)=cos(star_data(:,3)*pi/180).*cos(star_data(:,2)*pi/180);
star(:,2)=cos(star_data(:,3)*pi/180).*sin(star_data(:,2)*pi/180);
star(:,3)=sin(star_data(:,3)*pi/180);
D=randsample(star_data(:,1),1);%从星表中随机选取一个星点作为目标点
[A,I]=sort(star_data(:,3));%由于星表赤纬的顺序是打乱的，要重排序，以便比对和筛选
vi=12;%设置目标天区
V=24;%设置天区搜索范围，单位为度
N=20;%设定选取星点的数目
%% 不同条件的目标随机筛选
if (star_data(end,2)-star_data(D,2))<vi %如果随机选择的星点位于星表赤经一列的末端，且与星表最后一颗星赤经少于vi
    D1=vi-(star_data(end,2)-star_data(D,2));
    i=1;
    j=1;
    while star_data(i,2)<=D1
        i=i+1;
    end
    while (star_data(D,2)-star_data(D-j,2))<=vi
        j=j+1;
    end
    if (A(end,1)-star_data(D,3))<vi %如果随机选择的星点位于星表赤纬一列的末端，且与星表最后一颗星赤纬少于vi
        D2=V-(A(end,1)-star_data(D,3));
        m=1;
        [~,Da]=intersect(I,D);%找到星表赤纬信息排序之后，目标星点的位置
        n=4908-Da;
        while (A(Da,1)-A(Da-m,1))<=D2
            m=m+1;
        end
        G1=[D+j:4908,1:i];
        G=intersect(G1,I(Da-m:end,1));%寻找同时满足赤纬范围和赤经范围的星点，即寻找二者交集
        if size(G,1)>N %随机抽取一定量的星点
          G=randsample(G,N);
        end
        [p1,p2]=projection(D,G,star); %调用projection
        [pre,epsion,perm]=callo(G,D,star);%调用callo
        figure
        plot3(p1(:,1),p1(:,2),p1(:,3),'w.','MarkerSize',0.5)%绘制选取星点在天区的示意图
        set(gca,'color',[0,0,0])
        title('随机抽样星点三维示意图')
        figure
        plot(p2(:,1),p2(:,2),'w.','MarkerSize',0.5)%绘制选取星点投影在二维平面的示意图
        hold on
        plot(0,0,'r*') %画出目标点
        set(gca,'color',[0,0,0])
        title('随机抽样星点二维投影示意图')
    end
    if (star_data(D,3)-A(1,1))<vi %如果随机选择的星点位于星表赤纬一列的开端，且与星表第一颗星赤纬少于vi
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
       title('随机抽样星点三维示意图') 
       figure
       plot(p2(:,1),p2(:,2),'w.','MarkerSize',0.5)
       hold on
       plot(0,0,'r*') %画出目标点
       set(gca,'color',[0,0,0])
       title('随机抽样星点二维投影示意图')
    end
    % 以下同上面步骤对所有随机抽取的情况判断，并执行相同的操作
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
    title('随机抽样星点三维示意图') 
    figure
    plot(p2(:,1),p2(:,2),'w.','MarkerSize',0.5)
    hold on
    plot(0,0,'r*') %画出目标点
    set(gca,'color',[0,0,0])
    title('随机抽样星点二维投影示意图')
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
        title('随机抽样星点三维示意图') 
        figure
        plot(p2(:,1),p2(:,2),'w.','MarkerSize',0.5)
        hold on
        plot(0,0,'r*') %画出目标点
        set(gca,'color',[0,0,0])
        title('随机抽样星点二维投影示意图')
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
        title('随机抽样星点三维示意图') 
        figure
        plot(p2(:,1),p2(:,2),'w.','MarkerSize',0.5)
        hold on
        plot(0,0,'r*') %画出目标点
        set(gca,'color',[0,0,0])
        title('随机抽样星点二维投影示意图')
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
    title('随机抽样星点三维示意图') 
    figure
    plot(p2(:,1),p2(:,2),'w.','MarkerSize',0.5)
    hold on 
    plot(0,0,'r*') %画出目标点
    set(gca,'color',[0,0,0])
    title('随机抽样星点二维投影示意图')
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
    title('随机抽样星点三维示意图')
    figure
    plot(p2(:,1),p2(:,2),'w.','MarkerSize',0.5)
    hold on
    plot(0,0,'r*') %画出目标点
    set(gca,'color',[0,0,0])
    title('随机抽样星点二维投影示意图')
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
   title('随机抽样星点三维示意图')
   figure
   plot(p2(:,1),p2(:,2),'w.','MarkerSize',0.5)
   hold on 
   plot(0,0,'r*') %画出目标点
   set(gca,'color',[0,0,0])
   title('随机抽样星点二维投影示意图')
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
title('随机抽样星点三维示意图')
figure
plot(p2(:,1),p2(:,2),'w.','MarkerSize',0.5)
hold on
plot(0,0,'r*') %画出目标点
set(gca,'color',[0,0,0])
title('随机抽样星点二维投影示意图')



    
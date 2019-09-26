%我们从特定视场的天区筛选了一定数量的星点
%从这些星点中做排列组合
%遍历组合寻找计算精度最高的位置
function [pre,epsion,Perm1,pre1,pre2]=callo(G,D,star)
load('附件2  简易星表.mat')
epsion=zeros(5,2);
Perm=nchoosek(G,3);%从特定天区中的星点中做组合
m=size(Perm,1);
Ep=zeros(m,2);
pre=zeros(m,2);
G(:,2)=star_data(G(:,1),2);
G(:,3)=star_data(G(:,1),3);

for n=1:m
    P=Perm(n,:);%循环遍历所有组合，找到合适位置的组合
    pre(n,:)=trirev(P,star,D);% 位置计算
    
% for n=1:m
%     Ex1=cos(star_data(Perm(n,1),2))*cos(star_data(Perm(n,1),3));
%     Ex2=cos(star_data(Perm(n,2),2))*cos(star_data(Perm(n,2),3));
%     Ex3=cos(star_data(Perm(n,3),2))*cos(star_data(Perm(n,3),3));
%     Ey1=sin(star_data(Perm(n,1),2))*cos(star_data(Perm(n,1),3));
%     Ey2=sin(star_data(Perm(n,2),2))*cos(star_data(Perm(n,2),3));
%     Ey3=sin(star_data(Perm(n,3),2))*cos(star_data(Perm(n,3),3));
%     Ez1=sin(star_data(Perm(n,1),3));
%     Ez2=sin(star_data(Perm(n,2),3));
%     Ez3=sin(star_data(Perm(n,3),3));
%     E=[Ex2-Ex1,Ey2-Ey1,Ez2-Ez1;Ex3-Ex1,Ey3-Ey1,Ez3-Ez1;Ex3-Ex2,Ey3-Ey2,Ez3-Ez2];
%     x=linsolve(E,[0;0;0]);
%     f=sqrt(x(1,1)^2+x(2,1)^2+x(3,1)^2);
%     A(n,1)=acos(x(1,1)/sqrt(x(1,1)^2+x(2,1)^2));
%     A(n,2)=asin(x(3,1)/f);
    Ep(n,1)=abs((pre(n,1)-star_data(D,2))/star_data(D,2));
    Ep(n,2)=abs((pre(n,2)-star_data(D,3))/star_data(D,3));  
end
[Et,It]=sort(Ep(:,1));%将相对误差值做排序，并输出index
[Et1,~]=sort(Ep(:,2));
pre1=star_data(It(1:5),2);
pre2=star_data(It(1:5),3);
epsion(:,1)=Et(1:5,1);%输出相对误差值
epsion(:,2)=Et1(1:5,1);
Perm1=Perm(It(1:5,1),:);%输出最佳的十种星点位置组合







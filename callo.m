%���Ǵ��ض��ӳ�������ɸѡ��һ���������ǵ�
%����Щ�ǵ������������
%�������Ѱ�Ҽ��㾫����ߵ�λ��
function [pre,epsion,Perm1,pre1,pre2]=callo(G,D,star)
load('����2  �����Ǳ�.mat')
epsion=zeros(5,2);
Perm=nchoosek(G,3);%���ض������е��ǵ��������
m=size(Perm,1);
Ep=zeros(m,2);
pre=zeros(m,2);
G(:,2)=star_data(G(:,1),2);
G(:,3)=star_data(G(:,1),3);

for n=1:m
    P=Perm(n,:);%ѭ������������ϣ��ҵ�����λ�õ����
    pre(n,:)=trirev(P,star,D);% λ�ü���
    
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
[Et,It]=sort(Ep(:,1));%��������ֵ�����򣬲����index
[Et1,~]=sort(Ep(:,2));
pre1=star_data(It(1:5),2);
pre2=star_data(It(1:5),3);
epsion(:,1)=Et(1:5,1);%���������ֵ
epsion(:,2)=Et1(1:5,1);
Perm1=Perm(It(1:5,1),:);%�����ѵ�ʮ���ǵ�λ�����







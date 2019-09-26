%计算投影到视场平面的三维和二维星点位置
function [p1,p2]=projection(D,G,star)
a1=zeros(size(G,1),1);
for n=1:size(G,1)
    a1(n,1)=sqrt(1-(2-(star(D,1)-star(G(n),1))^2-(star(D,2)-star(G(n),2))^2-...
        (star(D,3)-star(G(n),3))^2)/2)/((2-(star(D,1)-star(G(n),1))^2-(star(D,2)-star(G(n),2))^2-...
        (star(D,3)-star(G(n),3))^2)/2);
end

p1=zeros(size(G,1),3);
for n=1:size(G,1)
    p1(n,1)=sqrt(1+a1(n,1)^2)*star(G(n),1)-star(D,1);
    p1(n,2)=sqrt(1+a1(n,1)^2)*star(G(n),2)-star(D,2);
    p1(n,3)=sqrt(1+a1(n,1)^2)*star(G(n),3)-star(D,3);
end
p2=zeros(size(G,1),2);
tem=zeros(3,size(G,1));
%% 投影矩阵
% syms a b
% M=[cos(a)*cos(b),-sin(a),-sin(b)*cos(a);sin(a)*cos(b),cos(a),-sin(a)*sin(b);sin(a),0,cos(a)];

for n=1:size(G,1)
    if star(D,2)<0
       alpha=360-acosd(star(D,1)/sqrt(star(D,1)^2+star(D,2)^2));
       delta=asind(star(D,3));
    end
    if star(D,2)>0
       alpha=acosd(star(D,1)/sqrt(star(D,1)^2+star(D,2)^2));
       delta=asind(star(D,3));
    end
    Ma=[cosd(alpha),-sind(delta),0;sind(alpha),cosd(alpha),0;0,0,1];%赤经的旋转矩阵，即以z轴为旋转轴
    Mb=[cosd(delta),0,-sind(delta);0,1,0;sind(delta),0,cosd(delta)];%赤纬的旋转矩阵，以y轴为旋转轴
    M=Ma*Mb;
    tem(:,n)=M\star(G(n),:)';
    p2(n,:)=tem(2:end,n);
end

    
    



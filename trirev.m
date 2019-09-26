%截距三棱锥修正定位法
%为了方便讨论，这里我们假设f=1
%P为特定视场中的一种组合
%算法计算主要用到余弦定理
%算法的实现在于每次循环的判断和实时修正更新
function pre=trirev(P,star,D)
% P=G(1:3,1)';
a=zeros(3,1);
for n=1:3
    a(n,1)=sqrt((1/((2-(star(D,1)-star(P(n),1))^2-...
        (star(D,2)-star(P(n),2))^2-(star(D,3)-star(P(n),3))^2)/2))^2-1);
%     a(n,1)=sqrt(1-(2-(star(D,1)-star(P(n),1))^2-(star(D,2)-star(P(n),2))^2-...
%         (star(D,3)-star(P(n),3))^2)/2)/((2-(star(D,1)-star(P(n),1))^2-...
%         (star(D,2)-star(P(n),2))^2-(star(D,3)-star(P(n),3))^2)/2);
end
for n=1:size(a,1) %加入随机扰动
     if round(rand(1))==1
        a(n,1)=a(n,1)+1e-6*rand(1);
     end
     a(n,1)=a(n,1)-1e-6*rand(1);
end
pre=zeros(1,2);
l1=sqrt(a(1,1)^2+1);
l2=sqrt(a(2,1)^2+1);
l3=sqrt(a(3,1)^2+1);
Ev1=l1*star(P(1,1),:);
Ev2=l2*star(P(1,2),:);
Ev3=l3*star(P(1,3),:);
C12=Ev2-Ev1;
C23=Ev3-Ev2;
C31=Ev1-Ev3;
n=0;
%% 主动修正，循环试探过程
while n<=100 %做足够次数的循环，使得三棱锥成为截距三棱锥
if (sum(Ev1.^2)+sum(Ev2.^2))~=sum(C12.^2) %判断三棱锥不同侧面棱之间是否垂直
    Cos=((sum(Ev1.^2)+sum(C12.^2)-sum(Ev2.^2))/(2*sqrt(sum(Ev1.^2))*sqrt(sum(C12.^2))));%通过余弦定理计算棱1和底边夹角余弦值
    Cos1=((sum(Ev2.^2)+sum(C12.^2)-sum(Ev1.^2))/(2*sqrt(sum(Ev2.^2))*sqrt(sum(C12.^2))));%通过余弦定理计算棱2和底边夹角余弦值
    theta=acosd(Cos);
    if theta<90 %判断棱和相邻底边夹角是锐角还是钝角，锐角执行该步骤
        C12=sqrt(sum(Ev1.^2))*(C12/(sqrt(sum(C12.^2))))/Cos;
        Ev2=Ev1+C12; %更新棱的位置信息
        C23=Ev3-Ev2; %更新相应底边的信息
    end
    if theta>90 %钝角，执行该步骤
        C12=sqrt(sum(Ev2.^2))*(C12/(sqrt(sum(C12.^2))))/Cos1;
        Ev1=Ev2-C12;
        C31=Ev1-Ev3;
    end
%     y1=sqrt(sum(Ev1.^2))*sqrt(1-((sum(Ev1.^2)+sum(C12.^2)-sum(Ev2.^2))/(2*sqrt(sum(Ev1.^2))*sqrt(sum(C12.^2))))^2)...
%             /((sum(Ev1.^2)+sum(C12.^2)-sum(Ev2.^2))/(2*sqrt(sum(Ev1.^2))*sqrt(sum(C12.^2))));
%     C12=(sqrt(y1^2+sum(Ev1.^2)))*C12/sqrt(sum(C12.^2));

%分别执行三棱锥三个侧面的判断和修正，步骤同上
end
if (sum(Ev2.^2)+sum(Ev3.^2))~=sum(C23.^2)
    Cos=(sum(Ev2.^2)+sum(C23.^2)-sum(Ev3.^2))/(2*sqrt(sum(Ev2.^2))*sqrt(sum(C23.^2)));
    Cos1=(sum(Ev3.^2)+sum(C23.^2)-sum(Ev2.^2))/(2*sqrt(sum(Ev3.^2))*sqrt(sum(C23.^2)));
    theta=acosd(Cos);
    if theta<90
        C23=sqrt(sum(Ev2.^2))*(C23/(sqrt(sum(C23.^2))))/Cos;
        Ev3=Ev2+C23;
        C31=Ev1-Ev3;
    end
    if theta>90
        C23=sqrt(sum(Ev3.^2))*(C23/(sqrt(sum(C23.^2))))/Cos1;
        Ev2=Ev3-C23;
        C12=Ev2-Ev1;
    end
    
end
if (sum(Ev3.^2)+sum(Ev1.^2))~=sum(C31.^2)
    Cos=((sum(Ev3.^2)+sum(C31.^2)-sum(Ev1.^2))/(2*sqrt(sum(Ev3.^2))*sqrt(sum(C31.^2))));
    Cos1=((sum(Ev1.^2)+sum(C31.^2)-sum(Ev3.^2))/(2*sqrt(sum(Ev1.^2))*sqrt(sum(C31.^2))));
    theta=acosd(Cos);
    if theta<90
        C31=sqrt(sum(Ev3.^2))*(C31/(sqrt(sum(C31.^2))))/Cos;
        Ev1=Ev3+C31;
        C12=Ev2-Ev1;
    end
    if theta>90
        C31=sqrt(sum(Ev1.^2))*(C31/(sqrt(sum(C31.^2))))/Cos1;
        Ev3=Ev1-C31;
        C23=Ev3-Ev2;
    end
end
%通过修正后的三条相互垂直的侧棱表示三棱锥的垂线，即目标星的位置
OD(1,1)=1/sqrt(sum(Ev1.^2));
OD(2,1)=1/sqrt(sum(Ev2.^2));
OD(3,1)=1/sqrt(sum(Ev3.^2));
OD=OD/sqrt(sum(OD.^2));
%需要提到的是这里转换计算的结果是基于当前三侧棱为坐标轴的球坐标经纬，需要做一个旋转变换才能求得在天球坐标系下的目标点位置
%即我们设Ev1为x轴，Ev2方向为y轴。Ev3方向为z轴
%此时我们可以很容易通过矢量叠加找到该坐标系和天球坐标系的转换关系
ex=Ev1/sqrt(sum(Ev1.^2));
ey=Ev2/sqrt(sum(Ev2.^2));
ez=Ev3/sqrt(sum(Ev3.^2));
M=[ex;ey;ez]';%三条侧棱在天球坐标系下的表示矩阵
OD=M*OD; %目标星点在天球坐标系下的直角坐标表示
%将直角坐标转换为球坐标
if OD(2,1)<0
   pre(1,1)=360-acosd(OD(1,1)/sqrt(OD(1,1)^2+OD(2,1)^2));
   pre(1,2)=asind(OD(3,1));
end
if OD(2,1)>0
   pre(1,1)=acosd(OD(1,1)/sqrt(OD(1,1)^2+OD(2,1)^2));
   pre(1,2)=asind(OD(3,1));
end
n=n+1;
end


    

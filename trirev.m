%�ؾ�����׶������λ��
%Ϊ�˷������ۣ��������Ǽ���f=1
%PΪ�ض��ӳ��е�һ�����
%�㷨������Ҫ�õ����Ҷ���
%�㷨��ʵ������ÿ��ѭ�����жϺ�ʵʱ��������
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
for n=1:size(a,1) %��������Ŷ�
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
%% ����������ѭ����̽����
while n<=100 %���㹻������ѭ����ʹ������׶��Ϊ�ؾ�����׶
if (sum(Ev1.^2)+sum(Ev2.^2))~=sum(C12.^2) %�ж�����׶��ͬ������֮���Ƿ�ֱ
    Cos=((sum(Ev1.^2)+sum(C12.^2)-sum(Ev2.^2))/(2*sqrt(sum(Ev1.^2))*sqrt(sum(C12.^2))));%ͨ�����Ҷ��������1�͵ױ߼н�����ֵ
    Cos1=((sum(Ev2.^2)+sum(C12.^2)-sum(Ev1.^2))/(2*sqrt(sum(Ev2.^2))*sqrt(sum(C12.^2))));%ͨ�����Ҷ��������2�͵ױ߼н�����ֵ
    theta=acosd(Cos);
    if theta<90 %�ж�������ڵױ߼н�����ǻ��Ƕ۽ǣ����ִ�иò���
        C12=sqrt(sum(Ev1.^2))*(C12/(sqrt(sum(C12.^2))))/Cos;
        Ev2=Ev1+C12; %�������λ����Ϣ
        C23=Ev3-Ev2; %������Ӧ�ױߵ���Ϣ
    end
    if theta>90 %�۽ǣ�ִ�иò���
        C12=sqrt(sum(Ev2.^2))*(C12/(sqrt(sum(C12.^2))))/Cos1;
        Ev1=Ev2-C12;
        C31=Ev1-Ev3;
    end
%     y1=sqrt(sum(Ev1.^2))*sqrt(1-((sum(Ev1.^2)+sum(C12.^2)-sum(Ev2.^2))/(2*sqrt(sum(Ev1.^2))*sqrt(sum(C12.^2))))^2)...
%             /((sum(Ev1.^2)+sum(C12.^2)-sum(Ev2.^2))/(2*sqrt(sum(Ev1.^2))*sqrt(sum(C12.^2))));
%     C12=(sqrt(y1^2+sum(Ev1.^2)))*C12/sqrt(sum(C12.^2));

%�ֱ�ִ������׶����������жϺ�����������ͬ��
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
%ͨ��������������໥��ֱ�Ĳ����ʾ����׶�Ĵ��ߣ���Ŀ���ǵ�λ��
OD(1,1)=1/sqrt(sum(Ev1.^2));
OD(2,1)=1/sqrt(sum(Ev2.^2));
OD(3,1)=1/sqrt(sum(Ev3.^2));
OD=OD/sqrt(sum(OD.^2));
%��Ҫ�ᵽ��������ת������Ľ���ǻ��ڵ�ǰ������Ϊ������������꾭γ����Ҫ��һ����ת�任�����������������ϵ�µ�Ŀ���λ��
%��������Ev1Ϊx�ᣬEv2����Ϊy�ᡣEv3����Ϊz��
%��ʱ���ǿ��Ժ�����ͨ��ʸ�������ҵ�������ϵ����������ϵ��ת����ϵ
ex=Ev1/sqrt(sum(Ev1.^2));
ey=Ev2/sqrt(sum(Ev2.^2));
ez=Ev3/sqrt(sum(Ev3.^2));
M=[ex;ey;ez]';%������������������ϵ�µı�ʾ����
OD=M*OD; %Ŀ���ǵ�����������ϵ�µ�ֱ�������ʾ
%��ֱ������ת��Ϊ������
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


    

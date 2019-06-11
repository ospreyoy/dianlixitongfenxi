%*********************************************
%Name:H.Y. Wang
%No:3016203104
%Data:2019.06.11
%*********************************************
%ע��������������ݾ��洢����OUT.xlsx���ļ�֮��
%*********************************************
%�����γɽڵ㵼�ɾ���Y�������ļ��б�
A=[1 4 0 0.0576 0; 
   2 7 0 0.0625 0;
   3 9 0 0.0586 0;
   4 5 0.010 0.085 0.088;
   4 6 0.017 0.092 0.079; 
   5 7 0.032 0.161 0.153; 
   6 9 0.039 0.170 0.179; 
   7 8 0.0085 0.0720 0.0745; 
   8 9 0.0119 0.1008 0.1045]; 
%�ڵ㵼�ɾ���Y
Y=zeros(9,9); 
%����ڵ㵼�ɾ���Y
for i=1:9
    m=A(i,1);
    n=A(i,2);
    Y(m,m)=Y(m,m)+1/(A(i,3)+1i*A(i,4))+1i*A(i,5); 
    Y(n,n)=Y(n,n)+1/(A(i,3)+1i*A(i,4))+1i*A(i,5); 
    Y(m,n)=Y(m,n)-1/(A(i,3)+1i*A(i,4)); 
    Y(n,m)=Y(n,m)-1/(A(i,3)+1i*A(i,4));
end

U=[1.04,1.025,1.025,1,1,1,1,1,1]; 
P=[1,1.63,0.85,0,-1.25,-0.9,0,-1,0]; 
Q=[0,0,0,0,-0.5,-0.3,0,-0.35,0]; 
dP=[0,0,0,0,0,0,0,0,0]; 
dQ=[0,0,0,0,0,0,0,0,0]; 
a=[0,0,0,0,0,0,0,0,0]; 
%��ѹ����U���й���������P���޹���������Q���������a���й�����ֵdP���޹�����ֵdQ

G=real(Y);
B=imag(Y); 

Pt=zeros(1,9); 
Qt=zeros(1,9); 
%�����Ÿ��Ⱦ����Ӿ���H��N��M��L
H=zeros(8,8); 
N=zeros(8,6); 
M=zeros(6,8); 
L=zeros(6,6);

Ai=zeros(1,9); 
Bi=zeros(1,9); 
JJ=zeros(14,14); 
I=zeros(1,9); 
%����ϳ��Ÿ��Ⱦ���JJ
%�����������k
k=0;
jingdu=1;
%ţ��-����ѷ����ѭ��
while jingdu>0.0000001
	for m=2:9 
    	for n=1:9 
        	pt(n)=U(m)*U(n)*(G(m,n)*cos(a(m)-a(n))+B(m,n)*sin(a(m)-a(n))); 
        end 
    	dP(m)=P(m)-sum(pt); 
     end 
     for m=4:9 
        for n=1:9
            qt(n)=U(m)*U(n)*(G(m,n)*sin(a(m)-a(n))-B(m,n)*cos(a(m)-a(n))); 
        end
        dQ(m)=Q(m)-sum(qt);
     end
     
     for m=1:8
         for n=1:8
             if m==n
             else
                 H(m,n)=-U(m+1)*U(n+1)*(G(m+1,n+1)*sin(a(m+1)-a(n+1))-B(m+1,n+1)*cos(a(m+1)-a(n+1))); 
             end
         end
     end
     for m=1:8
         for n=1:6
             if m==n+2
             else
                 N(m,n)=-U(m+1)*U(n+3)*(G(m+1,n+3)*cos(a(m+1)-a(n+3))+B(m+1,n+3)*sin(a(m+1)-a(n+3))); 
             end
         end
     end
     for m=1:6
         for n=1:8
             if m+2==n
             else
                M(m,n)=U(m+3)*U(n+1)*(G(m+3,n+1)*cos(a(m+3)-a(n+1))+B(m+3,n+1)*sin(a(m+3)-a(n+1)));  
             end
         end
     end
     for m=1:6
         for n=1:6
             if m==n
             else
                L(m,n)=-U(m+3)*U(n+3)*(G(m+3,n+3)*sin(a(m+3)-a(n+3))-B(m+3,n+3)*cos(a(m+3)-a(n+3))); 
             end
         end
     end
     for m=1:8
         for n=1:9
             Ai(n)=U(n)*(G(m+1,n)*sin(a(m+1)-a(n))-B(m+1,n)*cos(a(m+1)-a(n))); 
         end
         H(m,m)=U(m+1)*U(m+1)*B(m+1,m+1)+U(m+1)*sum(Ai); 
     end
     for m=1:6 
         for n=1:9 
             Bi(n)=U(n)*(G(m+1,n)*cos(a(m+1)-a(n))+B(m+1,n)*sin(a(m+1)-a(n))); 
         end
         N(m+2,m)=-U(m+1)*U(m+1)*G(m+1,m+1)-U(m+1)*sum(Bi); 
     end
     for m=1:6 
         for n=1:9 
             Bi(n)=U(n)*(G(m+3,n)*cos(a(m+3)-a(n))+B(m+3,n)*sin(a(m+3)-a(n))); 
         end
         M(m,m+2)=U(m+3)*U(m+3)*G(m+3,m+3)-U(m+3)*sum(Bi); 
     end
     for m=1:6 
         for n=1:9
             Ai(n)=U(n)*(G(m+3,n)*sin(a(m+3)-a(n))-B(m+3,n)*cos(a(m+3)-a(n))); 
         end
         L(m,m)=U(m+3)*U(m+3)*B(m+3,m+3)-U(m+3)*sum(Ai); 
     end
     for m=1:8
         for n=1:8
             JJ(m,n)=H(m,n); 
         end
     end
     for m=1:8 
         for n=1:6
             JJ(m,n+8)=N(m,n);
         end
     end
     for m=1:6
         for n=1:8
             JJ(m+8,n)=M(m,n);
         end
     end
     for m=1:6
         for n=1:6
             JJ(m+8,n+8)=L(m,n);
         end
     end
     for m=1:8
         PQ(m)=dP(m+1);
     end
     for m=9:14
         PQ(m)=dQ(m-5);
     end
     
     dUa=-inv(JJ)*PQ';
     jingdu=max(abs(dUa));
     for m=1:8
         a(m+1)=a(m+1)+dUa(m);
     end
     for m=9:14
         U(m-5)=U(m-5)+dUa(m);
     end
     k=k+1;
end

u=U.*cos(a)+1i*U.*sin(a); 
i=Y*u.'; 
I=conj(i); 
pq=u.*I.'; 
P=real(pq); 
Q=imag(pq);

%�ڵ㵼�ɾ���Y���
for qq=1:9
   for ww=1:9
       Yzifu(qq,ww)={num2str(Y(qq,ww))};
   end
end
xlswrite('OUT.xlsx',Yzifu,'�ڵ㵼�ɾ���Y');

%��ѹ��ֵU���
for qq=1:9
    Uzifu(qq)={num2str(U(qq))};
end
xlswrite('OUT.xlsx',Uzifu,'��ѹ��ֵ');

%��ѹ���a���
for qq=1:9
    azifu(qq)={num2str(a(qq)*360/(2*pi))};
end
xlswrite('OUT.xlsx',azifu,'��ѹ���');

%�ڵ��й�����P���
for qq=1:9
    Pzifu(qq)={num2str(P(qq))};
end
xlswrite('OUT.xlsx',Pzifu,'�ڵ��й�����');

%�ڵ��޹�����Q���
for qq=1:9
    Qzifu(qq)={num2str(Q(qq))};
end
xlswrite('OUT.xlsx',Qzifu,'�ڵ��޹�����');

%*******************�ָ���********************
%����������븺�ɽڵ��Ե���
fadianji=[1 0.3 1.137; 
   2 0.3 1.211; 
   3 0.3 1.043; 
   4 0 0; 
   5 0 0; 
   6 0 0; 
   7 0 0;  
   8 0 0; 
   9 0 0];
for m=1:9 
    if abs(P(m))>0.000001 && abs(Q(m))>0.000001
        if P(m)>0
            Y(m,m)=Y(m,m)-1i*(1/fadianji(m,2));
        else
            Y(m,m)=Y(m,m)-conj(pq(m))/U(m)/U(m);
        end
    end
end
%������Ľڵ��Ե������
for qq=1:9
    Yxiuzhengzifu(qq)={num2str(Y(qq,qq))};
end
xlswrite('OUT.xlsx',Yxiuzhengzifu,'�����ڵ��Ե���');



%*******************�ָ���********************
Z=inv(Y);
Z4lie=Z(:,4);
%����迹����ĵ�����
for qq=1:9
    Z4liezifu(qq)={num2str(Z4lie(qq))};
end
xlswrite('OUT.xlsx',Z4liezifu.','�迹���������');

%*******************�ָ���********************
%����ڵ��ѹ���·����
If=u(4)/Z(4,4);
for m=1:9
    uu(m)=u(m)-Z4lie(m)*If;
end
If1=1/Z(4,4);
for m=1:9
    uu1(m)=1-Z4lie(m)/Z(4,4);
end
%*******************�ָ���********************
for m=1:9 
    hang=A(m,:);
    y=zeros(2,2);
    if hang(1)>3 && hang (2)>3
        y(1,1)=y(1,1)+1/(A(m,3)+1i*A(m,4))+1i*A(m,5); 
        y(1,2)=y(1,2)-1/(A(m,3)+1i*A(m,4)); 
        y(2,1)=y(2,1)-1/(A(m,3)+1i*A(m,4)); 
        y(2,2)=y(2,2)+1/(A(m,3)+1i*A(m,4))+1i*A(m,5); 
        ii(m,1)=y(1,1)*uu(hang(1))+y(1,2)*uu(hang(2)); 
        ii(m,2)=y(2,1)*uu(hang(1))+y(2,2)*uu(hang(2));
    else
        y(1,1)=y(1,1)+1/(A(m,3)+1i*A(m,4)); 
        y(1,2)=y(1,2)-1/(A(m,3)+1i*A(m,4)); 
        y(2,1)=y(2,1)-1/(A(m,3)+1i*A(m,4)); 
        y(2,2)=y(2,2)+1/(A(m,3)+1i*A(m,4)); 
        ii(m,1)=(uu(hang(2))-uu(hang(1)))*y(1,2); 
        ii(m,2)=(uu(hang(1))-uu(hang(2)))*y(2,1);
    end
end


a=zeros(9,6);
for m=1:9
    a(m,1)=A(m,1); 
    a(m,2)=A(m,2); 
    a(m,3)=real(ii(m,1)); 
    a(m,4)=imag(ii(m,1)); 
    a(m,5)=real(ii(m,2)); 
    a(m,6)=imag(ii(m,2)); 
end

b=zeros(9,6);
for m=1:9 
    b(m,1)=A(m,1); 
    b(m,2)=A(m,2); 
    b(m,3)=abs(ii(m,1)); 
    b(m,4)=angle(ii(m,1))*360/(2*pi); 
    b(m,5)=abs(ii(m,2)); 
    b(m,6)=angle(ii(m,2))*360/(2*pi); 
end
for n=1:3
    for m=1:5
        s=a(m,:); 
        a(m,:)=a(m+1,:); 
        a(m+1,:)=s; 
        m=m+1; 
    end
    n=n+1;
end
for n=1:3
    for m=1:5
        s=a(m+3,:); 
        a(m+3,:)=a(m+4,:); 
        a(m+4,:)=s; 
        m=m+1;
    end
    n=n+1;
end
for n=1:3 
    for m=1:5 
        s=b(m,:); 
        b(m,:)=b(m+1,:); 
        b(m+1,:)=s; 
        m=m+1; 
    end 
    n=n+1; 
end
for n=1:3 
    for m=1:5 
        s=b(m+3,:); 
        b(m+3,:)=b(m+4,:); 
        b(m+4,:)=s; 
        m=m+1; 
    end 
    n=n+1; 
end

%��ȷ�������
xlswrite('OUT.xlsx',abs(If),'��·������ֵ����ȷ��');
xlswrite('OUT.xlsx',angle(If)*360/(2*pi),'��·������ǣ���ȷ��');
for qq=1:9
    Ujingque(qq)={num2str(abs(uu(qq)))};
end
xlswrite('OUT.xlsx',Ujingque,'�ڵ��ѹ��ֵ����ȷ��');

%���Ƽ������
xlswrite('OUT.xlsx',abs(If1),'��·������ֵ�����ƣ�');
xlswrite('OUT.xlsx',angle(If1)*360/(2*pi),'��·������ǣ����ƣ�');
for qq=1:9
    Ujinsi(qq)={num2str(abs(uu1(qq)))};
end
xlswrite('OUT.xlsx',Ujinsi,'�ڵ��ѹ��ֵ�����ƣ�');

%������
for qq=1:9
    Uwuchazifu(qq)={num2str(abs(uu1(qq))-abs(uu(qq)))};
end
xlswrite('OUT.xlsx',Uwuchazifu,'�ڵ��ѹ��ֵ��ֵ������ֵ-��ȷֵ��');
for qq=1:9
    Awuchazifu(qq)={num2str((angle(uu1(qq))-angle(uu(qq)))*360/(2*pi))};
end
Awuchazifu{1,4}=[0.0000];
xlswrite('OUT.xlsx',Awuchazifu,'�ڵ��ѹ��ǲ�ֵ������ֵ-��ȷֵ��');
xlswrite('OUT.xlsx',(abs(If1)-abs(If)),'��·������ֵ��ֵ������ֵ-��ȷֵ��');
xlswrite('OUT.xlsx',((angle(If1)-angle(If))*360/(2*pi)),'��·������ǲ�ֵ������ֵ-��ȷֵ��');
Iduanlu=cell(10,4);
Iduanlu(1,1)={'�ڵ�i'};
Iduanlu(1,2)={'�ڵ�j'};
Iduanlu(1,3)={'Iij'};
Iduanlu(1,4)={'Iji'};
for qq=2:10
    Iduanlu(qq,1)={num2str(a(qq-1,1))};
    Iduanlu(qq,2)={num2str(a(qq-1,2))};
    Iduanlu(qq,3)={num2str((a(qq-1,3))+1i*(a(qq-1,4)))};
    Iduanlu(qq,4)={num2str((a(qq-1,5))+1i*(a(qq-1,6)))};
end
xlswrite('OUT.xlsx',Iduanlu,'��ȷ�����µ�֧·����');
%*****************************************
%END
%*****************************************
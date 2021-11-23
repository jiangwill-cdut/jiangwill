clear
load XieModel.mat;
model21=model2(1:2:end,1:2:1254);
model31=model21(15:414,15:614);
[m,n]=size(model31);
a=[6057.669009	3331.001632	0.24
    5622.715772	2897.542884	0.28
    6172.839506	3113.131187	0.32
    ]
for i=1:m
    for j=1:n
        if model31(i,j)==1.5
            vp(i,j)=a(1,1);
            vs(i,j)=a(1,2);
            den(i,j)=a(1,3);
        end
        if model31(i,j)==1
            vp(i,j)=a(2,1);
            vs(i,j)=a(2,2);
            den(i,j)=a(2,3);
        end
        if model31(i,j)==0
            vp(i,j)=a(3,1);
            vs(i,j)=a(3,2);
            den(i,j)=a(3,3);
        end
    end
end

%正演
angle=45;
pre=[];

k=500;
model(:,1)=vp(:,k);
model(:,2)=vs(:,k);
model(:,3)=den(:,k);
for j=1:m-1
    for i=1:angle
        Rp(j,i)=avopp(model(j,1),model(j,2),model(j,3),model(j+1,1),model(j+1,2),model(j+1,3),i,3);
        %Rs(j,i)=avops(model(j,1),model(j,2),model(j,3),model(j+1,1),model(j+1,2),model(j+1,3),i,3);
    end
end

%Rp=abs(Rp);
Rp(find(isnan(Rp)==1))=0;
%Rp=abs(Rp);

%set wavelat
fm=30;
dt=0.001;
lengthwalvelate=100;
i=-lengthwalvelate/2:1:lengthwalvelate/2;
t=i*dt;
w=exp(-(pi*fm*t).^2).*(1-2*(pi*fm*t).^2);

%sesimic
[m1,n1]=size(Rp);
for i=1:n1
    Pseismic1(:,i)=conv(Rp(:,i),w);
end
Pseismic=Pseismic1(lengthwalvelate/2:end-lengthwalvelate/2,:);

P=Pseismic;


[n1,n2]=size(P);
randn('seed',20200309)
a=randn(1,n2);
a= floor(0.001*a/max(a)/dt);

y=0.026*wgn(n1,n2,0.01);

Pnoise=dither(P,a)+y;

[In1,in2,in3,in4]=BEMD(Pnoise);
figure(3)
subplot(2,2,1);wigb(In1);subplot(2,2,2);wigb(in2);
subplot(2,2,3);wigb(in3);subplot(2,2,4);wigb(in4);
figure(31)
wigb([In1;in2;in3;in4])
wigb([In1(200:250,:);in4])

figure(1);
subplot(1,4,1);wigb(P);subplot(1,4,3);wigb(Pnoise)
subplot(1,4,2);plot(StackWave(P));view(90,90);
subplot(1,4,4);plot(StackWave(Pnoise));view(90,90)
suptitle('含噪对比图')

[I1,i2,i3,i4]=BEMD(P);
[In1,in2,in3,in4]=BEMD(Pnoise);

figure(2)
subplot(2,2,1);wigb(I1);subplot(2,2,2);wigb(i2);
subplot(2,2,3);wigb(i3);subplot(2,2,4);wigb(i4);
figure(3)
subplot(2,2,1);wigb(In1);subplot(2,2,2);wigb(in2);
subplot(2,2,3);wigb(in3);subplot(2,2,4);wigb(in4);

m1=max(max([In1,in2,in3,in4]))
m2=min(min([In1,in2,in3,in4]))
m11(1:45)=m1;
m22(1:45)=m2;
wigb([m22;In1;m11])
wigb([m22;in2;m11])
wigb([m22;in3;m11])
wigb([m22;in4;m11])


figure(4)
subplot(1,4,1);plot(StackWave(P));view(90,90);subplot(1,4,2);plot(StackWave(I1));view(90,90);
subplot(1,4,3);plot(StackWave(Pnoise));view(90,90);subplot(1,4,4);plot(StackWave(In1));view(90,90);

Ptk=Teager2D(I1);Ptk0=Teager2D(P);


figure(5)
subplot(1,5,1);hht(P(:,1),500,'FrequencyLimits',[0 50]);title("Ori");view(90,90)
subplot(1,5,2);hht(I1(:,1),500,'FrequencyLimits',[0 50]);title('BEMD2');view(90,90)
subplot(1,5,3);hht(Ptk(:,1),500,'FrequencyLimits',[0 50]);title('BEMDtk');view(90,90)
subplot(1,5,4);hht(Ptk0(:,1),500,'FrequencyLimits',[0 50]);title('tk');view(90,90)
subplot(1,5,5);plot(StackWave(P));view(90,90)

figure(6)
subplot(1,5,1);hht(P,500,'FrequencyLimits',[0 50]);title("Ori");view(90,90)
subplot(1,5,2);hht(I1,500,'FrequencyLimits',[0 50]);title('BEMD2');view(90,90)
subplot(1,5,3);hht(Ptk,500,'FrequencyLimits',[0 50]);title('BEMDtk');view(90,90)
subplot(1,5,4);hht(Ptk0,500,'FrequencyLimits',[0 50]);title('tk');view(90,90)
subplot(1,5,5);plot(StackWave(P));view(90,90)


Ptkn=Teager2D(in2);Ptk0n=Teager2D(Pnoise);

for i=1:45
    Ptkn2(:,i)=smoothdata(in2(:,i),'gaussian',5);
end
Ptkn21=Teager2D(Ptkn2);
figure(2)
wigb(Ptkn21)
wigb(Ptkn)
wigb(Ptk0n)

figure(7)
subplot(1,5,1);hht(Pnoise(:,1),500,'FrequencyLimits',[0 50]);title("Ori");view(90,90)
subplot(1,5,2);hht(In1(:,1),500,'FrequencyLimits',[0 50]);title('BEMD2');view(90,90)
subplot(1,5,3);hht(Ptkn(:,1),500,'FrequencyLimits',[0 50]);title('BEMDtk');view(90,90)
subplot(1,5,4);hht(Ptk0n(:,1),500,'FrequencyLimits',[0 50]);title('tk');view(90,90)
subplot(1,5,5);plot(StackWave(Pnoise));view(90,90)

figure(8)
subplot(1,5,1);hht(Pnoise,500,'FrequencyLimits',[0 100]);title("Ori");view(90,90)
subplot(1,5,2);hht(In1,500,'FrequencyLimits',[0 100]);title('BEMD2');view(90,90)
subplot(1,5,3);hht(Ptkn,500,'FrequencyLimits',[0 100]);title('BEMDtk');view(90,90)
subplot(1,5,4);hht(Ptk0n,500,'FrequencyLimits',[0 100]);title('tk');view(90,90)
subplot(1,5,5);plot(Rp(:,1));view(90,90)

figure(8)
subplot(1,5,1);hht(P,500,'FrequencyLimits',[0 100]);title("Ori");view(90,90)
subplot(1,5,2);hht(Pnoise,500,'FrequencyLimits',[0 100]);title('BEMD2');view(90,90)
subplot(1,5,3);hht(Ptkn,500,'FrequencyLimits',[0 100]);title('BEMDtk');view(90,90)
subplot(1,5,4);hht(Ptk0n,500,'FrequencyLimits',[0 100]);title('tk');view(90,90)
subplot(1,5,5);plot(Rp(:,1));view(90,90)


A=hht(Ptk0n,500,'FrequencyLimits',[0 100]);
A1=full(A);
A2=sum(A1);
A3 =smoothdata(A2,'gaussian',10);
figure(13)
subplot(1,2,1);plot(Rp(:,1));view(90,90)
subplot(1,2,2);plot(A3);view(90,90)
findpeaks(A3);view(90,90)

figure(9)
subplot(1,9,1);plot(Rp(:,1));view(90,90)
subplot(1,9,2);plot(StackWave(P));title("Ori");view(90,90)
subplot(1,9,3);hht(StackWave(P),500,'FrequencyLimits',[0 50]);title("Ori");view(90,90)
subplot(1,9,4);plot(StackWave(I1));title('BEMD');view(90,90)
subplot(1,9,5);hht(StackWave(I1),500,'FrequencyLimits',[0 50]);title('BEMD');view(90,90)
subplot(1,9,6);plot(StackWave(Ptk));title('BEMDtk');view(90,90)
subplot(1,9,7);hht(StackWave(Ptk),500,'FrequencyLimits',[0 50]);title('BEMDtk');view(90,90)
subplot(1,9,8);plot(StackWave(Ptk0));title('tk');view(90,90)
subplot(1,9,9);hht(StackWave(Ptk0),500,'FrequencyLimits',[0 50]);title('tk');view(90,90)

figure(10)
subplot(1,9,1);plot(Rp(:,1));view(90,90)
subplot(1,9,2);plot(StackWave(Pnoise));title("Ori");view(90,90)
subplot(1,9,3);hht(StackWave(Pnoise),500,'FrequencyLimits',[0 50]);title("Ori");view(90,90)
subplot(1,9,4);plot(StackWave(In1));title('BEMD2');view(90,90)
subplot(1,9,5);hht(StackWave(In1),500,'FrequencyLimits',[0 50]);title('BEMD2');view(90,90)
subplot(1,9,6);plot(StackWave(Ptkn));title('BEMDtk');view(90,90)
subplot(1,9,7);hht(StackWave(Ptkn),500,'FrequencyLimits',[0 50]);title('BEMDtk');view(90,90)
subplot(1,9,8);plot(StackWave(Ptk0n));title('tk');view(90,90)
subplot(1,9,9);hht(StackWave(Ptk0n),500,'FrequencyLimits',[0 50]);title('tk');view(90,90)

Pstack=StackWave(P);PstackT=Teager(Pstack',1);

figure(11)
subplot(1,5,1);plot(Rp(:,1));view(90,90)
subplot(1,5,2);plot(Pstack);title("Ori");view(90,90)
subplot(1,5,3);hht(Pstack,500,'FrequencyLimits',[0 50]);title("Ori");view(90,90)
subplot(1,5,4);plot(PstackT);title('BEMD2');view(90,90)
subplot(1,5,5);hht(PstackT,500,'FrequencyLimits',[0 50]);title('BEMD2');view(90,90)

Pstackn=StackWave(Pnoise);PstackTn=Teager(Pstackn',1);

figure(12)
subplot(1,5,1);plot(Rp(:,1));view(90,90)
subplot(1,5,2);plot(Pstackn);title("Ori");view(90,90)
subplot(1,5,3);hht(Pstackn,500,'FrequencyLimits',[0 50]);title("Ori");view(90,90)
subplot(1,5,4);plot(PstackTn);title('BEMD2');view(90,90)
subplot(1,5,5);hht(PstackTn,500,'FrequencyLimits',[0 50]);title('BEMD2');view(90,90)



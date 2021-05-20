clear
clc
close all

N=12;%номер варианта 

%data
Lx=5080; %meter
Ly=4000; %meter
Lz1=20+2*N; %meter, толщина I пласта
Lz2=45+5*N;%meter, толщина II пласта

%координаты скважины 
%H1 - разность давлений в метрах

%Коэффициенты передачи  скорректированной модели  заданы в виде:
K1=0.02+0.01*N;
K2=0.002+0.001*N;

%object descroption
%ыерхний пласт
kx1=148.4/86400; % выбором kx можно обеспечить скорость 
ky1=24.64/86400; %  
kz1=14.64/86400;
eta1=0.000128;

%нижний пласт %размерность
kx2=148.4/86400;
ky2=24.64/86400;
kz2=14.64/86400;
eta2=0.00012;
b1=0.00128/86400;%коэффициент перетекания (перевести в секунды)

%E1,E2,E3;%усилительные звенья

%значения граничных условий из условий Дарси
tau1=2; % период вывода графика в сек * tau1
tau2=3600; %период вывода графика в сек (=час)
tau=tau2*tau1; %время моделирования в секундах
dtau=1; %шаг 
%два цикла сделны для вывода графиков =

Q=100; %m^3/сек дебет ()  входное воздействие    в секундах
Q=Q/86400;
Kf=10; % 5 в начальный момент задаём. 
%по заданному понижению меняем, т.е. выходу .
V=Kf*Q;
%заданное понижение

%170.75
%кол-во точки дискртизации %Nx,Ny,Nz по 15
Nx=10;
Ny=10;
Nz=9; %кол-во точек на каждый пласт
%добывающая скважина в нижнем пласте 

%координаты добывающей скважины
Xs1d=(2000+100+5*N)/Lx;%координаты скважины в долях от Lx
Ys1d=(2000)/Ly;%координаты скважины в долях от Ly
Xs1=ceil(Nx*Xs1d);
Ys1=ceil(Ny*Ys1d);
%кордината измерительного устройства относительно  нижнего валанжина
Zs2=ceil(5);
%координаты контрольной скважины
Xs2d=(2000)/Lx;%координаты скважины в долях от Lx
Ys2d=(2000)/Ly;%координаты скважины в долях от Ly
Xs2=ceil(Nx*Xs2d);
Ys2=ceil(Ny*Ys2d);
Zs2k=ceil(5);

%шаг
dx=Lx/(Nx-1); %шаг
dy=Ly/(Ny-1);
%dz=(Lz1+Lz2)/(Nz-1);%одинаковый шаг для каждого пласта
dz1=Lz1/(Nz-1);%равны
dz2=Lz2/(Nz-1);

%% инициализация 
for x_i=1:Nx
    for y_i=1:Ny

% обнуление массивов
for z_i=1:Nz
H1(x_i,y_i,z_i)=0;
H2(x_i,y_i,z_i)=0;
end
        
vF=50;
%начальное состояние пласта
for z_i=1:Nz
H1(x_i,y_i,z_i)=193-vF*(x_i-1)/(Nx-1);
end


for z_i=1:Nz
H2(x_i,y_i,z_i)=193-vF*(x_i-1)/(Nx-1);
TK(x_i,y_i,z_i)=H2(x_i,y_i,z_i);%начальное состояние
end

    end
end

tau=tau2*tau1; %время моделирования в секундах
dtau=tau/tau; %шаг 
%два цикла сделны для вывода графиков 

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');%шрифт
fig=figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
%h=plot(0,0);
sub1=subplot(2,1,1);
cla(sub1);
h = animatedline(sub1);
ax1=gca;
title("процессы в скважинах");
ylabel( [{"рассогласования уровня"};{" в точке измерения"};{" относительно начального"};{" состояния, м"}]);
xlabel("                       возмущющая скважина,   время, час");
grid on
drawnow
%контрольная 
sub2=subplot(2,1,2);
cla(sub2);
hk = animatedline(sub2);
ax2=gca;
ylabel( [{"рассогласования уровня"};{" в точке измерения"};{" относительно начального"};{" состояния, м"}]);
xlabel("контрольная скважина,   время, час");
grid on
drawnow
tau=0;
del1=[];
a = tic; % start timer

while 1 
for i=1:tau1
    for j=1:tau2
        %dH
        for x_i=2:Nx-1
            for y_i=2:Ny-1
             %проверка попадания на скважину
             del=0;
                 if abs(x_i-Xs1)<1 && abs(y_i-Ys1)<1 
             del=1;%
             end
                for z_i=2:Nz-1
                    dH1(x_i,y_i,z_i) = single( dtau*(1/eta1)*(kx1*(H1(x_i-1,y_i,z_i)-2*H1(x_i,y_i,z_i)+H1(x_i+1,y_i,z_i))/(dx^2)+...
                        ky1*(H1(x_i-1,y_i,z_i)-2*H1(x_i,y_i,z_i)+H1(x_i+1,y_i,z_i))/(dy^2)+...
                        kz1*(H1(x_i-1,y_i,z_i)-2*H1(x_i,y_i,z_i)+H1(x_i+1,y_i,z_i))/(dz1^2)) );
                    dH2(x_i,y_i,z_i) = single( dtau*(1/eta2)*(kx2*(H2(x_i-1,y_i,z_i)-2*H2(x_i,y_i,z_i)+H2(x_i+1,y_i,z_i))/(dx^2) + ...
                        ky2*(H2(x_i-1,y_i,z_i)-2*H2(x_i,y_i,z_i)+H2(x_i+1,y_i,z_i))/(dy^2) + ...
                        kz2*(H2(x_i-1,y_i,z_i)-2*H2(x_i,y_i,z_i)+H2(x_i+1,y_i,z_i))/(dz2^2)) + V*dtau*del );
                end
            end        
        end
        %расчёт H
        for x_i=2:Nx-1
            for y_i=2:Ny-1
                for z_i=2:Nz-1
                    H1(x_i,y_i,z_i)=single( H1(x_i,y_i,z_i)+dH1(x_i,y_i,z_i) );
                    H2(x_i,y_i,z_i)=single( H2(x_i,y_i,z_i)+dH2(x_i,y_i,z_i) );
                end
            end        
        end
        %для того, чтобы производняа на границах была равна нулю
        for x_i=1:Nx
            for y_i=1:Ny
            H1(x_i,y_i,1)=H1(x_i,y_i,2);
            H2(x_i,y_i,Nz)=H2(x_i,y_i,Nz-1);
            end
        end
        %прописать граничные условия
        for x_i=2:Nx-1
            for y_i=2:Ny-1
                H1(x_i,y_i,Nz)=single(H1(x_i,y_i,Nz)+b1*dtau*(H2(x_i,y_i,1)-H1(x_i,y_i,Nz)) ); %давлеие в верхнем уменьшиться и произойдё перетекание
                H2(x_i,y_i,1)=single( H2(x_i,y_i,1)-b1*dtau*(H2(x_i,y_i,1)-H1(x_i,y_i,Nz)) );                
            end
        end

        addpoints(h,tau,H2(Xs1,Ys1,Zs2)-TK(Xs1,Ys1,Zs2)); %рассогласование относительно начального состояния в добывающей скважине
        addpoints(hk,tau,H2(Xs2,Ys2,Zs2k)-TK(Xs1,Ys1,Zs2)); %рассогласование относительно начального состояния в контрольеой скважине       
        b = toc(a); % check timer
        if b > (1/30)
            drawnow % update screen every 1/30 seconds
            a = tic; % reset timer after updating
        end
        tau=double(tau+dtau);
    end
end
new_time_in_hour=[0:3600:tau];
new_time_in_hour2=[0:1:ceil(tau/3600)];
set(ax1, 'XTick', new_time_in_hour,'XTickLabel',new_time_in_hour2);
set(ax2, 'XTick', new_time_in_hour,'XTickLabel',new_time_in_hour2);
drawnow % draw final frame

break
    
    
end

%h2=plot(fig,new_time_in_hour,H2d); 

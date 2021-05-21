clear
clc
close all

N=12;%номер варианта 

Lx=5080; %meter
Ly=4000; %meter
Lz1=20+2*N; %meter, толщина I пласта
Lz2=45+5*N;%meter, толщина II пласта

%Коэффициенты передачи, которые должны получиться в скорректированной модели:
K1=0.02+0.01*N;
K2=0.002+0.001*N;

%object description
%ыерхний пласт
kx1=148.4/86400; % выбором kx можно обеспечить скорость 
ky1=24.64/86400; %  
kz1=14.64/86400;
eta1=0.000128;

%нижний пласт %размерность
kx2=148.4/86400;%0.1054  
ky2=24.64/86400;
kz2=14.64/86400;
eta2=0.00012;
b1=0.00128/86400;%коэффициент перетекания (перевести в секунды)
%для  дискретной модели (можно путём изменения добиться нужного k2 но это
%доллго

Q=100; %m^3/сек дебет, входное воздействие, в секундах
Qs=Q/86400;
Kf=0.0001113987771965157; %изменением этого коэффициента мы  получаем требу
%емый коэффициент передачи
V=Kf*Qs;

%кол-во точкек дискретизации %Nx,Ny,Nz по 15
Nx=50;
Ny=10;
Nz=9; %кол-во точек на каждый пласт
%добывающая скважина в нижнем пласте 

%выбор этапа расчёта (ОБЯЗАТЕЛЬНО выбрать что-то)
KOEF_EN = 0; %вывести значение коэффициентов передачи 
CONTROL_EN = 2;%изобразить контрольную скважину на графике под графиком добывающей
PID_EN=0;%добавить ПИД регулятор
    %усилительные звенья для ПИД регулятора:
    %dphi=pi/3 %запас по фазе, использованный при расчёте
    E1=495;
    E2=0.8795;
    E4=0.01137;
    Zad=-0.2*N;%задающее изменение понижение уровня
%включить расчёт времени запаздывания
PHASE_EN=0;
AXES_IN_HOUR=1;%оси в часах , если нет, то в секундах
%частота колебаний добавляемого сигнала
w=0.001;

%коэффицента аппроксимирующей динамической системы
%необходима для расчёта параметров ПИД-регуятора
Dz=1;
a1=0.00765;
K=0.1422;

if PID_EN
tau1=1; % период вывода графика в сек * tau1
tau2=1500; %период вывода графика в сек (=час)
else
%моделирование процесса в скважине
tau1=48; % период вывода графика в сек * tau2 (=часов)
tau2=3600; %период вывода графика в сек (=1час)
end
dtau=1; %шаг, равен  1 секунде
%два цикла сделны для вывода графиков

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

%расчёт шага в пространственных координатах
dx=Lx/(Nx-1); %шаг
dy=Ly/(Ny-1);
dz1=Lz1/(Nz-1);%равны
dz2=Lz2/(Nz-1);

%% инициализация массивов
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

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');%шрифт
fig=figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
if CONTROL_EN
sub1=subplot(2,1,1);
cla(sub1);
h = animatedline(sub1);
else
h = animatedline;
end
ax1=gca;
title("процессы в скважинах");
if CONTROL_EN
ylabel( [{"рассогласования уровня"};{" в точке измерения"};{" относительно начального"};{" состояния, м"}]);
else
ylabel('error=-\Delta H_{x_0,y_0,z_0,t} +\Delta H_{zad}, м.');
end
if AXES_IN_HOUR
xlabel("возмущающая скважина,   время, час");
else
xlabel("возмущающая скважина,   время, секунды");
end
grid on
drawnow
if CONTROL_EN
    %контрольная 
    sub2=subplot(2,1,2);
    cla(sub2);
    hk = animatedline(sub2);
    ax2=gca;
    ylabel( [{"рассогласования уровня"};{" в точке измерения"};{" относительно начального"};{" состояния, м"}]);
    if AXES_IN_HOUR
    xlabel("контрольная скважина,   время, час");
    else
    xlabel("контрольная скважина,   время, секунды");
    end
    grid on
    drawnow
end

%инициализация статовых значений
tau=0; %в сек
a = tic; % start timer

%инициализация значений  для регулятора
error = [0 0];
S=0;

for i=1:tau1
    for j=1:tau2
        %расчёт рипращения уровня понижения в скважине dH
        for x_i=2:Nx-1
            for y_i=2:Ny-1
             
             if PID_EN
                 %реализация ПИД регулятора
                 %рассогласование 
                 error(2)=(TK(Xs1,Ys1,Zs2)+Zad)-H2(Xs1,Ys1,Zs2);
                 
                 %интегратор 
                 S=double(S+error(2)*dtau);
                 
                 %формирование задающего воздействия Qs
                 Qs=(E1*error(2)+E2*(error(2)-error(1))/dtau+E4*S);

                 %сохранение текущего значения ошибки 
                 error(1)=error(2);
             end
             
             %проверка попадания на скважину
             del=0;
             if abs(x_i-Xs1)<1 && abs(y_i-Ys1)<1 
             del=1;%
             end
             
                for z_i=2:Nz-1
                    dH1(x_i,y_i,z_i) = double( dtau*(1/eta1)*(kx1*(H1(x_i-1,y_i,z_i)-2*H1(x_i,y_i,z_i)+H1(x_i+1,y_i,z_i))/(dx^2)+...
                        ky1*(H1(x_i-1,y_i,z_i)-2*H1(x_i,y_i,z_i)+H1(x_i+1,y_i,z_i))/(dy^2)+...
                        kz1*(H1(x_i-1,y_i,z_i)-2*H1(x_i,y_i,z_i)+H1(x_i+1,y_i,z_i))/(dz1^2)) );
                    if ~PHASE_EN
                    dH2(x_i,y_i,z_i) = double( dtau*(1/eta2)*(kx2*(H2(x_i-1,y_i,z_i)-2*H2(x_i,y_i,z_i)+H2(x_i+1,y_i,z_i))/(dx^2) + ...
                        ky2*(H2(x_i-1,y_i,z_i)-2*H2(x_i,y_i,z_i)+H2(x_i+1,y_i,z_i))/(dy^2) + ...
                        kz2*(H2(x_i-1,y_i,z_i)-2*H2(x_i,y_i,z_i)+H2(x_i+1,y_i,z_i))/(dz2^2)) + Kf*Qs*dtau*del);
                    else
                    dH2(x_i,y_i,z_i) = double( dtau*(1/eta2)*(kx2*(H2(x_i-1,y_i,z_i)-2*H2(x_i,y_i,z_i)+H2(x_i+1,y_i,z_i))/(dx^2) + ...
                        ky2*(H2(x_i-1,y_i,z_i)-2*H2(x_i,y_i,z_i)+H2(x_i+1,y_i,z_i))/(dy^2) + ...
                        kz2*(H2(x_i-1,y_i,z_i)-2*H2(x_i,y_i,z_i)+H2(x_i+1,y_i,z_i))/(dz2^2)) + Kf*Qs*sin(w*tau)*dtau*del);    
                    end
                end
            end        
        end
        %расчёт H
        for x_i=2:Nx-1
            for y_i=2:Ny-1
                for z_i=2:Nz-1
                    H1(x_i,y_i,z_i)=double( H1(x_i,y_i,z_i)+dH1(x_i,y_i,z_i) );
                    H2(x_i,y_i,z_i)=double( H2(x_i,y_i,z_i)+dH2(x_i,y_i,z_i) );
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
        %прописать граничные условия %значения граничных условий из условий Дарси
        for x_i=2:Nx-1
            for y_i=2:Ny-1
                H1(x_i,y_i,Nz)=double(H1(x_i,y_i,Nz)+b1*dtau*(H2(x_i,y_i,1)-H1(x_i,y_i,Nz)) ); %давлеие в верхнем уменьшиться и произойдё перетекание
                H2(x_i,y_i,1)=double( H2(x_i,y_i,1)-b1*dtau*(H2(x_i,y_i,1)-H1(x_i,y_i,Nz)) );                
            end
        end
    
    if PID_EN
    addpoints(h,tau,error(2)); %рассогласование относительно начального состояния в добывающей скважине
    else
    addpoints(h,tau,H2(Xs1,Ys1,Zs2)-TK(Xs1,Ys1,Zs2));    
    end
    if CONTROL_EN
    addpoints(hk,tau,H2(Xs2,Ys2,Zs2k)-TK(Xs2,Ys2,Zs2k)); %рассогласование относительно начального состояния в контрольеой скважине       
    end
    
    b = toc(a); % check timer
    if b > (1/30)
        drawnow % update screen every 1/30 seconds
        a = tic; % reset timer after updating
    end
    tau=double(tau+dtau); 
    end

end
if KOEF_EN
H1ust=H2(Xs1,Ys1,Zs2)-TK(Xs1,Ys1,Zs2)
H2ust=H2(Xs2,Ys2,Zs2k)-TK(Xs2,Ys2,Zs2k)
Qs
K1=abs((H2(Xs1,Ys1,Zs2)-TK(Xs1,Ys1,Zs2)))/Qs%коэффициент передачи по возм. скважине
K2=abs((H2(Xs2,Ys2,Zs2k)-TK(Xs2,Ys2,Zs2k)))/Qs%коэффициент передачи по контр. скважине
end
if AXES_IN_HOUR
new_time_in_hour=[0:3600:tau];
new_time_in_hour2=[0:1:ceil(tau/3600)];
set(ax1, 'XTick', new_time_in_hour,'XTickLabel',new_time_in_hour2);
    if CONTROL_EN
    set(ax2, 'XTick', new_time_in_hour,'XTickLabel',new_time_in_hour2);
    end
end
drawnow % draw final frame
if PHASE_EN
    dt=169329-168873
    dphi=dt*w
end
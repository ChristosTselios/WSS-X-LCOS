clear all; clc
Eout = 0; pixel = 3.4e-6;
omega=10e-6;lamda=1550e-9; 
[phase_set] = phase_set_creation(100);
iter = 0; theta=0.00700701; fiber_port=1.25e-3; fl = 0.025;
theta = linspace(-10,10,1000); iter = 0;
for cou = 1:length(theta)
for pix_count = -50:49
    iter = iter + 1;
    Eout = Eout + exp(-1*j*(phase_set(iter))*(erfz((lamda*pixel*(2*pix_count+1)-1*j*pi*omega^2*sin(theta(cou)))/(sqrt(2)*lamda*omega))-...
    erfz((lamda*pixel*(2*pix_count-1)-1*j*pi*omega^2*sin(theta(cou)))/(sqrt(2)*lamda*omega))));
end

output_efficiency = 0.5*exp(-0.5*(pi*omega*sin(theta(cou))/lamda)^2)*Eout;

cp(cou) = -10*log(abs(output_efficiency).^2);

Eout = 0; iter = 0; output_efficiency = 0;
end
figure; plot(theta,cp)
%%
clear all; clc

Eout = 0; Eout2 = 0; pixel = 4.4e-6;
omega=10e-6;lamda=1550e-9; 
%[phase_st] = phase_set_creation(100);
iter = 0; theta=3; fiber_port=1.25e-3; fl = 0.025;

ports=-20:20;
jj=1;
for jj = 1:length(ports)

    %lamda_factor(jj) = 6;
    lamda_factor(jj) = (lamda*fl)/((fiber_port*sin(3))*pixel);

for pix_count = -500:500
    
    iter = iter + 1;
    
    A_term = sqrt(2/pi)*(pixel/omega)*sinc((pixel/lamda)*tan(3));
    
    phase_d2(iter) = exp(-1*j*(2*pi/lamda_factor(jj))*pix_count);
    phase_d = exp(-1*j*(2*pi/lamda_factor(jj))*pix_count);
    C_factor2(iter)=erfi((lamda*pixel*(2*pix_count+1)-1*j*pi*omega^2*sin(3))/(sqrt(2)*lamda*omega));
    C_factor3(iter)=erfi((lamda*pixel*(2*pix_count-1)-1*j*pi*omega^2*sin(3))/(sqrt(2)*lamda*omega));
    C_factor=erfi((lamda*pixel*(2*pix_count+1)-1*j*pi*omega^2*sin(3))/(sqrt(2)*lamda*omega))-...
    erfz((lamda*pixel*(2*pix_count-1)-1*j*pi*omega^2*sin(3))/(sqrt(2)*lamda*omega));
    Eout = Eout + A_term*exp(-1*j*(2*pi/lamda_factor(jj))*pix_count)*C_factor;

    Eout2 = Eout2 + A_term*phase_d*exp(-(2*(pix_count*pixel)^2/omega^2)+1*j*(2*pi/lamda)*tan(3)*pix_count*pixel);

end
%output_efficiency = 0.5*exp(-0.5*(pi*omega*sin(theta)/lamda)^2)*Eout;
cp(jj) = 10*log(abs(Eout).^2);
cp2(jj) = 10*log(abs(Eout2).^2);
Eout = 0; iter = 0; output_efficiency = 0; Eout2 = 0;C_factor=0;
end

figure; plot(ports,cp,'linewidth',2)

xlabel ('Diffraction angle $(deg.)$','Interpreter','latex');
ylabel ('Tranmittance $(dB)$','Interpreter','latex');
set(gca,'Linewidth',2)
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))

figure; plot(ports,cp2,'linewidth',2)
ax=gca;
                    ax.FontSize=14;

                   pos=get(gca,'pos');
                   set(gca,'pos',[pos(1) pos(2) pos(3) pos(4)*.95]);
xlabel ('Diffraction angle $(deg.)$','Interpreter','latex');
ylabel ('Tranmittance $(dB)$','Interpreter','latex');
set(gca,'Linewidth',2)
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
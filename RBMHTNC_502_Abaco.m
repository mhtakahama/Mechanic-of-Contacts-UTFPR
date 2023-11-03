%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modelo de previs�o de temperatura para a bancada RBPLR
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Doutorado em Engenharia Mecanica PPGEM
%Universidade Tecnologica Federal do Parana (UTFPR) - Campus Curitiba
%Doutorando: Marcos Hiroshi Takahama
%Professor: Tiago Cousseau , Cezar Otaviano Ribeiro Negrao
%abr 2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4 Calculo de Ca e k

function [Ca_determinado,k_determinado,plt_numb]=RBMHTNC_502_Abaco(AB,plt_numb);
% Hertz_Factores
%Teoria de Hertz
%Determinação dos parametros K, Ca, Cdelta e Csigma
vx=(1.0:0.01:100.0); %(min:discretização:máximo)
M=zeros(size(vx,2),8);

for jx=1:size(vx,2);
    x=vx(jx);
    k0=x^(-2/3);
    k=k0;
    for i=1:30
        e2=1-k.*k;
        e=sqrt(e2);
        [Ke,Ee]=ellipke(e2,0.00001);
        kn=sqrt(Ee/((Ke-Ee)*x+Ke));
        erro=100*abs(kn-k);
        if erro<1e-5;
            break
        else
            k=kn;
        end
    end
    k=k;
    e2=1-k*k;
    e=sqrt(1-k*k);
    [Ke,Ee]=ellipke(e2,0.00001);
    Ca=(3*k*Ee/2/pi)^0.333333;
    Cd=3*k*Ke/2;
    Cs=(3*k/2/pi)/(Ca^3);
    M(jx,1)=x;
    M(jx,2)=k;
    M(jx,3)=e;
    M(jx,4)=Ke;
    M(jx,5)=Ee;
    M(jx,6)=Ca;
    M(jx,7)=Cd;
    M(jx,8)=Cs;
end


% Calculo do ponto individual
x=AB;
k0=x^(-2/3);
k=k0;
for i=1:30
    e2=1-k.*k;
    e=sqrt(e2);
    [Ke,Ee]=ellipke(e2,0.00001);
    kn=sqrt(Ee/((Ke-Ee)*x+Ke));
    erro=100*abs(kn-k);
    if erro<1e-5;
        break
    else
        k=kn;
    end
end
k_determinado=k;
e2=1-k_determinado*k_determinado;
e_determinado=sqrt(1-k_determinado*k_determinado);
[Ke_determinado,Ee_determinado]=ellipke(e2,0.00001);
Ca_determinado=(3*k_determinado*Ee_determinado/2/pi)^0.333333;

% figure
% set(gcf,'Name',['1 - Abaco Teoria de Hertz'],'NumberTitle','off');
% set(gcf, 'units','normalized','outerposition',[0 0 1 1]);%Maximiza a janela
% plot(M(:,1),M(:,6),'Color','#D95319','LineWidth',4); % Curva Ca do ábaco
% hold on
% plot(AB,Ca_determinado,'r*','LineWidth',5);
% 
% hold on
% plot(M(:,1),M(:,2),'Color','#0072BD','LineWidth',4); % Curva k do ábaco
% hold on
% plot(AB,k_determinado,'b*','LineWidth',5);
% 
% line([AB AB],[0 Ca_determinado],'Color','red')
% line([0 AB],[Ca_determinado Ca_determinado],'Color','red')
% line([AB AB],[0 k_determinado],'Color','blue')
% line([0 AB],[k_determinado k_determinado],'Color','blue')
% 
% grid on;
% xlim([1 20])
% ylim([0 1])
% legend('Curva Ca', 'Ponto calculado Ca','Curva k', 'Ponto calculado k')
% title(['Abaco de Ca e k para A/B = ' num2str(AB)],'FontSize',20) %Legend options
% xlabel('A/B')
% ylabel('Ca, k')

end
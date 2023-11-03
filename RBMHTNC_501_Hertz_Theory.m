%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modelo de previsão de temperatura para a bancada RBPLR
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Doutorado em Engenharia Mecanica PPGEM
%Universidade Tecnologica Federal do Parana (UTFPR) - Campus Curitiba
%Doutorando: Marcos Hiroshi Takahama
%Professor: Tiago Cousseau , Cezar Otaviano Ribeiro Negrao
%abr 2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5 Teoria de Hertz

function [Hertzproperties]=RBMHTNC_501_Hertz_Theory(F_carga,Rx1,Ry1,Rx2,Ry2,E1,E2)
% F_carga
% Rx1
% Ry1
% Rx2
% Ry2
% E1
% E2
%Calculo de a,b,p0 e pm

% Geometry constants [m]
Hertzproperties.Rx1=Rx1; % ball
Hertzproperties.Ry1=Ry1; % ball
Hertzproperties.Rx2=Rx2; % disc %buscar de livros
Hertzproperties.Ry2=Ry2; % disc %puxar do perfil geometrico do rolamento 6310
% Equivalent Radius
Hertzproperties.Rx=abs((0.5*(1/Hertzproperties.Rx1+1/Hertzproperties.Rx2))^-1);
Hertzproperties.Ry=abs((0.5*(1/Hertzproperties.Ry1+1/Hertzproperties.Ry2))^-1);
if Hertzproperties.Rx<=Hertzproperties.Ry
    Hertzproperties.A=1/Hertzproperties.Rx;
    Hertzproperties.B=1/Hertzproperties.Ry;
elseif Hertzproperties.Rx>Hertzproperties.Ry
    Hertzproperties.A=1/Hertzproperties.Ry;
    Hertzproperties.B=1/Hertzproperties.Rx;
end
% Material Properties
% Elastic Modulus [Pa]
Hertzproperties.E1=E1; % ball
Hertzproperties.E2=E2; % disc
% Poisson coefficient [-]
Hertzproperties.v1=3e-1; % ball
Hertzproperties.v2=3e-1; % disc
Hertzproperties.E=((1-Hertzproperties.v1^2)/Hertzproperties.E1+(1-Hertzproperties.v2^2)/Hertzproperties.E2)^-1;
% Hertzproperties.A
% Hertzproperties.B

% Gerador do Ã¡baco
[Hertzproperties.Ca,Hertzproperties.k]=RBMHTNC_502_Abaco(Hertzproperties.A/Hertzproperties.B);

% Calculo de a e b
Hertzproperties.a=Hertzproperties.Ca*(F_carga/((Hertzproperties.A+Hertzproperties.B)*Hertzproperties.E))^(1/3);
Hertzproperties.b=Hertzproperties.a/Hertzproperties.k;
Hertzproperties.p0=3/2*F_carga/(pi*Hertzproperties.a*Hertzproperties.b);
Hertzproperties.pm=F_carga/(pi*Hertzproperties.a*Hertzproperties.b);
end
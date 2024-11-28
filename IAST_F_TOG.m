function IAST_F_TOG

% Yu Ting 2024/11/28
% references for equations for IAST models: Adsorption technology in water
% Treatment: Fundamentals, Proccesses, and Modeling.  Eckhard Worch (2012) Section 4.5.3. 
% references for the IAST-F model establishment for binary solute: https://www.youtube.com/watch?v=BsCSurcpRyE&t=17s
%
% Comments:
% Using binary IAST_F model from Rahmat Sunarya to construct IAST-F model
% for our model. 

clc; 
clear; 

     %%%%%%%%%%%  Sorption isotherm/kinetic properties for TOG  %%%%%%%%%%%%%%%

     % partition coefficient for PCE [L/kg]
     KF_PCE_TOG= 10987.528;
     % n ratio for Freudlich isotherm for PCE
     N_PCE_TOG = 0.4884; 

     % partition coefficient for TCE [L/kg]
     KF_TCE_TOG= 6287.82;
     % n ratio for Freudlich isotherm for TCE
     N_TCE_TOG = 0.5227; 

%%%%%%%%%%%%%%%%%%%%%% Function area %%%%%%%%%%%%%%%%%%%%%%%%

    function Qe_IAST=IAST_F(Ce)

        Intfun_PCE=@(c) KF_PCE_TOG.*c.^(N_PCE_TOG)./c;
        Intfun_TCE=@(c) KF_TCE_TOG.*c.^(N_TCE_TOG)./c;

        F_Z=@(z) (integral(Intfun_PCE,0,Ce./z)-integral(Intfun_TCE,0,Ce./(1-z))).^2;
        
        z=fsolve(F_Z,0.5);
        Z_PCE=z; 
        Z_TCE=1-Z_PCE;

        C0_PCE=Ce./Z_PCE;
        C0_TCE=Ce./Z_TCE;
        
        Q0_PCE= KF_PCE_TOG.*C0_PCE.^(N_PCE_TOG); 
        Q0_TCE= KF_TCE_TOG.*C0_TCE.^(N_TCE_TOG);

        qT=1./(Z_PCE./Q0_PCE+Z_TCE./Q0_TCE);

        Qe_PCE_IAST=qT.*Z_PCE;
        Qe_TCE_IAST=qT.*Z_TCE;

        Qe_IAST= [Qe_PCE_IAST Qe_TCE_IAST];

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input varaible
Ce_IAST=linspace(1e-10,1.0,500);

n=length(Ce_IAST); 
for i=1:n
    Qe_IAST(i,:)=IAST_F(Ce_IAST(i));
end
Qe_IAST_PCE=Qe_IAST(:,1);
Qe_IAST_TCE=Qe_IAST(:,2);

figure(1)

plot(Ce_IAST,Qe_IAST_PCE,'-');
hold on 
grid on
plot(Ce_IAST, KF_PCE_TOG.*Ce_IAST.^N_PCE_TOG, '--')

text(0.03,11000,' (a) PCE isotherm')
xlabel('Ce (mmol/L)')
ylabel('Qe PCE (mmol/Kg)')
legend ('IAST-Freundlich model', 'Monosolute Freundlich model','Location','se')
fontsize(14,"points")

figure(2)

plot(Ce_IAST,Qe_IAST_TCE,'-');
hold on 
grid on
plot(Ce_IAST, KF_TCE_TOG.*Ce_IAST.^N_TCE_TOG, '--')
text(0.03,6500,' (b) TCE isotherm')
xlabel('Ce (mmol/L)')
ylabel('Qe TCE (mmol/Kg)')
legend ('IAST-Freundlich model', 'Monosolute Freundlich model','Location','se')
fontsize(14,"points")
 
writematrix([Ce_IAST' Qe_IAST_PCE],'IAST_F_PCE_isotherm_TOG.xlsx')
writematrix([Ce_IAST' Qe_IAST_TCE],'IAST_F_TCE_isotherm_TOG.xlsx')

end
function TOG_20X_vit_IAST
% 2024 11 28 Yu Ting
% 

% For bisection method, reference:  https://www.youtube.com/watch?v=JX47_5h2FdE
% For 4th rugga-kutta method:   https://www.youtube.com/watch?v=WCoKrOQdNg0
% For RMSE and r square calculation: https://www.mathworks.com/matlabcentral/answers/478999-how-to-show-r-square-correlation-and-rmse-on-a-scatterplot
% reference for the IAST-F model establishment for binary solute:
% https://www.youtube.com/watch?v=BsCSurcpRyE&t=17s
% Notes: 
% This model accounts for biodechlorination kinetic (Monod kinetic), sorption kinetic (1st order) and IAST model for PCE and TCE to account for multisolute impact for activated carbon sorption. 

clc
clear

fig = figure;
left_color = [0 0 0];
right_color = [0 0.6 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);


%%%%%%%%%%%% Set up Basic Parameters %%%%%%%%%%%%%%%

     % volume of solution  [L]
     Vw = 0.5;
     % volume of gas phase  [L]
     Vg = 0.22;
     % Henry constant of PCE  [atm*L/mol]; 
     H_PCE = 17.7;
     % Henry constant of TCE  [atm*L/mol]; 
     H_TCE = 9.02;
     % Henry constant of cis-DCE  [atm*L/mol]; 
     H_cDCE = 4.08;
     % Henry constant of trans-DCE  [atm*L/mol]; 
     H_tDCE = 9.28;
     % Ideal gas constant  [atm*L/(mol*K)]; 
     R = 0.082;
     % Temperature  [K]; 
     T = 301;    %T=28oC
     % maximum utilization rate for PCE [mmol/d/cell];
     qmPCE = 7.4E-12;
     % maximum utilization rate for TCE [mmol/d/cell];
     qmTCE = 0.43E-12;
     % half-saturation constant for PCE [mmol/L]
     KPCE = 4.0E-02;
     % half-saturation constant for TCE [mmol/L]
     KTCE = 8.5E-02;
     % cell Yield for PCE [cell/mmol]
     YPCE = 4.34E+10;
     % cell Yield for TCE [cell/mmol]
     YTCE = 4.34E+10;
     % transDCE/(trans+cis) ratio 
     tratio = 0.57;
     % cisDCE/(trans+cis) ratio 
     cratio = 0.42;
     % decay term   [1/d]
     b=0.06;
     % weight of the adsorbent [kg]
     ma= 0.00010; 
     
     %%%%%%%%%%%  Sorption isotherm (IAST-F)/kinetic properties for TOG  %%%%%%%%%%%%%%%

     % partition coefficient for PCE [L/kg]
     KF_PCE_TOG= 9745.407;
     % n ratio for Freudlich isotherm for PCE
     N_PCE_TOG = 0.4821; 

     % partition coefficient for TCE [L/kg]
     KF_TCE_TOG= 2991.575;
     % n ratio for Freudlich isotherm for TCE
     N_TCE_TOG = 0.5488; 

     % partition coefficient for DCE [L/kg]
     KF_DCE_TOG= 122.0956;
     % n ratio for Freudlich isotherm for DCE
     N_DCE_TOG = 0.9629; 

     % sorption first-order kinetic constant for PCE [1/d] 
     ksorp_PCE_TOG = 0.77; 

     % sorption first-order kinetic constant for TCE [1/d] 
     ksorp_TCE_TOG = 0.89; 
     
     % sorption first-order kinetic constant for cDCE [1/d] 
     ksorp_cDCE_TOG = 4.25; 
     
     % sorption first-order kinetic constant for tDCE [1/d] 
     ksorp_tDCE_TOG = 4.25; 

 %%%%%%%%%%%  Sorption isotherm (IAST-F)/kinetic properties for OLC  %%%%%%%%%%%%%%%
     
     % partition coefficient for PCE [L/kg]
     KF_PCE_OLC= 11439.314;
     % n ratio for Freudlich isotherm for PCE
     N_PCE_OLC = 0.2718; 

     % partition coefficient for TCE [L/kg]
     KF_TCE_OLC= 3807.150;
     % n ratio for Freudlich isotherm for TCE
     N_TCE_OLC = 0.7654; 

     % partition coefficient for DCE [L/kg]
     KF_DCE_OLC= 128.2331;
     % n ratio for Freudlich isotherm for DCE
     N_DCE_OLC = 1.0834; 
     
     % sorption first-order kinetic constant for PCE [1/d] 
     ksorp_PCE_OLC = 2.40; 

     % sorption first-order kinetic constant for TCE [1/d] 
     ksorp_TCE_OLC = 2.16; 
     
     % sorption first-order kinetic constant for cDCE [1/d] 
     ksorp_cDCE_OLC = 3.90; 
     
     % sorption first-order kinetic constant for tDCE [1/d] 
     ksorp_tDCE_OLC = 3.90; 

%%%%%%%%%%%%%%%%%%%%%% Cw data input %%%%%%%%%%%%%%%%%%%%

% t = incubation time (days)

t = [
2
11
15
21
25
37
64
];

% PCE, TCE, cDCE, and tDCEs in the aq mM for replicate 1

repli_1=[
0.102	0.005	0.000	0.000
0.093	0.019	0.000	0.000
0.077	0.078	0.000	0.000
0.009	0.187	0.017	0.017
0.006	0.181	0.031	0.031
0.003	0.169	0.038	0.038
0.004	0.228	0.058	0.058
];


% PCE, TCE, cDCE, and tDCEs in the aq mM for replicate 2

repli_2=[
0.096	0.002	0.000	0.000
0.112	0.019	0.000	0.000
0.078	0.065	0.000	0.000
0.011	0.209	0.010	0.010
0.004	0.219	0.024	0.023
0.003	0.288	0.044	0.045
0.002	0.333	0.063	0.064
];



%%%%%%%%%%%%%%%%%%%% average and std calculation %%%%%%%%%%%%%

PCE = mean([repli_1(:,1)';repli_2(:,1)']);
PCE_std = std([repli_1(:,1)';repli_2(:,1)']);

TCE = mean([repli_1(:,2)';repli_2(:,2)']);
TCE_std = std([repli_1(:,2)';repli_2(:,2)']);

cDCE = mean([repli_1(:,3)';repli_2(:,3)']);
cDCE_std = std([repli_1(:,3)';repli_2(:,3)']);

tDCE = mean([repli_1(:,4)';repli_2(:,4)']);
tDCE_std = std([repli_1(:,4)';repli_2(:,4)']);


%%%%%%%%%%%%%%%%%%%%%%%%%%% modeling algorithm by 4th order Runga Kutta method and bisection method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up inital Cw condition
t_sim=0; 
PCE_sim= 0.1322;  %%%%set up Initial PCE concentration from theoratical calculation
TCE_sim=0;
cDCE_sim=0;
tDCE_sim=0;
cell_sim=4.91E05;  % copies/ml
h=0.01;
t_end=80; 

% Assume sorption equlibrium at hte inital timepoint
% calculate initial CEs in the sorbed phase
QPCE = KF_PCE_TOG.*PCE_sim.^N_PCE_TOG; 
QTCE = 0; 
QcDCE = 0; 
QtDCE = 0;

% set up initial Cw,eq for CEs
PCE_eq = PCE_sim; 
TCE_eq = 0; 
cDCE_eq = 0; 
tDCE_eq = 0; 

% set up inital TOTCEs
TOTPCE0 = (Vw+(H_PCE*Vg/R/T)).*PCE_sim+ ma.*QPCE; 
TOTTCE0 = (Vw+(H_TCE*Vg/R/T)).*TCE_sim+ ma.*QTCE; 
TOTcDCE0 = (Vw+(H_cDCE*Vg/R/T)).*cDCE_sim+ ma.*QcDCE;
TOTtDCE0 = (Vw+(H_tDCE*Vg/R/T)).*tDCE_sim+ ma.*QtDCE;
TOTCEs0 = TOTPCE0 + TOTTCE0 + TOTcDCE0 +TOTtDCE0; 


% Define the ODE function for each variable (9 variables in total)
f_PCE = @(t_sim,PCE_sim,TCE_sim,cDCE_sim,tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE)   -ksorp_PCE_TOG.*(PCE_sim-PCE_eq)    -Vw./(Vw+(H_PCE*Vg/R/T)).*qmPCE.*(PCE_sim)./(PCE_sim+KPCE).*cell_sim.*1000;

f_TCE = @(t_sim,PCE_sim,TCE_sim,cDCE_sim,tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE)   -ksorp_TCE_TOG.*(TCE_sim-TCE_eq)    +Vw./(Vw+(H_TCE*Vg/R/T)).*qmPCE.*(PCE_sim)./(PCE_sim+KPCE).*cell_sim.*1000 ...
                                                                                                                                                          -Vw./(Vw+(H_TCE*Vg/R/T)).*qmTCE.*(TCE_sim)./(TCE_sim+KTCE).*cell_sim.*1000;
f_cDCE = @(t_sim,PCE_sim,TCE_sim,cDCE_sim,tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE)  -ksorp_cDCE_TOG.*(cDCE_sim-cDCE_eq) + cratio.*Vw./(Vw+(H_cDCE*Vg/R/T)).*qmTCE.*(TCE_sim)./(TCE_sim+KTCE).*cell_sim.*1000;

f_tDCE = @(t_sim,PCE_sim,TCE_sim,cDCE_sim,tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE)  -ksorp_tDCE_TOG.*(tDCE_sim-tDCE_eq) + tratio.*Vw./(Vw+(H_tDCE*Vg/R/T)).*qmTCE.*(TCE_sim)./(TCE_sim+KTCE).*cell_sim.*1000;

f_cell = @(t_sim,PCE_sim,TCE_sim,cDCE_sim,tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE)  YPCE.*qmPCE.*(PCE_sim)./(PCE_sim+KPCE).*cell_sim    +YTCE.*qmTCE.*(TCE_sim)./(TCE_sim+KTCE).*cell_sim   -b.*cell_sim ;

f_QPCE = @(t_sim,PCE_sim,TCE_sim,cDCE_sim,tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE)  (Vw+(H_PCE*Vg/R/T))./ma.*ksorp_PCE_TOG.*(PCE_sim - PCE_eq); 

f_QTCE = @(t_sim,PCE_sim,TCE_sim,cDCE_sim,tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE)  (Vw+(H_TCE*Vg/R/T))./ma.*ksorp_TCE_TOG.*(TCE_sim - TCE_eq); 

f_QcDCE = @(t_sim,PCE_sim,TCE_sim,cDCE_sim,tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE)  (Vw+(H_cDCE*Vg/R/T))./ma.*ksorp_cDCE_TOG.*(cDCE_sim - cDCE_eq); 

f_QtDCE = @(t_sim,PCE_sim,TCE_sim,cDCE_sim,tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE)  (Vw+(H_tDCE*Vg/R/T))./ma.*ksorp_tDCE_TOG.*(tDCE_sim - tDCE_eq); 



% 4th order Rugga-Kutta calculation for all variables
step=1;
while t_end-t_sim>=-10^(-10)

    % derivative of PCE calculation 
    k1=h.*f_PCE(t_sim,PCE_sim,TCE_sim,cDCE_sim, tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    k2=h.*f_PCE(t_sim+h/2,PCE_sim+k1/2,TCE_sim,cDCE_sim, tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    k3=h.*f_PCE(t_sim+h/2,PCE_sim+k2/2,TCE_sim,cDCE_sim, tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    k4=h.*f_PCE(t_sim+h,PCE_sim+k3,TCE_sim,cDCE_sim, tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    k=(k1+2.*k2+2.*k3+k4)./6;

    % derivative of TCE calculation 
    m1=h.*f_TCE(t_sim,PCE_sim,TCE_sim,cDCE_sim, tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    m2=h.*f_TCE(t_sim+h/2,PCE_sim,TCE_sim +m1/2,cDCE_sim, tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    m3=h.*f_TCE(t_sim+h/2,PCE_sim,TCE_sim +m2/2,cDCE_sim, tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    m4=h.*f_TCE(t_sim+h,PCE_sim,TCE_sim +m3,cDCE_sim, tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    m=(m1+2.*m2+2.*m3+m4)./6;

    % derivative of cDCE calculation 
    n1=h.*f_cDCE(t_sim,PCE_sim,TCE_sim,cDCE_sim, tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    n2=h.*f_cDCE(t_sim+h/2,PCE_sim,TCE_sim,cDCE_sim  +n1/2, tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    n3=h.*f_cDCE(t_sim+h/2,PCE_sim,TCE_sim,cDCE_sim  +n2/2, tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    n4=h.*f_cDCE(t_sim+h,PCE_sim,TCE_sim,cDCE_sim  +n3, tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    n=(n1+2.*n2+2.*n3+n4)./6;

    % derivative of tDCE calculation 
    o1=h.*f_tDCE(t_sim,PCE_sim,TCE_sim,cDCE_sim, tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    o2=h.*f_tDCE(t_sim+h/2,PCE_sim,TCE_sim,cDCE_sim,tDCE_sim  +o1/2,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    o3=h.*f_tDCE(t_sim+h/2,PCE_sim,TCE_sim,cDCE_sim,tDCE_sim   +o2/2,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    o4=h.*f_tDCE(t_sim+h,PCE_sim,TCE_sim,cDCE_sim, tDCE_sim   +o3,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    o=(o1+2.*o2+2.*o3+o4)./6;    

    % derivative of cell calculation 
    p1=h.*f_cell(t_sim,PCE_sim,TCE_sim,cDCE_sim, tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    p2=h.*f_cell(t_sim+h/2,PCE_sim,TCE_sim,cDCE_sim,tDCE_sim,cell_sim   +p1/2,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    p3=h.*f_cell(t_sim+h/2,PCE_sim,TCE_sim,cDCE_sim,tDCE_sim,cell_sim    +p2/2,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    p4=h.*f_cell(t_sim+h,PCE_sim,TCE_sim,cDCE_sim, tDCE_sim,cell_sim    +p3,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    p=(p1+2.*p2+2.*p3+p4)./6;    

    % derivative of QPCE calculation 
    q1=h.*f_QPCE(t_sim,PCE_sim,TCE_sim,cDCE_sim, tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    q2=h.*f_QPCE(t_sim+h/2,PCE_sim,TCE_sim,cDCE_sim,tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE +q1/2,QTCE,QcDCE,QtDCE);
    q3=h.*f_QPCE(t_sim+h/2,PCE_sim,TCE_sim,cDCE_sim,tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE +q2/2,QTCE,QcDCE,QtDCE);
    q4=h.*f_QPCE(t_sim+h,PCE_sim,TCE_sim,cDCE_sim, tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE +q3,QTCE,QcDCE,QtDCE);
    q=(q1+2.*q2+2.*q3+q4)./6;    


    % derivative of QTCE calculation 
    r1=h.*f_QTCE(t_sim,PCE_sim,TCE_sim,cDCE_sim, tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    r2=h.*f_QTCE(t_sim+h/2,PCE_sim,TCE_sim,cDCE_sim,tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE  +r1/2,QcDCE,QtDCE);
    r3=h.*f_QTCE(t_sim+h/2,PCE_sim,TCE_sim,cDCE_sim,tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE  +r2/2,QcDCE,QtDCE);
    r4=h.*f_QTCE(t_sim+h,PCE_sim,TCE_sim,cDCE_sim, tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE  +r3,QcDCE,QtDCE);
    r=(r1+2.*r2+2.*r3+r4)./6;      

    % derivative of QcDCE calculation 
    s1=h.*f_QcDCE(t_sim,PCE_sim,TCE_sim,cDCE_sim, tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    s2=h.*f_QcDCE(t_sim+h/2,PCE_sim,TCE_sim,cDCE_sim,tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE   +s1/2,QtDCE);
    s3=h.*f_QcDCE(t_sim+h/2,PCE_sim,TCE_sim,cDCE_sim,tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE   +s2/2,QtDCE);
    s4=h.*f_QcDCE(t_sim+h,PCE_sim,TCE_sim,cDCE_sim, tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE  +s3,QtDCE);
    s=(s1+2.*s2+2.*s3+s4)./6;      
    
    % derivative of QtDCE calculation 
    u1=h.*f_QtDCE(t_sim,PCE_sim,TCE_sim,cDCE_sim, tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE);
    u2=h.*f_QtDCE(t_sim+h/2,PCE_sim,TCE_sim,cDCE_sim,tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE    +u1/2);
    u3=h.*f_QtDCE(t_sim+h/2,PCE_sim,TCE_sim,cDCE_sim,tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE    +u2/2);
    u4=h.*f_QtDCE(t_sim+h,PCE_sim,TCE_sim,cDCE_sim, tDCE_sim,cell_sim,PCE_eq,TCE_eq,cDCE_eq,tDCE_eq,QPCE,QTCE,QcDCE,QtDCE +u3);
    u=(u1+2.*u2+2.*u3+u4)./6;      

    % calculate mass balance for Total PCE, TCE, DCE, and CEs
    TOTPCE = (Vw+(H_PCE*Vg/R/T)).*PCE_sim+ ma.*QPCE; 
    TOTTCE = (Vw+(H_TCE*Vg/R/T)).*TCE_sim+ ma.*QTCE; 
    TOTcDCE = (Vw+(H_cDCE*Vg/R/T)).*cDCE_sim+ ma.*QcDCE;
    TOTtDCE = (Vw+(H_tDCE*Vg/R/T)).*tDCE_sim+ ma.*QtDCE;
    TOTCEs = TOTPCE + TOTTCE + TOTcDCE +TOTtDCE; 

    % record for printing porpose
     fork(step,:)=[t_sim PCE_sim TCE_sim cDCE_sim tDCE_sim cell_sim PCE_eq TCE_eq cDCE_eq tDCE_eq QPCE QTCE QcDCE QtDCE h k m TOTPCE TOTTCE TOTcDCE TOTtDCE TOTCEs ]; 
    Fork_result = array2table(fork,'VariableNames', {'t','PCE', 'TCE','cDCE','tDCE','cell','PCE_eq', 'TCE_eq','cDCE_eq','tDCE_eq','QPCE', 'QTCE','QcDCE','QtDCE','h', 'k','m','TOTPCE', 'TOTTCE','TOTcDCE', 'TOTtDCE', 'TOTCEs ' }); 

    % update each of nine variables for the next step
    step= step+1; 
    t_sim = t_sim+h; 
    PCE_sim = PCE_sim+k; 
    TCE_sim = TCE_sim+m; 
    cDCE_sim = cDCE_sim+n; 
    tDCE_sim = tDCE_sim+o; 
    cell_sim = cell_sim+p; 
    QPCE = QPCE+q; 
    QTCE = QTCE+r;
    QcDCE = QcDCE+s; 
    QtDCE = QtDCE+u; 

    % mass balance correction for maintaning the consistent of the total CEs
    PCE_sim = PCE_sim.*(TOTCEs0./TOTCEs);
    TCE_sim = TCE_sim.*(TOTCEs0./TOTCEs);
    cDCE_sim = cDCE_sim.*(TOTCEs0./TOTCEs);
    tDCE_sim = tDCE_sim.*(TOTCEs0./TOTCEs);
    QPCE = QPCE.*(TOTCEs0./TOTCEs);
    QTCE = QTCE.*(TOTCEs0./TOTCEs);
    QcDCE = QcDCE.*(TOTCEs0./TOTCEs);
    QtDCE = QtDCE.*(TOTCEs0./TOTCEs);


    % Use bisection method to find PCE_eq for the next step
    f = @(x) (Vw+(H_PCE*Vg/R/T)).*(PCE_sim-x)+ ma.*(QPCE-KF_PCE_TOG.*x.^N_PCE_TOG);
    left  = 0; 
    right = 1; 
    n = 60; 
    e = 0.000001;

      if f(left)*f(right)<0
        for i = 1:n
              PCE_eq=(left+right)/2;

              if abs( PCE_eq-right)<e || abs( PCE_eq-left)<e
                   break
              end
              if f(left)*f( PCE_eq)<0
                  right= PCE_eq;
             elseif f(right)*f( PCE_eq)<0
                  left= PCE_eq;
     end 
        end
     end


    % Use bisection method to find TCE_eq for the next step
    f = @(x) (Vw+(H_TCE*Vg/R/T)).*(TCE_sim-x)+ ma.*(QTCE-KF_TCE_TOG.*x.^N_TCE_TOG);
    left  = 0; 
    right = 1; 
    n = 60; 
    e = 0.000001;

    if f(left)*f(right)<0
        for i = 1:n
              TCE_eq=(left+right)/2;

              if abs(TCE_eq-right)<e || abs(TCE_eq-left)<e
                   break
              end
              if f(left)*f(TCE_eq)<0
                  right=TCE_eq;
             elseif f(right)*f( TCE_eq)<0
                  left= TCE_eq;
     end 
        end
    end


    % Use bisection method to find cDCE_eq for the next step
    f = @(x) (Vw+(H_cDCE*Vg/R/T)).*(cDCE_sim-x)+ ma.*(QcDCE-KF_DCE_TOG.*x.^N_DCE_TOG);
    left  = 0; 
    right = 1; 
    n = 60; 
    e = 0.000001;

    if f(left)*f(right)<0
        for i = 1:n
              cDCE_eq=(left+right)/2;

              if abs(cDCE_eq-right)<e || abs(cDCE_eq-left)<e
                   break
              end
              if f(left)*f(cDCE_eq)<0
                  right=cDCE_eq;
             elseif f(right)*f(cDCE_eq)<0
                  left=cDCE_eq;
     end 
        end
    end

    % Use bisection method to find tDCE_eq for the next step
    f = @(x) (Vw+(H_tDCE*Vg/R/T)).*(tDCE_sim-x)+ ma.*(QtDCE-KF_DCE_TOG.*x.^N_DCE_TOG);
    left  = 0; 
    right = 1; 
    n = 60; 
    e = 0.000001;

    if f(left)*f(right)<0
        for i = 1:n
              tDCE_eq=(left+right)/2;

              if abs(tDCE_eq-right)<e || abs(tDCE_eq-left)<e
                   break
              end
              if f(left)*f(tDCE_eq)<0
                  right=tDCE_eq;
             elseif f(right)*f(tDCE_eq)<0
                  left=tDCE_eq;
     end 
        end
    end

end

Fork_result;
writetable(Fork_result,'Forth order Runga Kutta method.xlsx')


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculating error %%%%%%%%%%%%%%%%%%%%%%%%%
C(:,1)=PCE;
C(:,2)=TCE;
C(:,3)=cDCE;
C(:,4)=tDCE;
Step_sim=t/h+1;
RegressionLine = Fork_result(Step_sim,2:5);

% RMSE between regression line and Cw
RMSE = sqrt(mean((C-RegressionLine).^2));

% R2 between regression line and y
SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
SS_Y = sum((C-mean(C)).^2);
SS_XY = sum((RegressionLine-mean(RegressionLine)).*(C-mean(C)));
R_squared = SS_XY./sqrt(SS_X.*SS_Y);

RMSE=table2array(RMSE);
R_squared=table2array(R_squared);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ploting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


errorbar(t,PCE,PCE_std,'-bs','MarkerFaceColor','b', 'MarkerSize', 12,'MarkerEdgeColor',[0 0 0]);
hold on
errorbar(t,TCE,TCE_std,'-md','MarkerFaceColor','m', 'MarkerSize', 12, 'MarkerEdgeColor',[0 0 0]);
errorbar(t,cDCE,cDCE_std,'-kd','MarkerFaceColor', [0.6 0.6 0.6], 'MarkerSize', 12, 'MarkerEdgeColor',[0 0 0]);
errorbar(t,tDCE,tDCE_std ,'-kd','MarkerFaceColor', 'k', 'MarkerSize',12, 'MarkerEdgeColor',[0 0 0]);
plot5=plot(fork(:,1), fork(:,2),'--b','LineWidth',2);
plot6=plot(fork(:,1), fork(:,3),'--m','LineWidth',2);
plot7=plot(fork(:,1), fork(:,4),':','Color',[0.6,0.6,0.6],'LineWidth',2);
plot8=plot(fork(:,1), fork(:,5),'-.','Color',[0,0,0],'LineWidth',2);
hold off

xlim([0 80])
ylim([0 0.4])
%text(1,0.37,'(b) Bio-TOG 20X vit')
text(1,0.37,'Coal AC + vit')
xlabel('Incubation time (days)')
%ylabel('Aqueous concentration (mM)')
fontsize(20,"points")
%legend('PCE','TCE','cDCE','tDCE')

fprintf(1,'\tFitting errors:\n')
fprintf(1, '\t\t PCE RMSE: %0.4f | R2: %0.4f\n',RMSE(1),R_squared(1))
fprintf(1, '\t\t TCE RMSE: %0.4f | R2: %0.4f\n',RMSE(2),R_squared(2))
fprintf(1, '\t\t cDCE RMSE: %0.4f | R2: %0.4f\n',RMSE(3),R_squared(3))
fprintf(1, '\t\t tDCE RMSE: %0.4f | R2: %0.4f\n',RMSE(4),R_squared(4))

end


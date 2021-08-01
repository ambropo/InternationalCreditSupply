%--------------------------------------------------------------------------%
%                         Replication codes for                            %
%                   "International Credit Supply Shocks"                   %
%            by A. Cesa-Bianchi, A. Ferrero, and A. Rebucci                %
%                             December, 2017                               %
%                                                                          %
%--------------------------------------------------------------------------%
% This code replicates the results from the theoretical model in A. 
% Cesa-Bianchi, A. Ferrero, and A. Rebucci (2018) "International Credit 
% Supply Shocks," published in the the Journal of International Economics. 
%--------------------------------------------------------------------------%
% BEFORE RUNNING THE CODE: 
% Please make sure to add the folder "Codes" and all its subfolders to your 
% Matlab path. The code has been tested with Matlab R2020a on a Mac.
%--------------------------------------------------------------------------%

clear; clc; 
close all; warning ('off','all');
global f yH1 yF1 lambda yH2 yF2 kappa theta beta betas zeta zeta_nb chi chi_new eta phi
lambda=0.79;
yH1=1;
yF1=1;
yH2=1;
yF2=1;
kappa = 0.85;
theta = .92;
beta = 0.90;
betas = 0.99;
zeta = 0.02; % cost of adjusting equity (note: zeta=psi/ebar)
zeta_nb = 0.03;
chi = 0.1;
chi_new = 0.02;
eta = 0.43;  % target: median of fc liab (~0.7) => eta = 1/(share fc liab)-1=0.43
phi = 0.0001; % cost of hedging

support = (0.566:0.000025:.572)';
q_nb(1:length(support),1) = kappa;



%% COMPUTE EQUILIBRIUM
%==========================================================================
% Exchange rate (s1) 
s1 = ((lambda*yH1)./(lambda*yF1+(1-lambda).*(1+eta).*support)).^(1-lambda);

% Collateral constraint
CC = s1.*(1+eta).*support-theta.*q_nb; % CC=0 when CC is binding with equality
NB = find(CC<0);
B  = find(CC>0);

% Demand: Defined by the implicit function Rb(b,s1(b),s2(b,Rb)))
for ii=1:length(support)
    f = support(ii);
    RbD(ii,1) = real(fsolve(@f_Rb,1));
    RbD_nb(ii,1) = real(fsolve(@f_Rb_nb,1));
end

% Supply (in binding region)
RbS(1:length(support),1) = (1/betas)*(1+chi*zeta*(1+eta).*support) + phi/(1+eta).*support;
eq = find(abs(RbD-RbS)==min(abs(RbD-RbS)));

% Supply Binding (shock in binding region)
RbS_new(1:length(support),1) = (1/betas)*(1+chi_new*zeta*(1+eta).*support) + phi/(1+eta).*support;
eq_new = find(abs(RbD-RbS_new)==min(abs(RbD-RbS_new)));

% Supply (in non non-binding region)
RbS_nb(1:length(support),1) = (1/betas)*(1+chi*zeta_nb*(1+eta).*support) + phi/(1+eta).*support;
eq_nb = find(abs(RbD_nb-RbS_nb)==min(abs(RbD_nb-RbS_nb)));

% Exchange rate (s2)
s2 =    ((lambda*yH2)./(lambda*yF2-(1-lambda).*RbD.*(1+eta).*support)).^(1-lambda);
s2_nb = ((lambda*yH2)./(lambda*yF2-(1-lambda).*RbD_nb.*(1+eta).*support)).^(1-lambda);

% House Price
q = kappa./(1-theta+(theta*beta.*RbD.*s2)./s1);

% Consumption
c1 = s1.*(1+eta).*support + s1.^(-lambda/(1-lambda)).*yH1;


%% PLOT EQUILIBRIUM
%==========================================================================
% Annualize interest rates
RbD     = 100*(RbD.^4-1);
RbD_nb  = 100*(RbD_nb.^4-1);
RbS     = 100*(RbS.^4-1);
RbS_nb  = 100*(RbS_nb.^4-1);
RbS_new = 100*(RbS_new.^4-1);

figure
FigSize(14,9)
% Real rate (demand)
h1 = plot(support(B),RbD(B),'LineWidth',2.5,'Color',cmap(1)); hold on; 
plot(support(NB),RbD_nb(NB),'LineWidth',2.5,'Color',cmap(1)); hold on;
grid on; xlim([support(1) support(end)]); %axis tight;
% Real Rate (supply non binding)
h2 = plot(support,RbS_nb,'LineStyle',':','LineWidth',2.5,'Color',cmap(3)); hold on; 
plot(support(eq_nb),RbS_nb(eq_nb),'kp','LineWidth',0.01,'MarkerSize',12,'MarkerFaceColor',rgb('red')); hold on;
t = text(support(eq_nb),RbS_nb(eq_nb),{'A'}); t.FontSize = 13; t.Position(2) = 1.02*(t.Position(2)); t.FontWeight = 'bold';
grid on; xlim([support(1) support(end)]); %axis tight;
% aux = ylim; ylim([0 aux(2)]);
xlabel('Credit ($f$)','Interpreter','latex'); ylabel('Lending Rate ($R$)','Interpreter','latex');
legend([h1 h2],{'Demand','Supply (A)'},'Location','SouthWest')
% opt = LegOption; opt.handle = [h1 h2 h3]; LegPlot({'Demand','Supply (A)','Supply (B)'},opt);
SaveFigure('equil_nb',1)
clf('reset')

FigSize(14,9)
% Real rate (demand)
h1 = plot(support(B),RbD(B),'LineWidth',2.5,'Color',cmap(1)); hold on; 
plot(support(NB),RbD_nb(NB),'LineWidth',2.5,'Color',cmap(1)); hold on;
grid on; xlim([support(1) support(end)]); %axis tight;
% Real Rate (supply binding)
h3 = plot(support,RbS,'LineWidth',2.5,'Color',cmap(2)); hold on; 
plot(support(eq),RbS(eq),'kp','LineWidth',0.01,'MarkerSize',12,'MarkerFaceColor',rgb('red')); hold on; 
t = text(support(eq),RbS(eq),{'B'}); t.FontSize = 13; t.Position(2) = 1.02*(t.Position(2)); t.FontWeight = 'bold';
grid on; xlim([support(1) support(end)]); %axis tight;
% aux = ylim; ylim([0 aux(2)]);
xlabel('Credit ($f$)','Interpreter','latex'); ylabel('Lending Rate ($R$)','Interpreter','latex');
legend([h1 h3],{'Demand','Supply (Binding)'},'Location','SouthWest')
% opt = LegOption; opt.handle = [h1 h2 h3]; LegPlot({'Demand','Supply (A)','Supply (B)'},opt);
SaveFigure('equil_b',1)
clf('reset')

FigSize(14,9)
% Real rate (demand)
h1 = plot(support(B),RbD(B),'LineWidth',2.5,'Color',cmap(1)); hold on; 
plot(support(NB),RbD_nb(NB),'LineWidth',2.5,'Color',cmap(1)); hold on;
grid on; xlim([support(1) support(end)]); %axis tight;
% Real Rate (supply non binding)
h2 = plot(support,RbS_nb,'LineStyle',':','LineWidth',2.5,'Color',cmap(3)); hold on; 
plot(support(eq_nb),RbS_nb(eq_nb),'kp','LineWidth',0.01,'MarkerSize',12,'MarkerFaceColor',rgb('red')); hold on;
t = text(support(eq_nb),RbS_nb(eq_nb),{'A'}); t.FontSize = 13; t.Position(2) = 1.02*(t.Position(2)); t.FontWeight = 'bold';
% Real Rate (supply binding)
h3 = plot(support,RbS,'LineWidth',2.5,'Color',cmap(2)); hold on; 
plot(support(eq),RbS(eq),'kp','LineWidth',0.01,'MarkerSize',12,'MarkerFaceColor',rgb('red')); hold on; 
t = text(support(eq),RbS(eq),{'B'}); t.FontSize = 13; t.Position(2) = 1.02*(t.Position(2)); t.FontWeight = 'bold';
grid on; xlim([support(1) support(end)]); %axis tight;
% aux = ylim; ylim([0 aux(2)]);
xlabel('Credit ($f$)','Interpreter','latex'); ylabel('Lending Rate ($R$)','Interpreter','latex');
legend([h1 h2 h3],{'Demand','Supply (A)','Supply (B)'},'Location','SouthWest')
% opt = LegOption; opt.handle = [h1 h2 h3]; LegPlot({'Demand','Supply (A)','Supply (B)'},opt);
SaveFigure('equil',1)
clf('reset')

Rd = 100*((1/betas)^4-1);
Re = 100*(((1+zeta*(1+eta)*support(eq))/betas)^4-1);
Re_nb = 100*(((1+zeta_nb*(1+eta)*support(eq))/betas)^4-1);
Rb = chi*Re + (1-chi)*Rd;
Rb_nb = chi*Re_nb + (1-chi)*Rd;
Rb_lc = s2(eq)*Rb/s1(eq);
Rb_lc_nb = s2(eq_nb)*Rb_nb/s1(eq_nb);
EP = Re-Rd;
EP_nb = Re_nb-Rd;
disp(['Rd is equal to : ' num2str(Rd)])
disp(['Rb and Rb_nb are equal to : ' num2str(Rb) ' and ' num2str(Rb_nb)])
disp(['Rb in LC and Rb_nb in LC are equal to : ' num2str(Rb_lc) ' and ' num2str(Rb_lc_nb)])
disp(['Re and Re_nb are equal to : ' num2str(Re) ' and ' num2str(Re_nb)])
disp(['EP and EP_nb are equal to : ' num2str(EP) ' and ' num2str(EP_nb)])


%% PLOT ALL
%==========================================================================
figure 
FigSize(22,12)
% Real rate (demand)
subplot(2,2,1); 
plot(support(B),RbD(B),'LineWidth',2.5,'Color',cmap(1)); hold on; 
plot(support(NB),RbD_nb(NB),'LineWidth',2.5,'Color',cmap(1)); hold on;
title('(a) Lending Rate','Fontweight','Normal'); grid on; xlim([support(1) support(end)]); %axis tight;
% Real Rate (supply)
subplot(2,2,1); 
plot(support,RbS,'LineWidth',2.5,'Color',cmap(2));
hold on; plot(support(eq),RbS(eq),'kp','LineWidth',0.01,'MarkerSize',12,'MarkerFaceColor',rgb('red'));
t = text(support(eq),RbS(eq),{'B'}); t.FontSize = 13; t.Position(2) = 1.05*(t.Position(2)); t.FontWeight = 'bold';
hold on; plot(support,RbS_new,'LineWidth',2.5,'LineStyle',':','Color',cmap(2));
hold on; plot(support(eq_new),RbS_new(eq_new),'ko','LineWidth',0.01,'MarkerSize',8,'MarkerFaceColor',rgb('red'));
t = text(support(eq_new),RbS_new(eq_new),{'B'''}); t.FontSize = 13; t.Position(2) = 1.05*(t.Position(2)); t.FontWeight = 'bold';
ylabel('$R$','Interpreter','latex')
% aux = ylim; ylim([0 aux(2)]);
% Exchange rate
subplot(2,2,2); 
plot(support,s1,'LineWidth',2.5,'Color',cmap(1)); 
hold on; plot(support(eq),s1(eq),'kp','LineWidth',0.01,'MarkerSize',12,'MarkerFaceColor',rgb('red'));
t = text(support(eq),s1(eq),{'B'}); t.FontSize = 13; t.Position(2) = 1.00006*(t.Position(2)); t.FontWeight = 'bold';
hold on; plot(support(eq_new),s1(eq_new),'ko','LineWidth',0.01,'MarkerSize',8,'MarkerFaceColor',rgb('red'));
t = text(support(eq_new),s1(eq_new),{'B'''}); t.FontSize = 13; t.Position(2) = 1.00006*(t.Position(2)); t.FontWeight = 'bold';
title('(b) Exch. Rate','Fontweight','Normal'); grid on; xlim([support(1) support(end)]); %axis tight;
ylabel('$s_1$','Interpreter','latex')
% House Price
subplot(2,2,3); 
plot(support(B),q(B),'LineWidth',2.5,'Color',cmap(1)); hold on; 
plot(support(NB),q_nb(NB),'LineWidth',2.5,'Color',cmap(1));
hold on; plot(support(eq),q(eq),'kp','LineWidth',0.01,'MarkerSize',12,'MarkerFaceColor',rgb('red'));
t = text(support(eq),q(eq),{'B'}); t.FontSize = 13; t.Position(2) = .9997*(t.Position(2)); t.FontWeight = 'bold';
hold on; plot(support(eq_new),q(eq_new),'ko','LineWidth',0.01,'MarkerSize',8,'MarkerFaceColor',rgb('red'));
t = text(support(eq_new),q(eq_new),{'B'''}); t.FontSize = 13; t.Position(2) = .9997*(t.Position(2)); t.FontWeight = 'bold';
title('(c) House Price','Fontweight','Normal'); grid on; xlim([support(1) support(end)]); %axis tight;
ylabel('$q$','Interpreter','latex')
% Consumption
subplot(2,2,4); 
plot(support,c1,'LineWidth',2.5,'Color',cmap(1));
hold on; plot(support(eq),c1(eq),'kp','LineWidth',0.01,'MarkerSize',12,'MarkerFaceColor',rgb('red'));
t = text(support(eq),c1(eq),{'B'}); t.FontSize = 13; t.Position(2) = 0.9995*(t.Position(2)); t.FontWeight = 'bold';
hold on; plot(support(eq_new),c1(eq_new),'ko','LineWidth',0.01,'MarkerSize',8,'MarkerFaceColor',rgb('red'));
t = text(support(eq_new),c1(eq_new),{'B'''}); t.FontSize = 13; t.Position(2) = 0.9995*(t.Position(2)); t.FontWeight = 'bold';
title('(d) Consumption','Fontweight','Normal'); grid on; xlim([support(1) support(end)]); %axis tight;
ylabel('$c_1$','Interpreter','latex')
SaveFigure('equil_shock',1)
close all

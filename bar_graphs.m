%%
%mutual repression with background production
clear all;
clc;


maxrange=3; %parameters vary logarithmically from 10^-maxrange to 10^maxrange
npts=1e7; %total sets of randomised parameters checked for antagonism
k=1;

parammax=1e3;
jp=0;
for i=1:npts


    beta1=10^(3*(-1+2*rand)); %same as beta_T
    beta2=10^(3*(-1+2*rand)); %same as beta_E
    Theta1=10^(3*(-1+2*rand)); %This is the same as mu_T in calculations
    Theta2=10^(3*(-1+2*rand)); %This is the same as mu_E in calculations
    eta1=10^(3*(-1+2*rand)); %Same as eta_T
    eta2=10^(3*(-1+2*rand)); %Same as eta_E
    

%\Delta M calculations
SIGnabT=eta1*Theta1-eta2*beta1/(1+beta1)^2; %only TGF-\beta gradient present 
SIGnabE=eta2*Theta2-eta1*beta2/(1+beta2)^2; %only EGF gradient present
SIGnabTnabE=eta1*(Theta1-beta2)/(1+beta2)^2+eta2*(Theta2-beta1)/(1+beta1)^2; %both TGF-\beta
%and EGF gradient present
SIGnabTE0=eta1*Theta1/(1+beta2)-eta2*(1+Theta2)*beta1/(1+beta1)^2; %TGF-\beta gradient
%on EGF background
SIGnabTT0=(eta1*Theta1-eta2*beta1/(1+10*beta1)^2); %TGF-\beta gradient on large TGF-\beta 
%background


%M calculation
vCtrl=eta1+eta2; %No gradient present
vnabT=eta1*(1+Theta1)+eta2/(1+beta1);%only TGF-\beta gradient present 
vnabE=eta2*(1+Theta2)+eta1/(1+beta2);%only EGF gradient present
vnabTnabE=eta1*(1+Theta1)/(1+beta2)+eta2*(1+Theta2)/(1+beta1); %both TGF-\beta
%and EGF gradient present
vnabTE0=vnabTnabE; %TGF-\beta gradient
%on EGF background
vnabTT0=eta1*(1+10*Theta1)+eta2/(1+10*beta1); %TGF-\beta gradient on large TGF-\beta 
%background



if SIGnabT>0 && SIGnabE>0 && SIGnabTnabE<min(SIGnabT,SIGnabE) %select parameters if
    %antagonism in \Delta M is shown
    betx(k)=beta1;
    bety(k)=beta2;
    phx(k)=Theta1;
    phy(k)=Theta2;
    et1(k)=eta1;
    et2(k)=eta2;
    
    %\Delta M (\equiv CI) values
    CIT(k)=SIGnabT;
    CIE(k)=SIGnabE;
    CITE(k)=SIGnabTnabE;
    CITE0(k)=SIGnabTE0;
    CITT0(k)=SIGnabTT0;
    
    %M (\equiv Speed) values
    vC(k)=vCtrl;
    vT(k)=vnabT;
    vE(k)=vnabE;
    vTE(k)=vnabTnabE;
    vTE0(k)=vnabTE0;
    vTT0(k)=vnabTT0;
    k=k+1;
end

end


%Normalising the \Delta M to the TGF-\beta gradient only case
CITbar=mean(CIT)/mean(CIT);
CIEbar=mean(CIE)/mean(CIT);
CITEbar=mean(CITE)/mean(CIT);
CITE0bar=mean(CITE0)/mean(CIT);
CITT0bar=mean(CITT0)/mean(CIT);

%Normalising the M to the TGF-\beta gradient only case
vTbar=mean(vT)/mean(vT);
vCbar=mean(vC)/mean(vT);
vEbar=mean(vE)/mean(vT);
vTEbar=mean(vTE)/mean(vT);
vTE0bar=mean(vTE0)/mean(vT);
vTT0bar=mean(vTT0)/mean(vT);

%printing values
fprintf("Normalized Delta M(~CI) : Ctrl=%f, Delta T=%f, Delta E=%f, Delta T+Delta E=%f, Delta T+E0=%f, Delta T+T0=%f\n", 0,CITbar,CIEbar,CITEbar,CITE0bar,CITT0bar);
fprintf("Normalized M(~speed) : Ctrl=%f, Delta T=%f, Delta E=%f, Delta T+Delta E=%f, Delta T+E0=%f, Delta T+T0=%f\n", vCbar,vTbar,vEbar,vTEbar,vTE0bar,vTT0bar);

prp=[.5 0 .5];
   brn=[0.5 0 0];
   lgy=[0.75 0.75 0.75];
    



%plotting \Delta M
hold off;
figure(1); clf;
fig1 =bar(0,0,'FaceColor','r','LineWidth',2); hold on; 
bar(1,CITbar,'FaceColor','r','LineWidth',2); hold on; 
bar(2,CIEbar,'FaceColor','b','LineWidth',2); hold on; 
bar(3,CITEbar,'FaceColor',prp,'LineWidth',2); 
ylim([0 1.2*CITbar]);
xlim([-1 4]);
xticks([0 1 2 3])
xticklabels({'Ctrl','\nabla T','\nabla E','\nabla T+\nabla E'})
ylabel('Normalized \Delta M')
title('Mutual Repression')
xtickangle(45)
set(gca,'fontsize',20);
set(gca,'fontweight','bold','fontname','arial')  

%plotting M
hold off;
figure(2); clf;
fig1 =bar(0,vCbar,'FaceColor',lgy,'LineWidth',2); hold on; 
bar(1,vTbar,'FaceColor','r','LineWidth',2); hold on; 
bar(2,vEbar,'FaceColor','b','LineWidth',2); hold on; 
bar(3,vTEbar,'FaceColor',prp,'LineWidth',2); 
ylim([0 1.2*vTbar]);
xlim([-1 4]);
xticks([0 1 2 3])
xticklabels({'Ctrl','\nabla T','\nabla E','\nabla T+\nabla E'})
ylabel('Normalized M')
title('Mutual Repression')
xtickangle(45)
set(gca,'fontsize',20);
set(gca,'fontweight','bold','fontname','arial') 

%%
%shared pathway with background production
clear all;
clc;


maxrange=3; 
npts=1e7;
k=1;

parammax=1e3;
for i=1:npts

    
    
    rho=10^(3*(-1+2*rand));
    alpha=10^(3*(-1+2*rand));
    Theta=rand;
    
    
    etaT=alpha/(1+rho);
    etaE=alpha*rho/(1+rho);
    

SIGT=etaT/(1+etaT)^2;
   SIGEBYT=etaE/(1+etaE)^2;
   SIGTEBYT=(etaT+etaE)/(1+etaT+etaE)^2;
   SIGTE0BYT=etaT/(1+etaT+etaE)^2;
   SIGTT0BYT=etaT/(1+10*etaT)^2;
   
   vctrl=Theta;
   vnabT=(Theta+etaT)/(1+etaT);
   vEBYT=(Theta+etaE)/(1+etaE);
   vTEBYT=(Theta+etaT+etaE)/(1+etaT+etaE);
   vTE0BYT=vTEBYT;
   vTT0BYT=(Theta+10*etaT)/(1+10*etaT);


CIEfull(i)=vEBYT;
CITfull(i)=SIGT;

vEfull(i)=vEBYT;
vTfull(i)=vnabT;

if SIGT>0 && SIGEBYT>0 && SIGTEBYT<min(SIGT,SIGEBYT) 
    etx(k)=rho;
    ety(k)=alpha;
    CIT(k)=SIGT;
    CIE(k)=SIGEBYT;
    CITE(k)=SIGTEBYT;
    CITE0(k)=SIGTE0BYT;
    CITT0(k)=SIGTT0BYT;
    
    vC(k)=vctrl;
    vT(k)=vnabT;
    vE(k)=vEBYT;
    vTE(k)=vTEBYT;
    vTE0(k)=vTE0BYT;
    vTT0(k)=vTT0BYT;
    k=k+1;
end

end



clc;
CITbar=mean(CIT)/mean(CIT);
CIEbar=mean(CIE)/mean(CIT);
CITEbar=mean(CITE)/mean(CIT);
CITE0bar=mean(CITE0)/mean(CIT);
CITT0bar=mean(CITT0)/mean(CIT);


vTbar=mean(vT)/mean(vT);
vCbar=mean(vC)/mean(vT);
vEbar=mean(vE)/mean(vT);
vTEbar=mean(vTE)/mean(vT);
vTE0bar=mean(vTE0)/mean(vT);
vTT0bar=mean(vTT0)/mean(vT);


fprintf("Normalized Delta M(~CI) : Ctrl=%f, Delta T=%f, Delta E=%f, Delta T+Delta E=%f, Delta T+E0=%f, Delta T+T0=%f\n", 0,CITbar,CIEbar,CITEbar,CITE0bar,CITT0bar);
fprintf("Normalized M(~speed) : Ctrl=%f, Delta T=%f, Delta E=%f, Delta T+Delta E=%f, Delta T+E0=%f, Delta T+T0=%f\n", vCbar,vTbar,vEbar,vTEbar,vTE0bar,vTT0bar);

prp=[.5 0 .5];
   brn=[0.5 0 0];
   lgy=[0.75 0.75 0.75];
    
    hold off;
figure(1); clf;
fig1 =bar(0,0,'FaceColor','r','LineWidth',2); hold on; 
bar(1,CITbar,'FaceColor','r','LineWidth',2); hold on; 
bar(2,CIEbar,'FaceColor','b','LineWidth',2); hold on; 
bar(3,CITEbar,'FaceColor',prp,'LineWidth',2);   hold on; 
bar(4,CITE0bar,'FaceColor','k','LineWidth',2);   hold on; 
bar(5,CITT0bar,'FaceColor',brn,'LineWidth',2);  
xlim([-1 6]);
xticks([0 1 2 3 4 5])
title('Shared Pathway')
xticklabels({'Ctrl','\nabla T','\nabla E','\nabla T+\nabla E','\nabla T+E_0','\nabla T+T_0'})
ylabel('Normalized \Delta M')
xtickangle(45)
set(gca,'fontsize',20);
set(gca,'fontweight','bold','fontname','arial')

hold off;
figure(2); clf;
fig1 =bar(0,vCbar,'FaceColor',lgy,'LineWidth',2); hold on; 
bar(1,vTbar,'FaceColor','r','LineWidth',2); hold on; 
bar(2,vEbar,'FaceColor','b','LineWidth',2); hold on; 
bar(3,vTEbar,'FaceColor',prp,'LineWidth',2);   hold on; 
bar(4,vTE0bar,'FaceColor','k','LineWidth',2);   hold on; 
bar(5,vTT0bar,'FaceColor',brn,'LineWidth',2);  
xlim([-1 6]);
xticks([0 1 2 3 4 5])
title('Shared Pathway')
xticklabels({'Ctrl','\nabla T','\nabla E','\nabla T+\nabla E','\nabla T+E_0','\nabla T+T_0'})
ylabel('Normalized M')
xtickangle(45)
set(gca,'fontsize',20);
set(gca,'fontweight','bold','fontname','arial');
    
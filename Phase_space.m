%This code derives the phase space as a function of two parameters, \alpha
%and \rho. It divides the phase space into three regions : Antagonistic
%(purple), syngergistic (pink) and neither (light gray)
clear all;
clc;


for i=1:200
    rhoc(i)=0.01*10^(4*i/200); %range of rho values
    if rhoc(i)<1 %if rho<1 then
    alphacg(i)=(1/rhoc(i)+1)^0.5; %The upper line is (1+1/rho)^0.5
    alphacl(i)=(rhoc(i)+1)^0.5; %The lower line is (1+rho)^0.5
    else %if rho>1
    alphacl(i)=(1/rhoc(i)+1)^0.5; %The lower line is (1+1/rho)^0.5
    alphacg(i)=(rhoc(i)+1)^0.5; %The upper line is (1+rho)^0.5
    end
    
    alphacT(i)=(1/rhoc(i)+1)^0.5; %The line (1+1/rho)^0.5
    alphacE(i)=(rhoc(i)+1)^0.5; %The line (1+rho)^0.5
    alphac2(i)=1e3; %Some constant high value
end

%Defining 3 colors
c1=[0.5 0 1];
c2=[1 0.75 0.75];
lgy=[0.85 0.85 0.85];


%Plotting : values of \alpha>\alphacg are antagonistic, \alpha<alphacl are
%synergistic

hold off;
figure(3); clf;
h1=area(rhoc,alphac2,'FaceColor',c1); hold on;
h2=area(rhoc,alphacg,'FaceColor',lgy); hold on;
h3=area(rhoc,alphacl,'FaceColor',c2); hold on;
h4=plot(rhoc,alphacT,'k','LineWidth',4); hold on;
h5=plot(rhoc,alphacE,'k','LineWidth',4);
xlabel('Relative pathway strength, \rho');
ylabel('Amplification factor, \alpha');
set(gca,'fontsize',25);
set(gca, 'LineWidth', 3.5)
% Choose a number between 0 (invisible) and 1 (opaque) for facealpha.  
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlim([1e-2 1e2])
ylim([1e-1 1e2])
xticks([1e-2 1 1e2])
yticks([1e-1 1 1e1 1e2])
saveas(h1,'phase_space_V3.png')


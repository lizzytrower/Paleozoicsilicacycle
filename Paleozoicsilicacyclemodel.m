%Paleozoic silica cycle Si isotope models

%These models were designed by Lizzy Trower (CU Boulder) using Matlab 2018b
%as simple tools to estimate fluxes and Si isotope compositions (d30Si) for
%the Paleozoic marine silica cycle. The model was modified after versions
%described in De La Rocha and Bickle (2015) and Frings et al. (2016). Each
%section of the code runs a different version described in the text of
%Trower et al. (submitted manuscript). The first section is commented; the
%remaining sections are structured identically, and just differ in terms of
%the ODE math and some additional parameters.

%This version was last updated on March 23, 2020.

%%
%Run this section for modern-like diatom + sponge model. 

clear

%This solves the primary ODE for the silica fluxes.
Fin = 5*10^12; %[mol/yr]
[t,C] = ode45(@(t,C) sibox_diatomsponge(t,C,Fin),0:.1:50000,[.5; .5]);
%units for timesteps (t) are years; units for DSi concentration (C) are mM

%Plot the evolution of concentrations in surface and deep ocean reservoirs
%over time.
figure
subplot(1,3,1)
semilogy(t,C(:,1)*10^3,t,C(:,2)*10^3)
legend('surface','deep')
ylim([1 10000])
xlim([0 5*10^4])
ylabel('[DSi] \muM')

%More convenient labels for concentration vectors.
Csurf = C(:,1);
Cdeep = C(:,2);

%Set up variables; these must match the variables set up in the
%'sibox_diatomsponge' function. If you change one of these values (e.g.,
%Fin), it must be changed in both places. These could be rewritten as
%function inputs to avoid that.
volocean = 1.35*10^18; %[m^3]
SAocean = 3.6*10^14; %[m^2]
Feuphotic = 0.0263;
volsurf = volocean*Feuphotic; %[m^3]
voldeep = volocean-volsurf; %[m^3]
Wex = 1.37*10^15; %[m^3/yr]
Vmax = 500*10^12; %[mol/yr]
Km = 9; %[uM]
Km = Km/1000; %[mol/m^3]
Velsurf = 1800; %[m/yr]
Veldeep = 73000; %[m/yr]
Dsurf = volsurf/SAocean; %[m]
Ddeep = volocean/SAocean - Dsurf; %[m]
ksurf = 9; %[mol/mol/yr]
kdeep = 24; %[mol/mol/yr]
Ceq = 350; %[uM]
Ceq = Ceq*10^-3; %[mol/m^3]
Km_sponge = 75;
Km_sponge = Km_sponge/1000;
Vmax_sponge = 0.13; % [umol/hr/mL sponge]
Vmax_sponge = Vmax_sponge*10^-6*24*365*1000; %[mol/yr/L sponge]
spongeseafloordensity = 0.34; % [L sponge/m^2]
Vmax_sponge = Vmax_sponge*spongeseafloordensity*SAocean/50; %[mol/yr]

%Calculate BSi productivity by diatoms and sponges.
BSi_prod_diatom = Vmax*Csurf./(Km + Csurf);
f_BSi_diatom = 1 - BSi_prod_diatom./(Csurf.*volsurf);
BSi_prod_sponge = Vmax_sponge*Cdeep./(Km_sponge + Cdeep);
f_BSi_sponge = 1 - BSi_prod_sponge./(Cdeep.*voldeep);

%Set d30Si of input (Fin) and eps30Si for diatoms. eps30Si for sponges is
%set by a function (here, 'spongefrac', which is based on the modern
%calibration) that depends on [DSi].
d30Si_in = 1;
eps_diatom = -1.1;

%Set up vectors for d30Si outputs.
d30Sisurf = zeros(length(t),1);
d30Sideep = zeros(length(t),1);
d30BSi_diatom = zeros(length(t),1);
d30BSi_sponge = zeros(length(t),1);
burialBSi = zeros(length(t),1);

%Run for loop to calculate d30Si values for each of the components for each
%timestep.
for n = 2:length(t)
    %Diatom and sponge d30Si are determined using a Rayleigh model
    d30BSi_diatom(n) = d30Sisurf(n-1) - eps_diatom*log(f_BSi_diatom(n-1))*...
        f_BSi_diatom(n-1)/(1-f_BSi_diatom(n-1));
    d30BSi_sponge(n) = d30Sideep(n-1) - spongefrac(Cdeep(n-1)*10^3)*...
        log(f_BSi_sponge(n-1))*f_BSi_sponge(n-1)/(1-f_BSi_sponge(n-1));
    d30Sisurf(n) = d30Sisurf(n-1) + 1/(volsurf*Csurf(n-1))*(Fin*d30Si_in ...
        - Wex*Csurf(n-1)*d30Sisurf(n-1) + Wex*Cdeep(n-1)*d30Sideep(n-1) ...
        - (BSi_prod_diatom(n-1) - max([ksurf*(1-Csurf(n-1)/Ceq),0]).*BSi_prod_diatom(n-1)*...
        Dsurf/Velsurf)*d30BSi_diatom(n-1));
    d30Sideep(n) = d30Sideep(n-1) + 1/(voldeep*Cdeep(n-1))*(Wex*Csurf(n-1)*...
        d30Sisurf(n-1) - Wex*Cdeep(n-1)*d30Sideep(n-1) + ...
        max([kdeep*(1-Cdeep(n-1)/Ceq),0]).*(BSi_prod_diatom(n-1) - ...
        max([ksurf*(1-Csurf(n-1)/Ceq),0]).*BSi_prod_diatom(n-1)*...
        Dsurf/Velsurf)*Ddeep/Veldeep*...
        d30BSi_diatom(n-1) - BSi_prod_sponge(n-1)*d30BSi_sponge(n-1));
    burialBSi(n) = BSi_prod_diatom(n-1) - max([ksurf*(1-Csurf(n-1)/Ceq),0]).*...
        BSi_prod_diatom(n-1)*Dsurf/Velsurf - max([kdeep*(1-Cdeep(n-1)/Ceq),0]).*...
        (BSi_prod_diatom(n-1) - max([ksurf*(1-Csurf(n-1)/Ceq),0]).*BSi_prod_diatom(n-1)*...
        Dsurf/Velsurf)*Ddeep/Veldeep;
end

%Add subplots of d30Si and BSi fluxes over time.
subplot(1,3,2)
plot(t,d30Sisurf,t,d30Sideep,t,d30BSi_diatom,t,d30BSi_sponge)
%Add horizontal line to indicate input (d30Si of Fin)
yline(d30Si_in)
legend('surface','deep','diatom','sponge','input')
ylabel('\delta^3^0Si')
ylim([-5 6])
xlim([0 5*10^4])

subplot(1,3,3)
semilogy(t,burialBSi,t,BSi_prod_sponge,t,burialBSi+BSi_prod_sponge)
%Add horizontal line to indicate input flux (Fin)
yline(Fin)
ylim([10^12 10^15])
xlim([0 5*10^4])
ylabel('BSi output (mol/yr)')
legend('diatom','sponge','sum','input')

%%
%Run this section for sponge-only model (uses modern DSi affinities and
%modern eps30Si-[DSi] calibration

clear

Fin = 5*10^12; %[mol/yr]
[t,C] = ode45(@(t,C) sibox_sponge(t,C,Fin),0:.1:200000,[.5; .5]);

figure
subplot(1,3,1)
semilogy(t,C(:,1)*10^3,t,C(:,2)*10^3)
legend('surface','deep')
ylim([1 10000])
ylabel('[DSi] \muM')

Csurf = C(:,1);
Cdeep = C(:,2);

volocean = 1.35*10^18; %[m^3]
SAocean = 3.6*10^14; %[m^2]
Feuphotic = 0.0263;
volsurf = volocean*Feuphotic; %[m^3]
voldeep = volocean-volsurf; %[m^3]
Wex = 1.37*10^15; %[m^3/yr]
Vmax = 0.13; % [umol/hr/mL sponge]
Vmax = Vmax*10^-6*24*365*1000; %[mol/yr/L sponge]
spongeseafloordensity = 0.34; % [L sponge/m^2]
Vmax = Vmax*spongeseafloordensity*SAocean/5; %[mol/yr]
Km = 75; %[uM]
Km = Km/1000; %[mol/m^3]

BSi_prod_sponge = Vmax*Cdeep./(Km + Cdeep);
f_BSi_sponge = 1 - BSi_prod_sponge./(Cdeep.*voldeep);

d30Si_in = 1;

d30Sisurf = zeros(length(t),1);
d30Sisurf(1) = 1;
d30Sideep = zeros(length(t),1);
d30Sideep(1) = 1;
d30BSi = zeros(length(t),1);
d30BSi(1) = d30Sideep(1) + spongefrac(Cdeep(1));

for n = 2:length(t)
    d30BSi(n) = d30Sideep(n-1) - spongefrac(Cdeep(n-1)*10^3)*...
        log(f_BSi_sponge(n-1))*f_BSi_sponge(n-1)/(1-f_BSi_sponge(n-1));
    d30Sisurf(n) = d30Sisurf(n-1) + 1/(volsurf*Csurf(n-1))*(Fin*d30Si_in ...
        - Wex*Csurf(n-1)*d30Sisurf(n-1) + Wex*Cdeep(n-1)*d30Sideep(n-1));
    d30Sideep(n) = d30Sideep(n-1) + 1/(voldeep*Cdeep(n-1))*(Wex*Csurf(n-1)*...
        d30Sisurf(n-1) - Wex*Cdeep(n-1)*d30Sideep(n-1) - ...
        BSi_prod_sponge(n-1)*d30BSi(n-1));

end
    
subplot(1,3,2)
plot(t,d30Sisurf,t,d30Sideep,t,d30BSi)
yline(d30Si_in)
legend('surface','deep','sponge','input')
ylabel('\delta^3^0Si')
ylim([-5 6])

subplot(1,3,3)
semilogy(t,BSi_prod_sponge)
yline(Fin)
ylim([10^12,10^15])
ylabel('BSi output (mol/yr)')
legend('sponge','input')

%%
%Run this section for sponge + radiolarian model w/ modern DSi affinities
%and modern eps30Si-[DSi] calibration

clear

Fin = 5*10^12; %[mol/yr]
[t,C] = ode45(@(t,C) sibox_spongerad(t,C,Fin),0:.1:100000,[.5; .5]);

figure
subplot(1,3,1)
semilogy(t,C(:,1)*10^3,t,C(:,2)*10^3)
legend('surface','deep')
ylim([1 10000])
ylabel('[DSi] \muM')

Csurf = C(:,1);
Cdeep = C(:,2);

volocean = 1.35*10^18; %[m^3]
SAocean = 3.6*10^14; %[m^2]
Feuphotic = 0.0263;
volsurf = volocean*Feuphotic; %[m^3]
voldeep = volocean-volsurf; %[m^3]
Wex = 1.37*10^15; %[m^3/yr]
Vmax_sponge = 0.13; % [umol/hr/mL sponge]
Vmax_sponge = Vmax_sponge*10^-6*24*365*1000; %[mol/yr/L sponge]
spongeseafloordensity = 0.34; % [L sponge/m^2]
Vmax_sponge = Vmax_sponge*spongeseafloordensity*SAocean/5; %[mol/yr]
Km_sponge = 75; %[uM]
Km_sponge = Km_sponge/1000; %[mol/m^3]
Vmax_rad = Vmax_sponge;
Km_rad = Km_sponge;
Velsurf = 1800; %[m/yr]
Veldeep = 73000; %[m/yr]
Dsurf = volsurf/SAocean; %[m]
Ddeep = volocean/SAocean - Dsurf; %[m]
ksurf = 9; %[mol/mol/yr]
kdeep = 3; %[mol/mol/yr]
Ceq = 100; %[uM]
Ceq = Ceq*10^-3; %[mol/m^3]

BSi_prod_rad = Vmax_rad*Csurf./(Km_rad + Csurf);
BSi_prod_sponge = Vmax_sponge*Cdeep./(Km_sponge + Cdeep);
f_BSi_rad = 1 - BSi_prod_rad./(Csurf.*volsurf);
f_BSi_sponge = 1 - BSi_prod_sponge./(Cdeep.*voldeep);

d30Si_in = 1;
eps_BSi_rad = -1.1;

d30Sisurf = zeros(length(t),1);
d30Sideep = zeros(length(t),1);
d30BSi_rad = zeros(length(t),1);
d30BSi_sponge = zeros(length(t),1);
burialBSi_rad = zeros(length(t),1);

for n = 2:length(t)
    d30BSi_rad(n) = d30Sisurf(n-1) - eps_BSi_rad*log(f_BSi_rad(n-1))*...
        f_BSi_rad(n-1)/(1-f_BSi_rad(n-1));
    d30BSi_sponge(n) = d30Sideep(n-1) - spongefrac(Cdeep(n-1)*10^3)*...
        log(f_BSi_sponge(n-1))*f_BSi_sponge(n-1)/(1-f_BSi_sponge(n-1));
    d30Sisurf(n) = d30Sisurf(n-1) + 1/(volsurf*Csurf(n-1))*(Fin*d30Si_in ...
        - Wex*Csurf(n-1)*d30Sisurf(n-1) + Wex*Cdeep(n-1)*d30Sideep(n-1) ...
        - (BSi_prod_rad(n-1) - max([ksurf*(1-Csurf(n-1)/Ceq),0]).*BSi_prod_rad(n-1)*...
        Dsurf/Velsurf)*d30BSi_rad(n-1));
    d30Sideep(n) = d30Sideep(n-1) + 1/(voldeep*Cdeep(n-1))*(Wex*Csurf(n-1)*...
        d30Sisurf(n-1) - Wex*Cdeep(n-1)*d30Sideep(n-1) + ...
        max([kdeep*(1-Cdeep(n-1)/Ceq),0]).*(BSi_prod_rad(n-1) - ...
        max([ksurf*(1-Csurf(n-1)/Ceq),0]).*BSi_prod_rad(n-1)*...
        Dsurf/Velsurf)*Ddeep/Veldeep*...
        d30BSi_rad(n-1) - BSi_prod_sponge(n-1)*d30BSi_sponge(n-1));
    burialBSi_rad(n) = BSi_prod_rad(n-1) - max([ksurf*(1-Csurf(n-1)/Ceq),0]).*...
        BSi_prod_rad(n-1)*Dsurf/Velsurf - max([kdeep*(1-Cdeep(n-1)/Ceq),0]).*...
        (BSi_prod_rad(n-1) - max([ksurf*(1-Csurf(n-1)/Ceq),0]).*BSi_prod_rad(n-1)*...
        Dsurf/Velsurf)*Ddeep/Veldeep;

end

subplot(1,3,2)
plot(t,d30Sisurf,t,d30Sideep,t,d30BSi_rad,t,d30BSi_sponge)
yline(d30Si_in)
legend('surface','deep','radiolarian','sponge','input')
ylabel('\delta^3^0Si')
ylim([-5 6])

subplot(1,3,3)
semilogy(t,burialBSi_rad,t,BSi_prod_sponge,t,burialBSi_rad+BSi_prod_sponge)
yline(Fin)
ylim([10^12,10^15])
ylabel('BSi output (mol/yr)')
legend('radiolarian','sponge','sum','input')

%%
%Run this section for sponge + radiolarian model with low DSi affinity and
%modern eps30Si-[DSi] calibration

clear

Fin = 5*10^12; %[mol/yr]
[t,C] = ode45(@(t,C)sibox_spongerad_lowDSiaffinity(t,C,Fin),0:.1:1500000,[.5; .5]);

figure
subplot(1,3,1)
semilogy(t,C(:,1)*10^3,t,C(:,2)*10^3)
legend('surface','deep')
ylim([1 10000])
ylabel('[DSi] \muM')

Csurf = C(:,1);
Cdeep = C(:,2);

volocean = 1.35*10^18; %[m^3]
SAocean = 3.6*10^14; %[m^2]
Feuphotic = 0.0263;
volsurf = volocean*Feuphotic; %[m^3]
voldeep = volocean-volsurf; %[m^3]
Wex = 1.37*10^15; %[m^3/yr]
Vmax_sponge = 0.05; % [umol/hr/mL sponge]
Vmax_sponge = Vmax_sponge*10^-6*24*365*1000; %[mol/yr/L sponge]
spongeseafloordensity = 0.34; % [L sponge/m^2]
Vmax_sponge = Vmax_sponge*spongeseafloordensity*SAocean/5; %[mol/yr]
Km_sponge = 500; %[uM]
Km_sponge = Km_sponge/1000; %[mol/m^3]
Vmax_rad = Vmax_sponge;
Km_rad = Km_sponge;
Velsurf = 1800; %[m/yr]
Veldeep = 73000; %[m/yr]
Dsurf = volsurf/SAocean; %[m]
Ddeep = volocean/SAocean - Dsurf; %[m]
ksurf = 9; %[mol/mol/yr]
kdeep = 3; %[mol/mol/yr]
Ceq = 350; %[uM]
Ceq = Ceq*10^-3; %[mol/m^3]

BSi_prod_rad = Vmax_rad*Csurf./(Km_rad + Csurf);
BSi_prod_sponge = Vmax_sponge*Cdeep./(Km_sponge + Cdeep);
f_BSi_rad = 1 - BSi_prod_rad./(Csurf.*volsurf);
f_BSi_sponge = 1 - BSi_prod_sponge./(Cdeep.*voldeep);

d30Si_in = 1;
eps_BSi_rad = -1.1;

d30Sisurf = zeros(length(t),1);
d30Sideep = zeros(length(t),1);
d30BSi_rad = zeros(length(t),1);
d30BSi_sponge = zeros(length(t),1);
burialBSi_rad = zeros(length(t),1);

for n = 2:length(t)
    d30BSi_rad(n) = d30Sisurf(n-1) - eps_BSi_rad*log(f_BSi_rad(n-1))*...
        f_BSi_rad(n-1)/(1-f_BSi_rad(n-1));
    d30BSi_sponge(n) = d30Sideep(n-1) - spongefrac(Cdeep(n-1)*10^3)*...
        log(f_BSi_sponge(n-1))*f_BSi_sponge(n-1)/(1-f_BSi_sponge(n-1));
    d30Sisurf(n) = d30Sisurf(n-1) + 1/(volsurf*Csurf(n-1))*(Fin*d30Si_in ...
        - Wex*Csurf(n-1)*d30Sisurf(n-1) + Wex*Cdeep(n-1)*d30Sideep(n-1) ...
        - (BSi_prod_rad(n-1) - max([ksurf*(1-Csurf(n-1)/Ceq),0]).*BSi_prod_rad(n-1)*...
        Dsurf/Velsurf)*d30BSi_rad(n-1));
    d30Sideep(n) = d30Sideep(n-1) + 1/(voldeep*Cdeep(n-1))*(Wex*Csurf(n-1)*...
        d30Sisurf(n-1) - Wex*Cdeep(n-1)*d30Sideep(n-1) + ...
        max([kdeep*(1-Cdeep(n-1)/Ceq),0]).*(BSi_prod_rad(n-1) - ...
        max([ksurf*(1-Csurf(n-1)/Ceq),0]).*BSi_prod_rad(n-1)*...
        Dsurf/Velsurf)*Ddeep/Veldeep*...
        d30BSi_rad(n-1) - BSi_prod_sponge(n-1)*d30BSi_sponge(n-1));
    burialBSi_rad(n) = BSi_prod_rad(n-1) - max([ksurf*(1-Csurf(n-1)/Ceq),0]).*...
        BSi_prod_rad(n-1)*Dsurf/Velsurf - max([kdeep*(1-Cdeep(n-1)/Ceq),0]).*...
        (BSi_prod_rad(n-1) - max([ksurf*(1-Csurf(n-1)/Ceq),0]).*BSi_prod_rad(n-1)*...
        Dsurf/Velsurf)*Ddeep/Veldeep;

end

subplot(1,3,2)
plot(t,d30Sisurf,t,d30Sideep,t,d30BSi_rad,t,d30BSi_sponge)
yline(d30Si_in)
legend('surface','deep','radiolarian','sponge','input')
ylabel('\delta^3^0Si')
ylim([-5 6])

subplot(1,3,3)
semilogy(t,burialBSi_rad,t,BSi_prod_sponge,t,burialBSi_rad+...
    BSi_prod_sponge)
yline(Fin)
ylim([10^12,10^15])
ylabel('BSi output (mol/yr)')
legend('radiolarian','sponge','sum','input')

%%
%Run this section for sponge + radiolarian model with low DSi affinity and
%Michaelis-Menten eps30Si-[DSi] calibration

clear

Fin = 5*10^12; %[mol/yr]
[t,C] = ode45(@(t,C)sibox_spongerad_lowDSiaffinity(t,C,Fin),0:.1:1500000,[.5; .5]);

figure
subplot(1,3,1)
semilogy(t,C(:,1)*10^3,t,C(:,2)*10^3)
legend('surface','deep')
ylim([1 10000])
ylabel('[DSi] \muM')

Csurf = C(:,1);
Cdeep = C(:,2);

volocean = 1.35*10^18; %[m^3]
SAocean = 3.6*10^14; %[m^2]
Feuphotic = 0.0263;
volsurf = volocean*Feuphotic; %[m^3]
voldeep = volocean-volsurf; %[m^3]
Wex = 1.37*10^15; %[m^3/yr]
Vmax_sponge = 0.05; % [umol/hr/mL sponge]
Vmax_sponge = Vmax_sponge*10^-6*24*365*1000; %[mol/yr/L sponge]
spongeseafloordensity = 0.34; % [L sponge/m^2]
Vmax_sponge = Vmax_sponge*spongeseafloordensity*SAocean/5; %[mol/yr]
Km_sponge = 500; %[uM]
Km_sponge = Km_sponge/1000; %[mol/m^3]
Vmax_rad = Vmax_sponge;
Km_rad = Km_sponge;
Velsurf = 1800; %[m/yr]
Veldeep = 73000; %[m/yr]
Dsurf = volsurf/SAocean; %[m]
Ddeep = volocean/SAocean - Dsurf; %[m]
ksurf = 9; %[mol/mol/yr]
kdeep = 3; %[mol/mol/yr]
Ceq = 350; %[uM]
Ceq = Ceq*10^-3; %[mol/m^3]

BSi_prod_rad = Vmax_rad*Csurf./(Km_rad + Csurf);
BSi_prod_sponge = Vmax_sponge*Cdeep./(Km_sponge + Cdeep);
f_BSi_rad = 1 - BSi_prod_rad./(Csurf.*volsurf);
f_BSi_sponge = 1 - BSi_prod_sponge./(Cdeep.*voldeep);

d30Si_in = 1;
eps_BSi_rad = -1.1;

d30Sisurf = zeros(length(t),1);
d30Sideep = zeros(length(t),1);
d30BSi_rad = zeros(length(t),1);
d30BSi_sponge = zeros(length(t),1);
burialBSi_rad = zeros(length(t),1);

for n = 2:length(t)
    d30BSi_rad(n) = d30Sisurf(n-1) - eps_BSi_rad*log(f_BSi_rad(n-1))*...
        f_BSi_rad(n-1)/(1-f_BSi_rad(n-1));
    d30BSi_sponge(n) = d30Sideep(n-1) - MMspongefrac(Cdeep(n-1)*10^3)*...
        log(f_BSi_sponge(n-1))*f_BSi_sponge(n-1)/(1-f_BSi_sponge(n-1));
    d30Sisurf(n) = d30Sisurf(n-1) + 1/(volsurf*Csurf(n-1))*(Fin*d30Si_in ...
        - Wex*Csurf(n-1)*d30Sisurf(n-1) + Wex*Cdeep(n-1)*d30Sideep(n-1) ...
        - (BSi_prod_rad(n-1) - max([ksurf*(1-Csurf(n-1)/Ceq),0]).*BSi_prod_rad(n-1)*...
        Dsurf/Velsurf)*d30BSi_rad(n-1));
    d30Sideep(n) = d30Sideep(n-1) + 1/(voldeep*Cdeep(n-1))*(Wex*Csurf(n-1)*...
        d30Sisurf(n-1) - Wex*Cdeep(n-1)*d30Sideep(n-1) + ...
        max([kdeep*(1-Cdeep(n-1)/Ceq),0]).*(BSi_prod_rad(n-1) - ...
        max([ksurf*(1-Csurf(n-1)/Ceq),0]).*BSi_prod_rad(n-1)*...
        Dsurf/Velsurf)*Ddeep/Veldeep*...
        d30BSi_rad(n-1) - BSi_prod_sponge(n-1)*d30BSi_sponge(n-1));
    burialBSi_rad(n) = BSi_prod_rad(n-1) - max([ksurf*(1-Csurf(n-1)/Ceq),0]).*...
        BSi_prod_rad(n-1)*Dsurf/Velsurf - max([kdeep*(1-Cdeep(n-1)/Ceq),0]).*...
        (BSi_prod_rad(n-1) - max([ksurf*(1-Csurf(n-1)/Ceq),0]).*BSi_prod_rad(n-1)*...
        Dsurf/Velsurf)*Ddeep/Veldeep;

end

subplot(1,3,2)
plot(t,d30Sisurf,t,d30Sideep,t,d30BSi_rad,t,d30BSi_sponge)
yline(d30Si_in)
legend('surface','deep','radiolarian','sponge','input')
ylabel('\delta^3^0Si')
ylim([-5 6])

subplot(1,3,3)
semilogy(t,burialBSi_rad,t,BSi_prod_sponge,t,burialBSi_rad+...
    BSi_prod_sponge)
yline(Fin)
ylim([10^12,10^15])
ylabel('BSi output (mol/yr)')
legend('radiolarian','sponge','sum','input')

%%
function dCdt =sibox_diatomsponge(t,C,Fin)

    volocean = 1.35*10^18; %[m^3]
    SAocean = 3.6*10^14; %[m^2]
    Feuphotic = 0.0263;
    volsurf = volocean*Feuphotic; %[m^3]
    voldeep = volocean-volsurf; %[m^3]
    Wex = 1.37*10^15; %[m^3/yr]
    Vmax = 500*10^12; %[mol/yr]
    Km = 9; %[uM]
    Km = Km/1000; %[mol/m^3 = mM]
    Velsurf = 1800; %[m/yr]
    Veldeep = 73000; %[m/yr]
    Dsurf = volsurf/SAocean; %[m]
    Ddeep = volocean/SAocean - Dsurf; %[m]
    ksurf = 9; %[mol/mol/yr]
    kdeep = 24; %[mol/mol/yr]
    Ceq = 350; %[uM]
    Ceq = Ceq*10^-3; %[mol/m^3 = mM]
    Km_sponge = 75;
    Km_sponge = Km_sponge/1000;
    Vmax_sponge = 0.13; % [umol/hr/mL sponge]
    Vmax_sponge = Vmax_sponge*10^-6*24*365*1000; %[mol/yr/L sponge]
    spongeseafloordensity = 0.34; % [L sponge/m^2]
    Vmax_sponge = Vmax_sponge*spongeseafloordensity*SAocean/50; %[mol/yr]
    
    dCdt = [1/volsurf*(Fin - Vmax*C(1)/(Km + C(1)) + ...
        max([ksurf*(1-C(1)/Ceq),0])*Vmax*C(1)/(Km + C(1))...
        *Dsurf/Velsurf - Wex*C(1) + Wex*C(2));...
        1/voldeep*(Wex*C(1) - Wex*C(2) + ...
        max([kdeep*(1-C(2)/Ceq),0])*(Vmax*C(1)/(Km + C(1)) - ...
        max([ksurf*(1-C(1)/Ceq),0])*...
        Vmax*C(1)/(Km + C(1))*Dsurf/Velsurf)...
        *Ddeep/Veldeep - Vmax_sponge*C(2)/(Km_sponge + C(2)))];
    
end

function dCdt =sibox_sponge(t,C,Fin)

    volocean = 1.35*10^18; %[m^3]
    SAocean = 3.6*10^14; %[m^2]
    Feuphotic = 0.0263;
    volsurf = volocean*Feuphotic; %[m^3]
    voldeep = volocean-volsurf; %[m^3]
    Wex = 1.37*10^15; %[m^3/yr]
    Vmax = 0.13; % [umol/hr/mL sponge]
    Vmax = Vmax*10^-6*24*365*1000; %[mol/yr/L sponge]
    spongeseafloordensity = 0.34; % [L sponge/m^2]
    Vmax = Vmax*spongeseafloordensity*SAocean/5; %[mol/yr]
    Km = 75; %[uM]
    Km = Km/1000; %[mol/m^3]
    
    dCdt = [1/volsurf*(Fin - Wex*C(1) + Wex*C(2));...
        1/voldeep*(Wex*C(1) - Wex*C(2) - Vmax*C(2)/(Km + C(2)))];
    
end

function dCdt =sibox_spongerad(t,C,Fin)

    volocean = 1.35*10^18; %[m^3]
    SAocean = 3.6*10^14; %[m^2]
    Feuphotic = 0.0263;
    volsurf = volocean*Feuphotic; %[m^3]
    voldeep = volocean-volsurf; %[m^3]
    Wex = 1.37*10^15; %[m^3/yr]
    Velsurf = 1800; %[m/yr]
    Veldeep = 73000; %[m/yr]
    Dsurf = volsurf/SAocean; %[m]
    Ddeep = volocean/SAocean - Dsurf; %[m]
    ksurf = 9; %[mol/mol/yr]
    kdeep = 3; %[mol/mol/yr]
    Ceq = 100; %[uM]
    Ceq = Ceq*10^-3; %[mol/m^3]
    Vmax_sponge = 0.13; % [umol/hr/mL sponge]
    Vmax_sponge = Vmax_sponge*10^-6*24*365*1000; %[mol/yr/L sponge]
    spongeseafloordensity = 0.34; % [L sponge/m^2]
    Vmax_sponge = Vmax_sponge*spongeseafloordensity*SAocean/5; %[mol/yr]
    Km_sponge = 75; %[uM]
    Km_sponge = Km_sponge/1000; %[mol/m^3]
    Vmax_rad = Vmax_sponge;
    Km_rad = Km_sponge;
    
    dCdt = [1/volsurf*(Fin - Vmax_rad*C(1)/(Km_rad + C(1)) + ...
        max([ksurf*(1-C(1)/Ceq),0])*Vmax_rad*C(1)/(Km_rad + C(1))...
        *Dsurf/Velsurf - Wex*C(1) + Wex*C(2));...
        1/voldeep*(Wex*C(1) - Wex*C(2) + ...
        max([kdeep*(1-C(2)/Ceq),0])*(Vmax_rad*C(1)/(Km_rad + C(1)) - ...
        max([ksurf*(1-C(1)/Ceq),0])*...
        Vmax_rad*C(1)/(Km_rad + C(1))*Dsurf/Velsurf)...
        *Ddeep/Veldeep - Vmax_sponge*C(2)/(Km_sponge + C(2)))];
    
end

function dCdt =sibox_spongerad_lowDSiaffinity(t,C,Fin)

    volocean = 1.35*10^18; %[m^3]
    SAocean = 3.6*10^14; %[m^2]
    Feuphotic = 0.0263;
    volsurf = volocean*Feuphotic; %[m^3]
    voldeep = volocean-volsurf; %[m^3]
    Wex = 1.37*10^15; %[m^3/yr]
    Velsurf = 1800; %[m/yr]
    Veldeep = 73000; %[m/yr]
    Dsurf = volsurf/SAocean; %[m]
    Ddeep = volocean/SAocean - Dsurf; %[m]
    ksurf = 9; %[mol/mol/yr]
    kdeep = 3; %[mol/mol/yr]
    Ceq = 350; %[uM]
    Ceq = Ceq*10^-3; %[mol/m^3]
    Vmax_sponge = 0.05; % [umol/hr/mL sponge]
    Vmax_sponge = Vmax_sponge*10^-6*24*365*1000; %[mol/yr/L sponge]
    spongeseafloordensity = 0.34; % [L sponge/m^2]
    Vmax_sponge = Vmax_sponge*spongeseafloordensity*SAocean/5; %[mol/yr]
    Km_sponge = 500; %[uM]
    Km_sponge = Km_sponge/1000; %[mol/m^3]
    Vmax_rad = Vmax_sponge;
    Km_rad = Km_sponge;
    
    dCdt = [1/volsurf*(Fin - Vmax_rad*C(1)/(Km_rad + C(1)) + ...
        max([ksurf*(1-C(1)/Ceq),0])*Vmax_rad*C(1)/(Km_rad + C(1))...
        *Dsurf/Velsurf - Wex*C(1) + Wex*C(2));...
        1/voldeep*(Wex*C(1) - Wex*C(2) + ...
        max([kdeep*(1-C(2)/Ceq),0])*(Vmax_rad*C(1)/(Km_rad + C(1)) - ...
        max([ksurf*(1-C(1)/Ceq),0])*...
        Vmax_rad*C(1)/(Km_rad + C(1))*Dsurf/Velsurf)...
        *Ddeep/Veldeep - Vmax_sponge*C(2)/(Km_sponge + C(2)))];
    
end

function eps_f = MMspongefrac(in)
    eps_tI = -1.34;
    eps_PE = -5.39;
    Vmax_P = 0.05; % [umol/hr/mL sponge]
    Vmax_P = Vmax_P*10^-6*24*365*1000; %[mol/yr/L sponge]
    SAocean = 3.6*10^14; %[m^2]
    spongeseafloordensity = 0.34; % [L sponge/m^2]
    Vmax_P = Vmax_P*spongeseafloordensity*SAocean/5; %[mol/yr]
    Km_P = 500; %[uM]
    Vmax_I = 9;% [umol/hr/mL sponge]
    Vmax_I = Vmax_I*10^-6*24*365*1000; %[mol/yr/L sponge]
    Vmax_I = Vmax_I*spongeseafloordensity*SAocean/5; %[mol/yr]
    Km_I = Vmax_I*Km_P/Vmax_P;
    eps_f = eps_tI + eps_PE*(1-Vmax_P./(Km_P./in+1)./(Vmax_I./(Km_I./in+1)));
end

function out = spongefrac(in)
    out = -4.6 + 27.6./(7.4+in);
end
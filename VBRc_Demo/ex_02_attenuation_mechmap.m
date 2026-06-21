
addpath(getenv('VBRpath'))
close all; clear all; clc
vbr_init


%set temperature
T = 1400+273; %K, temperature

%set frequency
f =  1e-2; %Hz

%  build 2D grid of grain size and stress
d   = logspace(1,    5,    40)      ; %um, grain sizes
sig = logspace(-1,   3,    41)      ; %MPa, differential stress

[SV.dg_um, SV.sig_MPa] = meshgrid(d, sig); %creates an array with dimensions in order of d, and then sig

sz = size(SV.dg_um);

% set constants
SV.T_K   = T * ones(sz); %K, Temperature
SV.f = f;
SV.phi = 0 * ones(sz);
SV.P_GPa = 3.5 * ones(sz); % pressure [GPa]
SV.Tsolidus_K = -5.104.*SV.P_GPa.^2 + 132.899.*SV.P_GPa + 1120.661 +273; %K , solidus temperature (for premelt model) from Hirschmann, 2000
SV.rho = 3300 * ones(sz); %kgm-3, density

VBR = struct();
VBR.in.SV = SV;
VBR.in.elastic.methods_list = {'anharmonic'}; % set methods list
VBR.in.elastic.anharmonic = Params_Elastic('anharmonic');
VBR.in.elastic.anharmonic.temperature_scaling = 'isaak';
VBR.in.elastic.anharmonic.pressure_scaling = 'abramson';


%% run VBR for each model
% Run VBR with linearized backstress model and save
disp("Computing linearized backstress model")
VBR_linBackstress = struct();
VBR_linBackstress.in = VBR.in;
VBR_linBackstress.in.anelastic.methods_list = {'backstress_linear'};
[VBR_linBackstress] = VBR_spine(VBR_linBackstress); % run VBR

% Run VBR with peak of premelt model only
disp("Computing xfit_premelt with no HTB")
VBR_peak = struct();
VBR_peak.in = VBR.in;
VBR_peak.in.anelastic.methods_list = {'xfit_premelt'};
VBR_peak.in.anelastic.xfit_premelt.A_B=0; % set amplitude of HTB to zero
[VBR_peak] = VBR_spine(VBR_peak); % run VBR

% Run VBR with HTB of premelt model only
disp("Computing xfit_premelt with no high frequency peak")
VBR_HTB = struct();
VBR_HTB.in = VBR.in;
VBR_HTB.in.anelastic.methods_list = {'xfit_premelt'};
VBR_HTB.in.anelastic.xfit_premelt.Ap_fac_1=0; % set amplitudes of high-frequency peak to zero
VBR_HTB.in.anelastic.xfit_premelt.Ap_fac_2=0;
VBR_HTB.in.anelastic.xfit_premelt.Ap_fac_3=0;
[VBR_HTB] = VBR_spine(VBR_HTB); % run VBR

% Run VBR with steady-state dislocation creep only
disp("Analytical maxwell with viscous backstress viscosity")
VBR_disl = struct();
VBR_disl.in = VBR.in;
VBR_disl.in.viscous.methods_list={'BKHK2023'};
VBR_disl.in.anelastic.methods_list={'maxwell_analytical'};
% select viscous method, in this case the backstress model with
% dislocation recovery by grain-boundary and pipe diffusion
VBR_disl.in.anelastic.maxwell_analytical.viscosity_method_mechanism = 'gbnp';
[VBR_disl] = VBR_spine(VBR_disl); % run VBR

% save VBR results for later use (optional)
% save('VBR_linBackstress_pmHTB_pmPeak_61','VBR_HTB',"VBR_peak","VBR_linBackstress","VBR_disl","VBR_in",'-v7.3');

%% or alternatively, load VBR results
% load VBR_XXX.mat

%% calculate complex viscosity and effective viscosity for each model

VBR.in.SV.f_reshape = reshape(VBR.in.SV.f, [1 1 1 length(f)]); %Hz, reshpae vector for multiplication with 4-D array

disp("Calculating complex viscosity, effective viscosity")
%complex viscosity
VBR_linBackstress.out.anelastic.backstress_linear.eta_complex = -1i*(VBR_linBackstress.out.anelastic.backstress_linear.J1+1i.*VBR_linBackstress.out.anelastic.backstress_linear.J2).^-1./VBR.in.SV.f_reshape; %complex viscosity
VBR_peak.out.anelastic.xfit_premelt.eta_complex               = -1i*(VBR_peak.out.anelastic.xfit_premelt.J1+1i.*VBR_peak.out.anelastic.xfit_premelt.J2)                            .^-1./VBR.in.SV.f_reshape; %complex viscosity
VBR_HTB.out.anelastic.xfit_premelt.eta_complex                = -1i*(VBR_HTB.out.anelastic.xfit_premelt.J1+1i.*VBR_HTB.out.anelastic.xfit_premelt.J2)                              .^-1./VBR.in.SV.f_reshape; %complex viscosity
VBR_disl.out.anelastic.maxwell_analytical.eta_complex         = -1i*(VBR_disl.out.anelastic.maxwell_analytical.J1+1i.*VBR_disl.out.anelastic.maxwell_analytical.J2)                .^-1./VBR.in.SV.f_reshape; %complex viscosity

%effective viscosity
VBR_linBackstress.out.anelastic.backstress_linear.eta_eff     = abs(VBR_linBackstress.out.anelastic.backstress_linear.eta_complex);%Pas, magnitude of complex viscosity (i.e., effective viscosity)
VBR_peak.out.anelastic.xfit_premelt.eta_eff                   = abs(VBR_peak.out.anelastic.xfit_premelt.eta_complex);%Pas, magnitude of complex viscosity (i.e., effective viscosity)
VBR_HTB.out.anelastic.xfit_premelt.eta_eff                    = abs(VBR_HTB.out.anelastic.xfit_premelt.eta_complex);%Pas, magnitude of complex viscosity (i.e., effective viscosity)
VBR_disl.out.anelastic.maxwell_analytical.eta_eff             = abs(VBR_disl.out.anelastic.maxwell_analytical.eta_complex);%Pas, magnitude of complex viscosity (i.e., effective viscosity)


%% Combine complex compliances and calculate attenuation and effective modulus for all
disp("Calculating effective complex compliances, attenuation")
% find unrelaxed modulus
Ju = (1./VBR_linBackstress.out.elastic.anharmonic.Gu); %GPa, unrelaxed shear compliance, is equal for all three models

% Set NaNs to 0 in J2 solution of maxwell_analytical (e.g., those that commonly result from the backstress model breaking down at small grain sizes and stresses)
VBR_disl.out.anelastic.maxwell_analytical.J2(isnan(VBR_disl.out.anelastic.maxwell_analytical.J2))=0;

% Sum J1 and J2s to get invQ and Geff
J1_tot =  VBR_HTB.out.anelastic.xfit_premelt.J1 + VBR_peak.out.anelastic.xfit_premelt.J1 + VBR_linBackstress.out.anelastic.backstress_linear.J1 + VBR_disl.out.anelastic.maxwell_analytical.J1 - 3*Ju; %subtracting Ju twice as it is incorporated in J1 of both the HTB, peak, and backstress models%
J2_tot =  VBR_HTB.out.anelastic.xfit_premelt.J2 + VBR_peak.out.anelastic.xfit_premelt.J2 + VBR_linBackstress.out.anelastic.backstress_linear.J2 + VBR_disl.out.anelastic.maxwell_analytical.J2;
invQ_tot = J2_tot./J1_tot; %calculate attenuation, Q-1
G_tot = 1./(J1_tot+1i.*J2_tot); %calculate combined complex shear modulus
G_eff_tot = abs(G_tot); %calculate effective shear modulus
Vs_tot      =  (G_eff_tot./VBR.in.SV.rho).^0.5;%km/s, shear-wave velocity
Eta_complex_tot     = -1i*G_tot./VBR.in.SV.f_reshape; %complex viscosity
% Eta_complex_tot(isnan(Eta_complex_tot))=0; %replace NaNs with zeros
Eta_eff_tot =  abs(Eta_complex_tot);%Pas, magnitude of complex viscosity (i.e., effective viscosity)
%% Map: stress and grain size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("finding largest mechanism")

% find indices of conditions that you want to plot
iF = 1;

% figure out which mechanism predicts the biggest
% J2 (i.e., the dominant anelastic mechanism)
% and biggest J1 (i.e., the dominant modulus
% relaxation mechanism) and create a grid of
% indices (1, 2, or 3 for the HTB, peak, and
% linearized backstress model, respectively) as a
% function of grain size and stress.
for i = 1:length(sig)
    for j = 1:length(d)
        [~,ind] = max([VBR_HTB.out.anelastic.xfit_premelt.J2(i,j,iF) VBR_peak.out.anelastic.xfit_premelt.J2(i,j,iF) VBR_linBackstress.out.anelastic.backstress_linear.J2(i,j,iF) VBR_disl.out.anelastic.maxwell_analytical.J2(i,j,iF)]); %find out premelt (ind=1) or backstress (ind=2) mechanism predicts larger J2 (i.e., is the dominant anelastic mechanism)
        % [~,ind] = max([VBR_HTB.out.anelastic.xfit_premelt.Qinv(i,j,iF) VBR_peak.out.anelastic.xfit_premelt.Qinv(i,j,iF) VBR_linBackstress.out.anelastic.backstress_linear.Qinv(i,j,iF) VBR_disl.out.anelastic.maxwell_analytical.Qinv(i,j,iF)]); %find out premelt (ind=1) or backstress (ind=2) mechanism predicts larger invQ (i.e., is the dominant anelastic mechanism)
        Qmech(i,j) = ind;
        [~,ind] = max([VBR_HTB.out.anelastic.xfit_premelt.J1(i,j,iF) VBR_peak.out.anelastic.xfit_premelt.J1(i,j,iF) VBR_linBackstress.out.anelastic.backstress_linear.J1(i,j,iF) VBR_disl.out.anelastic.maxwell_analytical.J1(i,j,iF)]); %find out premelt (ind=1) or backstress (ind=2) mechanism predicts larger J1 (i.e., is the dominant relaxation mechanism)
        % [~,ind] = min([VBR_HTB.out.anelastic.xfit_premelt.M(i,j,iF) VBR_peak.out.anelastic.xfit_premelt.M(i,j,iF) VBR_linBackstress.out.anelastic.backstress_linear.M(i,j,iF) VBR_disl.out.anelastic.maxwell_analytical.M(i,j,iF)]); %find out premelt (ind=1) or backstress (ind=2) mechanism predicts smaller M (i.e., is the dominant relaxation mechanism)
        Mmech(i,j) = ind;
    end
end

figure()
contourf(log10(sig), log10(d/1e6*1e3), Qmech')
xlabel('log10(stress)')
ylabel('log10(grain size (mm))')
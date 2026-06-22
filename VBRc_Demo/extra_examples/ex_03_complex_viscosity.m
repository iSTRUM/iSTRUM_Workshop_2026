
addpath(getenv('VBRpath'))
addpath('helper_functions')
close all; clear all; clc
vbr_init

VBR = struct();
VBR.in.elastic.methods_list={'anharmonic';};
VBR.in.anelastic.methods_list={'andrade_analytical';};

% load in the parameter set then use set the viscosity method to use to
% a fixed, constant value for the steady state viscosity.
% the value here corresponds to the maxwell viscosity for a maxwell time of
% 1000 years and an unrelaxed modulus of 60 GPa.
VBR.in.anelastic.andrade_analytical = Params_Anelastic('andrade_analytical');
VBR.in.anelastic.andrade_analytical.viscosity_method = 'fixed';
VBR.in.anelastic.andrade_analytical.eta_ss = 1.888272e+21;

% set state variables
n1 = 1;
VBR.in.SV.rho = 3300 * ones(n1,1); % density [kg m^-3]
% VBR.in.SV.P_GPa = 2 * ones(n1,1); % pressure [GPa]
% VBR.in.SV.T_K = 1473 * ones(n1,1); % temperature [K]
VBR.in.SV.f = logspace(-13,1,100);
VBR.in.elastic.Gu_TP = 60*1e9;
VBR.in.elastic.Ku_TP = 110*1e9;
VBR.in.elastic.quiet = 1;

VBR = VBR_spine(VBR) ;

% extract variables for convenience
tau_M = VBR.out.anelastic.andrade_analytical.tau_M;
eta_ss = VBR.in.anelastic.andrade_analytical.eta_ss;
Gu = VBR.out.elastic.anharmonic.Gu;
J1 = VBR.out.anelastic.andrade_analytical.J1;
J2 = VBR.out.anelastic.andrade_analytical.J2;

% complex_viscosity function: in dev, copied to ./helper_functions/complex_viscosity.m
% see https://github.com/vbr-calc/vbr/pull/245
[eta_star, eta_normalized, eta_app] = complex_viscosity(J1, J2, VBR.in.SV.f, Gu, tau_M);

% maxwell frequency
tau_f = 1./ tau_M;


figure()
subplot(2,1,1)
loglog(VBR.in.SV.f, eta_app, 'linewidth', 2, 'displayname', 'apparent viscosity')
hold on
loglog([tau_f, tau_f], [1e12,1e24],'--r', 'linewidth', 2, 'displayname', 'maxwell frequency')
loglog([VBR.in.SV.f(1), VBR.in.SV.f(end)], [eta_ss, eta_ss],'linewidth',2, '--k', 'displayname', 'steady-state viscosity')
ylim([1e12,1e24])
legend()
ylabel('||\eta*||')

subplot(2,1,2)
Qinv = VBR.out.anelastic.andrade_analytical.Qinv;
loglog(VBR.in.SV.f, Qinv, 'linewidth', 2)
hold on
loglog([tau_f, tau_f], [min(Qinv), max(Qinv)],'--k')
ylabel('Q^{-1}')


addpath(getenv('VBRpath'))
addpath('helper_functions')
close all; clear all; clc
vbr_init

VBR = struct();
VBR.in.elastic.methods_list={'anharmonic';};

VBR.in.anelastic.methods_list={'eburgers_psp';};
VBR.in.anelastic.eburgers_psp = Params_Anelastic('eburgers_psp');
VBR.in.anelastic.eburgers_psp.eBurgerFit = 'bg_peak';
method_to_plot = 'eburgers_psp';

% set state variables
VBR.in.SV.f = logspace(-13,1,50);

T_K = linspace(700, 1300, 20) + 273; % temperature [K]
dg_um = logspace(-6,-2,25) * 1e6; % grain size [um]
[VBR.in.SV.T_K, VBR.in.SV.dg_um] = meshgrid(T_K, dg_um);
sz = size(VBR.in.SV.T_K);
VBR.in.SV.rho = 3300 * ones(sz); % density [kg m^-3]
VBR.in.SV.P_GPa = 2 * ones(sz); % pressure [GPa]
VBR.in.SV.phi = 0 * ones(sz); %
VBR.in.SV.Tsolidus_K = 1200 * ones(sz) + 273; % solidus temperature [K]


VBR = VBR_spine(VBR) ;

% extract variables for convenience
tau_M = VBR.out.anelastic.(method_to_plot).tau_M;
% eta_ss = VBR.in.anelastic.andrade_analytical.eta_ss;
Gu = VBR.out.elastic.anharmonic.Gu;
J1 = VBR.out.anelastic.(method_to_plot).J1;
J2 = VBR.out.anelastic.(method_to_plot).J2;

% complex_viscosity function: in dev, copied to ./helper_functions/complex_viscosity.m
% see https://github.com/vbr-calc/vbr/pull/245
[eta_star, eta_normalized, eta_app] = complex_viscosity(J1, J2, VBR.in.SV.f, Gu, tau_M);

% maxwell frequency
Qinv = VBR.out.anelastic.eburgers_psp.Qinv;

disp(size(eta_app))

figure()
idg = find_index_where(dg_um/1e6, 1e-3);
for idk = 1:numel(T_K)
    clr = [0, 0, 1] + [1, 0, -1] * idk/numel(T_K);
    clr(clr>1) = 1;
    clr(clr<0) = 0;

    subplot(2,1,1)
    hold on
    loglog(VBR.in.SV.f, Qinv(idg,idk,:), 'linewidth', 2, 'color', clr)

    subplot(2,1,2)
    hold on
    loglog(VBR.in.SV.f, eta_app(idg,idk,:), 'linewidth', 2, 'color', clr)

end
title('dg = 1 mm, varying T')
subplot(2,1,1)
ylabel('Q^{-1}')

subplot(2,1,2)
ylabel('||\eta*||')

figure()
idk = find_index_where(T_K - 273, 1200);
for idg = 1:numel(dg_um )
    clr = [0, 1, 0] + [0, -1, 1] * idg/numel(dg_um );
    clr(clr>1) = 1;
    clr(clr<0) = 0;

    subplot(2,1,1)
    hold on
    loglog(VBR.in.SV.f, Qinv(idg,idk,:), 'linewidth', 2, 'color', clr)

    subplot(2,1,2)
    hold on
    loglog(VBR.in.SV.f, eta_app(idg,idk,:), 'linewidth', 2, 'color', clr)

end
title('T = 1200 C, varying dg_um')

subplot(2,1,1)
ylabel('Q^{-1}')

subplot(2,1,2)
ylabel('||\eta*||')



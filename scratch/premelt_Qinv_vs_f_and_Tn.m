
  addpath(getenv('vbrdir'))
  vbr_init

  VBR = struct();
  %% write method list %%
  VBR.in.elastic.methods_list={'anharmonic'};
  VBR.in.anelastic.methods_list={'xfit_premelt'};
  VBR.in.anelastic.xfit_premelt.include_direct_melt_effect = 1;


  %% Define the Thermodynamic State %%
  VBR.in.SV.f = logspace(-3,-1,30); % 1 Hz

  Tn_cases = transpose(linspace(0.5, 1.05, 10));
  VBR.in.SV.Tsolidus_K = (1200+273);
  VBR.in.SV.T_K = Tn_cases .* VBR.in.SV.Tsolidus_K;
  sz=size(VBR.in.SV.T_K); % temperature [K]
  VBR.in.SV.P_GPa = full_nd(2.5, sz); % pressure [GPa]
  % remaining state variables (ISV)
  VBR.in.SV.dg_um=full_nd(3, sz); % grain size [um]
  VBR.in.SV.rho = full_nd(3300, sz); % density [kg m^-3]
  VBR.in.SV.sig_MPa = full_nd(1, sz); % differential stress [MPa]
  VBR.in.SV.phi = full_nd(0, sz);


  VBR = VBR_spine(VBR);
  figure()
  for it = 1:numel(Tn_cases)
    hold all
    clr = [0, 0, 1] + [1, 0, -1] * (Tn_cases(it)-Tn_cases(1)) ./ (Tn_cases(end) - Tn_cases(1));
    T_K = VBR.in.SV.T_K(it,1);
    leg = [num2str(round(Tn_cases(it)*100)/100), ' T_K: ',num2str(round(T_K))];
    loglog(VBR.in.SV.f, squeeze(VBR.out.anelastic.xfit_premelt.Qinv(it,:)), ...
            'displayname', leg,...
             'color', clr, 'linewidth',2)
  end
xlabel('f [hz]')
ylabel("Q^{-1}")

% Axes properties
set(gca, ...
    'fontsize', 16, ...
    'linewidth', 2, ...
    'box', 'on', ...         % surrounding box
    'xminortick', 'on', ...
    'yminortick', 'on');

title(['1200 C solidus, 3 {\mu}m '])

% Legend outside the axes
legend('location', 'eastoutside');

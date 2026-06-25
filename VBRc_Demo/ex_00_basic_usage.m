%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basic usage of the VBRc and skills-building
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. initialize the VBRc
path_to_top_level_vbr=getenv('vbrdir');  % or /path/to/top_level/vbr/
addpath(path_to_top_level_vbr)
addpath('helper_functions')
% use the VBRc a lot? add the above to
% startup.m for matlab (https://www.mathworks.com/help/matlab/ref/startup.html)
% or ~/.octaverc for octave (https://docs.octave.org/interpreter/Startup-Files.html)
vbr_init








% 2. initialize the VBR structure
%    - set the properties and methods
%    - set thermodynamic state variable arrays  (T, phi)

VBR_list_methods

VBR.in.elastic.methods_list = {'anharmonic'};
VBR.in.anelastic.methods_list = {'eburgers_psp', 'andrade_psp', 'xfit_premelt'};

VBR.in.SV.f = logspace(-5, -1, 10);

VBR.in.SV.T_K = transpose(linspace(800, 1500, 30));
sz_T = size(VBR.in.SV.T_K);
VBR.in.SV.P_GPa = 2.5 * ones(sz_T);

VBR.in.SV.phi = zeros(sz_T);
VBR.in.SV.rho = 3300 * ones(sz_T);
VBR.in.SV.dg_um = 0.01 * 1e6 * ones(sz_T);
VBR.in.SV.Tsolidus_K=1200*ones(sz_T) + 273;

% what SVs are required? check docs!
% https://vbr-calc.github.io/vbr/vbrmethods/viscous/
% https://vbr-calc.github.io/vbr/vbrmethods/elastic/
% https://vbr-calc.github.io/vbr/vbrmethods/anelastic/
% or just try and you should get useful error messages...

% 3. call the VBRc
VBR = VBR_spine(VBR);









% 4. inspect output
% tab-complete VBR.
% shape of arrays (frequency dependence)
% units of output








% for outputs: https://vbr-calc.github.io/vbr/ left hand method types

% units of outputs









% 5. method citations








%ans =
%{
%  [1,1] = Jackson and Faul, 2010, Phys. Earth Planet. Inter., https://doi.org/10.1016/j.pepi.2010.09.005
%}










% 6. Looping over structures to compare methods
% ane_fields = fieldnames(VBR.out.anelastic)  % cell array
% ane_fields(1) % single element of cell array
% ane_fields{2} % just the string

% iT = find_index_where(VBR.in.SV.T_K, 1200+273);
% nfields = numel(ane_fields);
% figure()
% for ifield = 1:nfields
%     Qinv = VBR.out.anelastic.(ane_fields{ifield}).Qinv;

%     hold all
%     loglog(VBR.in.SV.f, Qinv(iT, :), 'displayname', ane_fields{ifield})
% end
% legend()
% ylabel('Q^{-1}')
% xlabel('f [Hz]')
% set(gca, "fontsize", 20)

% 7. saving results
% help VBR_save
% VBR_save(VBR, 'myvbr_results.mat')  % alias to save(...)
% VBR = load('myvbr_results.mat');
% VBR_save(VBR, 'myvbr_results.mat', 1) % exclude the state variable arrays


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basic usage of the VBRc and skills-building
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. initialize the VBRc
path_to_top_level_vbr=getenv('vbrdir');  % or /path/to/top_level/vbr/
addpath(path_to_top_level_vbr)
% use the VBRc a lot? add the above to
% startup.m for matlab (https://www.mathworks.com/help/matlab/ref/startup.html)
% or ~/.octaverc for octave (https://docs.octave.org/interpreter/Startup-Files.html)
vbr_init








% 2. initialize the VBR structure
%    - set the properties and methods
%    - set thermodynamic state variable arrays  (T, phi)

VBR_list_methods



% what SVs are required? check docs!
% https://vbr-calc.github.io/vbr/vbrmethods/viscous/
% https://vbr-calc.github.io/vbr/vbrmethods/elastic/
% https://vbr-calc.github.io/vbr/vbrmethods/anelastic/
% or just try and you should get useful error messages...

% 3. call the VBRc
% VBR = VBR_spine(VBR);









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
% viscous_fields = fieldnames(VBR.out.viscous)  % cell array
% viscous_fields(1) % single element of cell array
% viscous_fields{2} % just the string

% nfields = numel(viscous_fields);
% figure()
% for ifield = 1:nfields
%     visc = VBR.out.viscous.(viscous_fields{ifield});
%     eta_total = visc.eta_total;
%     hold all
%     semilogy(VBR.in.SV.T_K-273, eta_total, 'displayname', viscous_fields{ifield})
% end
% legend()
% ylabel(['\eta_{ss} (', VBR.out.viscous.HZK2011.units.eta, ')'])
% xlabel('T (deg C)')
% set(gca, "fontsize", 20)

% 7. saving results
% help VBR_save
% VBR_save(VBR, 'myvbr_results.mat')  % alias to save(...)
% VBR = load('myvbr_results.mat');
% VBR_save(VBR, 'myvbr_results.mat', 1) % exclude the state variable arrays



function cout = proc_applySlidingClassi(epo,proc,C)

dt = 1000/epo.fs;
fv_wind = proc.ival;
fv_len = length(fv_wind(1):dt:fv_wind(end));
n_points = size(epo.x,1) - fv_len + 1;

cout.t = zeros(n_points,1);
cout.x = zeros(n_points,1);
T = epo.t(1) - fv_wind(1);
for ii = 1:n_points
    
    fv = proc_selectIval(epo,fv_wind+T);
    fv.t = fv.t-T;
    fv = proc_chain(fv,proc.apply);
    
    cout.t(ii) = fv_wind(2)+T;
    cout.x(ii) = apply_separatingHyperplane(C,fv.x(:));
    
    T = T+dt;
    
end

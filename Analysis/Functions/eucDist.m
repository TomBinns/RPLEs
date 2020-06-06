function dists = eucDist(epos,template)
% eucDist - Gets the Euclidean distances between the RP template and: the 
%     single-trial RPs; the single-event RPLEs (phase 1 and phase RT)
%     
% Argument(s):
%  epos -     epochs of all subjects, of a single threshold, for all 
%             three modes, containing the single-trial RPs and single-event 
%             RPLEs
%  template - RP template of the subject     
%  
% Returns:
%  dists -    Euclidean distances between the RP template and the RPs and 
%             RPLEs for all subjects and all modes, for a single threshold
%                 
% Author(s): Thomas Binns, 2020                


global opt


dists = cell(1,3);
for aa = 1:3
    dists{aa} = cell(1,length(opt.subjs_all));
    for bb = 1:length(opt.subjs_all)
        rpLoc = find(strcmp(epos{aa}{bb}.className,'movement onset'));
        rple_p1Loc = find(strcmp(epos{aa}{bb}.className,'rple p1'));
        rple_rtLoc = find(strcmp(epos{aa}{bb}.className,'rple rt'));
        dists{aa}{bb}.rp = nan(sum(epos{aa}{bb}.y(rpLoc,:)),1);
        dists{aa}{bb}.rple_p1 = nan(sum(epos{aa}{bb}.y(rple_p1Loc,:)),1);
        dists{aa}{bb}.rple_rt = nan(sum(epos{aa}{bb}.y(rple_rtLoc,:)),1);
        xx = 1;
        yy = 1;
        zz = 1;
        for cc = 1:size(epos{aa}{bb}.y,2)
            event = epos{aa}{bb}.x(:,:,cc);
            if epos{aa}{bb}.y(rpLoc,cc) == 1
                dists{aa}{bb}.rp(xx) = norm(template{aa}{bb}(:)-event(:));
                xx = xx+1;
            elseif epos{aa}{bb}.y(rple_p1Loc,cc) == 1
                dists{aa}{bb}.rple_p1(yy) = norm(template{aa}{bb}(:)-...
                    event(:));
                yy = yy+1;
            elseif epos{aa}{bb}.y(rple_rtLoc,cc) == 1
                dists{aa}{bb}.rple_rt(zz) = norm(template{aa}{bb}(:)-...
                    event(:));
                zz = zz+1;
            else
                error('An event marker is missing')
            end
        end
    end
end


end
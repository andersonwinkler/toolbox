function bayesdata = combatapply(dat,batch,mod,batch0,grand_mean,B_hat,var_pooled,gamma_star,delta_star)
% This applies the harmonization parameters found with ComBat to new data
% 
% The original ComBat can be found at https://github.com/Jfortin1/ComBatHarmonization
% 
% For information, see and cite:
% * Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray 
%   expression data using empirical Bayes methods.
%   Biostatistics. 2007 Jan;8(1):118-27.
% * Fortin JP, Parker D, Tun? B, Watanabe T, Elliott MA, Ruparel K,
%   Roalf DR, Satterthwaite TD, Gur RC, Gur RE, Schultz RT, Verma R, 
%   Shinohara RT. Harmonization of multi-site diffusion tensor imaging
%   data. Neuroimage. 2017 Nov 1;161:149-170.
% * Fortin JP, Cullen N, Sheline YI, Taylor WD, Aselcioglu I, Cook PA, 
%   Adams P, Cooper C, Fava M, McGrath PJ, McInnis M, Phillips ML, 
%   Trivedi MH, Weissman MM, Shinohara RT. Harmonization of cortical 
%   thickness measurements across scanners and sites.
%   Neuroimage. 2018 Feb 15;167:104-120.
% 
% _____________________________________
% Anderson M. Winkler
% UTRGV
% Jun/2023
% http://brainder.org

levels0   = unique(batch0);
levels    = unique(batch);
shared    = intersect(levels,levels0);
if numel(levels) ~= numel(shared) || any(levels ~= shared)
    error('Some batches were not present in the original run of ComBat');
end
n_batch   = numel(levels0);
bidx      = cell(n_batch,1);
for b = 1:n_batch
    bidx{b} = batch == b;
end
mod        = mod - mean(mod);
stand_mean = grand_mean + mod*B_hat(n_batch+1:end,:);
s_data     = (dat-stand_mean)./sqrt(var_pooled);
bayesdata  = zeros(size(s_data));
for b = 1:n_batch
    if ~ isempty(bidx{b})
        idx = bidx{b};
        bayesdata(idx,:) = (s_data(idx,:)-gamma_star(b,:))./sqrt(delta_star(b,:));
    end
end
bayesdata  = bayesdata.*sqrt(var_pooled) + stand_mean;
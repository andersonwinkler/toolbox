function bayesdata = combatapply(dat,batch,mod,batch0,grand_mean,B_hat,var_pooled,gamma_star,delta_star,mod_mean)
% Apply the harmonization parameters found with ComBat to new data.
%
% Usage:
% bayesdata = combatapply(dat,batch,mod,batch0,grand_mean, ...
%             B_hat,var_pooled,gamma_star,delta_star,mod_mean)
%
% Inputs:
% dat        : New data, N by P (samples by features), same
%              orientation as combatfit.
% batch      : Vector of length N with batch labels for dat.
% mod        : N by Q covariates for the new samples, with the
%              same columns as in training. Use [] if training
%              used no covariates.
% batch0     : Original training batch vector passed to combatfit
%              (not a recoded or subset vector). Unique levels
%              and their sort order must match the rows of
%              gamma_star, delta_star, and the batch block of B_hat.
% grand_mean : From combatfit.
% B_hat      : From combatfit.
% var_pooled : From combatfit.
% gamma_star : From combatfit.
% delta_star : From combatfit.
% mod_mean   : From combatfit.
%
% Outputs:
% bayesdata  : Harmonized data (N by P).
%
% Notes:
% * Batch location/scale, pooled variance, grand mean, and
%   covariate coefficients are not re-estimated. Fit on a
%   training fold with combatfit, then correct a test fold
%   with this function:
%   [yc,gm,Bh,vp,gs,ds,mm] = combatfit(yTrain,batchTrain,modTrain,true);
%   yTestH = combatapply(yTest,batchTest,modTest,batchTrain,gm,Bh,vp,gs,ds,mm);
%
% * New covariates are centered with the training mod_mean,
%   then standardized with the training grand_mean, covariate
%   coefficients in B_hat, and var_pooled. Each row is location
%   and scale corrected with the training gamma_star and
%   delta_star of its batch, and mapped back to the original
%   scale. Subject-specific covariates are removed before the
%   scale step and added back after, as in original ComBat.
%
% * If training used covariates, mod must have the same number
%   of columns. Passing [] in that case is an error. For
%   predictive cross-validation, put in mod only covariates
%   that are known at test time; including the outcome leaks
%   that label into the harmonized features.
%
% * Batches that were not present in batch0 are left unadjusted
%   (returned as in dat) and a warning is issued.
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

if numel(batch) == max(size(batch))
    batch = batch(:);
else
    error('"batch" must be a vector.')
end
levels0   = unique(batch0);
n_batch   = numel(levels0);
if any(~ismember(batch, levels0(:)))
    warning(['Some batches were not present in the original run of ComBat. ' ...
             'Those rows will be left unadjusted.']);
end
bidx      = cell(n_batch,1);
for b = 1:n_batch
    bidx{b} = batch == levels0(b);
end
if isempty(mod)
    mod = zeros(size(dat,1),0);
elseif numel(mod_mean) ~= size(mod,2)
    error('mod_mean must have as many columns as mod.');
else
    mod = mod - mod_mean;
end
stand_mean = grand_mean + mod*B_hat(n_batch+1:end,:);
s_data     = (dat-stand_mean)./sqrt(var_pooled);
bayesdata  = s_data;
for b = 1:n_batch
    idx = bidx{b};
    if any(idx)
        bayesdata(idx,:) = (s_data(idx,:)-gamma_star(b,:))./sqrt(delta_star(b,:));
    end
end
bayesdata  = bayesdata.*sqrt(var_pooled) + stand_mean;
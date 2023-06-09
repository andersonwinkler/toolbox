function [bayesdata,grand_mean,B_hat,var_pooled,gamma_star,delta_star] = combat(dat,batch,mod,parametric)
% This is an edit of the original combat.m such that the
% harmonization parameters are produced.
% Other edits make the parametric adjustment much faster for large
% datasets. The results should be identical to the original combat.m
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

% Ensure data has correct sizes
N = numel(batch);
if N == max(size(batch))
    batch = batch(:);
else
    error('"batch" must be a vector.')
end
if size(dat,1) ~= N
    error('"dat" must have as many rows as number of elements in "batch"');
end
if numel(mod) > 0 && size(mod,1) ~= N
    error('"mod" must have as many rows as number of elements in "batch"');
end
nD = size(dat,2);

% Remove constant columns from the data
icte_y = ~ sum(diff(dat,1,1).^2,1);
dat(:,icte_y) = [];

% Create regressors for the batches and design matrix
levels    = unique(batch);
n_batch   = numel(levels);
bidx      = cell(n_batch,1);
for b = 1:n_batch
    bidx{b} = batch == b;
end
fprintf('[combat] Found %d batches\n', n_batch);
batchmod  = double(batch == levels');
n_batches = sum(batchmod);
mod       = mod-mean(mod);
design    = [batchmod mod];
n_cov     = size(mod,2);

fprintf('[combat] Adjusting for %d covariate(s) of covariate level(s)\n',n_cov)
if rank(mod) < size(mod,2)
    error('Error. The covariates are confounded. Please remove one or more of the covariates so the design is not confounded.');
end
if rank(design) < size(design,2)
    error('Error. At least one covariate is confounded with batch. Please remove confounded covariates and rerun ComBat.')
end

fprintf('[combat] Standardizing Data across features\n')
B_hat = pinv(design)*dat;

% Standarization model
grand_mean = n_batches*B_hat(1:n_batch,:)/N;
var_pooled = sum((dat-design*B_hat).^2,1)/N;
idx = var_pooled == 0;
var_pooled(idx) = median(var_pooled(~idx));
stand_mean = grand_mean + mod*B_hat(n_batch+1:end,:);
s_data = (dat-stand_mean)./sqrt(var_pooled);

% Get regression batch effect parameters
fprintf('[combat] Fitting L/S model and finding priors\n')
gamma_hat = pinv(batchmod)*s_data;
delta_hat = zeros(n_batch,nD);
for b = 1:n_batch
    delta_hat(b,:) = var(s_data(bidx{b},:),[],1);
end

% Find parametric priors:
gamma_bar = mean(gamma_hat,2);
t2 = var(gamma_hat,[],2);
m  = mean(delta_hat,2);
s2 = var(delta_hat,[],2);
a_prior = (2.*s2+m.^2)./s2;
b_prior = (m.*s2+m.^3)./s2;

gamma_star = zeros(n_batch,nD);
delta_star = zeros(n_batch,nD);
if parametric
    fprintf('[combat] Finding parametric adjustments\n')
    for b=1:n_batch
        [gamma_star(b,:),delta_star(b,:)] = itSol(...
            s_data(bidx{b},:),...
            gamma_hat(b,:),...
            delta_hat(b,:),...
            gamma_bar(b),...
            t2(b),...
            a_prior(b),...
            b_prior(b),...
            0.001);
    end
else
    fprintf('[combat] Finding non-parametric adjustments\n')
    for b = 1:n_batch
        [gamma_star(b,:),delta_star(b,:)] = inteprior(...
            s_data(bidx{b},:),...
            gamma_hat(b,:),...
            delta_hat(b,:));
    end
end

fprintf('[combat] Adjusting the Data\n')
bayesdata = zeros(size(s_data));
for b = 1:n_batch
    idx = bidx{b};
    bayesdata(idx,:) = (s_data(idx,:)-batchmod(idx,:)*gamma_star)./sqrt(delta_star(b,:));
end
bayesdata = bayesdata.*sqrt(var_pooled) + stand_mean;

function [g_new,d_new] = itSol(sdat,g_hat,d_hat,g_bar,t2,a,b,conv)
g_old = g_hat;
d_old = d_hat;
change = Inf;
n = size(sdat,1);
while change > conv
    g_new = (t2*n.*g_hat + d_old.*g_bar)./(t2*n + d_old);
    sum2  = sum(((sdat-sum(g_new,1)).^2),1);
    d_new = (.5.*sum2+b)./(n/2+a-1);
    change = max(max(abs(g_new-g_old)./g_old), max(abs(d_new-d_old)./d_old));
    g_old = g_new;
    d_old = d_new;
end

function [gstar,dstar] = inteprior(sdat, ghat, dhat)
[n,nD] = size(sdat);
gstar  = zeros(1,nD);
dstar  = zeros(1,nD);
idx    = true(1,nD);
for i = 1:nD
    idx(i) = false;
    sum2 = sum((sdat(:,i)-ghat(idx)).^2,1);
    LH = exp(-sum2./(2.*dhat(idx)))./(2*pi*dhat(idx)).^(n/2);
    gstar(i) = sum(ghat(idx).*LH)./sum(LH);
    dstar(i) = sum(dhat(idx).*LH)./sum(LH);
    idx(i) = true;
end

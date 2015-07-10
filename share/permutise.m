% #!/usr/bin/octave -q
function permutise(datafile,mtxfile,confile,outprefix,outextens,nP)
% Fit a model and compute p-values based on permutation.
%
% Usage:
% permutise(datafile,mtxfile,confile,outprefix,nperm)
%
% datafile:  A table separated by spaces containing one subject
%            per column and one face, vertex or voxel per line.
%            The transposed datafile is "Y".
% mtxfile:   Design matrix for the GLM. One subject per line, one
%            regressor per column. Has to be space-separated values.
% confile:   Contrast file. Use one contrast per file. Also space-
%            separated values.
% outprefix: File prefix for the outputs. Can be a full path, and
%            ideally should contain the hemisphere name (lh/rh).
% outextens: File extension for the outputs. Suggested 'dpf' or 'dpv'.
% nperm:     Number of permutations.
%
% The outputs are files named as:
% *.tstat.*    : Contain the statistics for the non-permuted (correct) model
% *.parmpval.* : Contain the parametric p-values
% *.permpval.* : Contain the non-parametric (permutation) p-values
% *.counter.*  : Contain the number of occurrences of a statistic higher (more
%                positive) than the one observed with the correct value
% *.maxt.txt   : List of the maximum statistic observed for each realization.
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Jul/2011
% http://brainder.org

% Load the contrast matrix
C = dlmread(confile,' ');

% Load the design matrix
X  = dlmread(mtxfile,' ');   % Design matrix
Xi = X(:,C ~= 0);            % Columns of interest of X
Xn = X(:,C == 0);            % Nuisance columns of X
N  = size(X,1);

% Load the actual data
fid = fopen(datafile);
Y = fscanf(fid,'%f',Inf);
fclose(fid);
nD = numel(Y)/N;
Y = reshape(Y,[N nD]);

% Print some diagnostics
fprintf('- Observed data: %s\n',datafile);
fprintf('- Design matrix: %s\n',mtxfile);
fprintf('- Contrast: %s\n',confile);
fprintf(strcat('- Contrast vector: [ ',repmat('%g ',[ 1 numel(C)]),']\n'),C);
fprintf('- Number of observations (rows in X): %d\n',N);
fprintf('- Number of regressors (cols in X): %d\n',size(X,2));
fprintf('- Number of discrete points (cols in Y): %d\n',nD);
fprintf('- Number of variables of interest: %d\n',sum(C~=0));
fprintf('- Number of nuisance variables: %d\n\n',sum(C==0));
fprintf('- Number of permutations: %d\n\n',nP);

% GLM for the correct model (no permutation)
fprintf('Running correct model.\n');
b    = X\Y;                             % Parameter estimates (betas)
r    = Y-X*b;                           % Residuals (errors)
df   = N-rank(X);                       % Degrees of freedom
sig2 = sum(r.^2)/df;                    % Variance
varC = sig2.*sum((C*pinv(X'*X))'.*C');  % Contrast variance
T    = C*b./(varC.^(1/2));              % T statistic
parmpvals = 1-tcdf(T,df);    % Compute parametric p-values (bonus)

if nP,  % If nP > 0, then compute permutations, otherwise, it's just parametric

    % Initialise some variables
    iT   = zeros(size(T));                  % To increment T
    imaxT  = zeros(nP,1);                   % To store to max T
    iminT  = zeros(nP,1);                   % To store to min T (reverse contrast)

    % Fit the nuisance only and take the residuals
    bn = Xn\Y;
    r = Y-Xn*bn;

    % Loop over permutations (with replacement)
    for p = 1:nP,
        fprintf('Running permutation %d.\n',p);

        % Generate a random indexer, making sure that the permutation doesn't
        % happen only within groups of dicrete values in the part of the design
        % matrix that contains columns with the effects of interest.
        whiletrick = true;
        while whiletrick,
            [ignore,idx] = sort(rand(N,1));
            tmp = (Xi(idx,:) == Xi);
            whiletrick = all(tmp(:));
        end

        % Permute the residuals and put the nuisance back
        r = r(idx,:);
        Yp = r+Xn*bn;

        % Fit the full model and get the statistic
        b     = X\Yp;
        r     = Yp-X*b;
        sig2f = sum(r.^2)/df;
        varCf = sig2f.*sum((C*pinv(X'*X))'.*C');
        Tf    = C*b./(varCf.^(1/2));

        % Increment the cumulative counters
        iT     = iT  + (Tf >= T);
        imaxT(p) = max(Tf(:));
        iminT(p) = min(Tf(:));
    end
    fprintf('Finished permutations.\n');
end

% Prepare to save
fprintf('Saving the results.\n');
crd  = zeros(nD,3);
idx  = 0:nD-1;

% Save the T statistic (not permuted)
fid = fopen(sprintf('%s.tstat.%s',outprefix,outextens),'w');
fprintf(fid,'%0.3d %g %g %g %0.16f\n',[idx' crd T']');
fclose(fid);

% Save the parametric p-value
fid = fopen(sprintf('%s.parmpval.%s',outprefix,outextens),'w');
fprintf(fid,'%0.3d %g %g %g %0.16f\n',[idx' crd parmpvals']');
fclose(fid);

% Only save permutation results if permutations were really carried out
if nP,

    % Convert the counter to p=values (nonparametric)
    permpvals = iT/nP;

    % Save the counter
    fid = fopen(sprintf('%s.counter.%s',outprefix,outextens),'w');
    fprintf(fid,'%0.3d %g %g %g %d\n',[idx' crd iT']');
    fclose(fid);

    % Save the P-value statistic
    fid = fopen(sprintf('%s.permpval.%s',outprefix,outextens),'w');
    fprintf(fid,'%0.3d %g %g %g %0.16f\n',[idx' crd permpvals']');
    fclose(fid);

    % Save the max T distribution and the max T value
    fid = fopen(sprintf('%s.maxt.txt',outprefix),'w');
    fprintf(fid,'%0.16f\n',imaxT);
    fclose(fid);

    % Save the max T distribution and the min T value
    fid = fopen(sprintf('%s.mint.txt',outprefix),'w');
    fprintf(fid,'%0.16f\n',iminT);
    fclose(fid);
end
% Done!
fprintf('Done!.\n');

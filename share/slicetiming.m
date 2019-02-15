function slicetiming(varargin)
% Run slice timing correction using the algorithm originally from 
% Geoffrey Aguirre and Eric Zarahn, later integrated into SPM.
% 
% Usage:
% slicetiming(filein,TR,refslice,sliceorder,fileout)
% 
% - filein     : Input 4D file to be corrected.
% - TR         : Repetition time.
% - reference  : Reference slice or reference time.
% - sliceorder : Vector of slice acquisition order, or timings (in s)
%                or a "slice_code" as in the NIFTI specification.
%                This assumes the 1st slice is the bottom. If
%                unsure, open the file in MATLAB first to confirm.
% - fileout    : Output file name
% 
% This version uses the same implementation from SPM, originally
% written by Darren Gitelman, without all the GUI elements and 
% other SPM-centric features.
% This version requires PALM (for file i/o).
% 
% _____________________________________
% Anderson M. Winkler
% Hospital Israelita Albert Einstein
% May/2017
% http://brainder.org

% ======================================================================
%      From spm_slice_timing_.m
% ======================================================================
% sliceorder  - slice acquisition order, a vector of integers, each
%               integer referring the slice number in the image file
%               (1=first), and the order of integers representing their
%               temporal acquisition order
%               OR vector containig the acquisition time for each slice
%               in milliseconds
% refslice    - slice for time 0
%               OR time in milliseconds for the reference slice
%__________________________________________________________________________
%
%   Note: The sliceorder arg that specifies slice acquisition order is
%   a vector of N numbers, where N is the number of slices per volume.
%   Each number refers to the position of a slice within the image file.
%   The order of numbers within the vector is the temporal order in which
%   those slices were acquired.
%
%   To check the order of slices within an image file, use the SPM Display
%   option and move the crosshairs to a voxel co-ordinate of z=1.  This
%   corresponds to a point in the first slice of the volume.
%
%   The function corrects differences in slice acquisition times.
%   This routine is intended to correct for the staggered order of
%   slice acquisition that is used during echoplanar scanning. The
%   correction is necessary to make the data on each slice correspond
%   to the same point in time. Without correction, the data on one
%   slice will represent a point in time as far removed as 1/2 the TR
%   from an adjacent slice (in the case of an interleaved sequence).
%
%   This routine "shifts" a signal in time to provide an output
%   vector that represents the same (continuous) signal sampled
%   starting either later or earlier. This is accomplished by a simple
%   shift of the phase of the sines that make up the signal.
%
%   Recall that a Fourier transform allows for a representation of any
%   signal as the linear combination of sinusoids of different
%   frequencies and phases. Effectively, we will add a constant
%   to the phase of every frequency, shifting the data in time.
%
%   Shifter - This is the filter by which the signal will be convolved
%   to introduce the phase shift. It is constructed explicitly in
%   the Fourier domain. In the time domain, it may be described as
%   an impulse (delta function) that has been shifted in time the
%   amount described by TimeShift.
%
%   The correction works by lagging (shifting forward) the time-series
%   data on each slice using sinc-interpolation. This results in each
%   time series having the values that would have been obtained had
%   the slice been acquired at the same time as the reference slice.
%
%   To make this clear, consider a neural event (and ensuing hemodynamic
%   response) that occurs simultaneously on two adjacent slices. Values
%   from slice "A" are acquired starting at time zero, simultaneous to
%   the neural event, while values from slice "B" are acquired one
%   second later. Without corection, the "B" values will describe a
%   hemodynamic response that will appear to have began one second
%   EARLIER on the "B" slice than on slice "A". To correct for this,
%   the "B" values need to be shifted towards the Right, i.e., towards
%   the last value.
%
% Written by Darren Gitelman at Northwestern U., 1998
%
% Based (in large part) on ACQCORRECT.PRO from G. Aguirre and E. Zarahn
% at U. Penn.
%
% Modified by R. Henson, C. Buechel and J. Ashburner, FIL, to
% handle different reference slices and memory mapping.
%
% Modified by M. Erb, at U. Tuebingen, 1999, to ask for non-continuous
% slice timing and number of sessions.
%
% Modified by R. Henson for more general slice order and SPM2.
%
% Modified by A. Hoffmann, M. Woletz and C. Windischberger from Medical
% University of Vienna, Austria, to handle multi-band EPI sequences.
%__________________________________________________________________________
% Copyright (C) 1998-2014 Wellcome Trust Centre for Neuroimaging
% 
% Darren Gitelman et al.
% $Id: spm_slice_timing.m 6130 2014-08-01 17:41:18Z guillaume $

try
    % Get the inputs
    varargin = argv();

    % Disable memory dump on SIGTERM
    sigterm_dumps_octave_core(0);

    % Print usage if no inputs are given
    if isempty(varargin) || strcmp(varargin{1},'-q')
        fprintf('Run slice timing correction using the algorithm originally from\n');
        fprintf('Geoffrey Aguirre and Eric Zarahn, later integrated into SPM.\n');
        fprintf('\n');
        fprintf('Usage:\n');
        fprintf('slicetiming(filein,TR,refslice,sliceorder,fileout)\n');
        fprintf('\n');
        fprintf('- filein     : Input 4D file to be corrected.\n');
        fprintf('- TR         : Repetition time.\n');
        fprintf('- refslice   : Reference slice.\n');
        fprintf('- sliceorder : Vector of slice acquisition order, or timings (in ms).\n');
        fprintf('               or a "slice_code" as in the NIFTI specification.\n');
        fprintf('               This assumes the 1st slice (index 1 in MATLAB) is the\n');
        fprintf('               bottom of the brain. If unsure, open the file in\n');
        fprintf('               MATLAB first to confirm.\n');
        fprintf('- fileout    : Output file name.\n');
        fprintf('\n');
        fprintf('This version uses the same implementation from SPM, originally\n');
        fprintf('written by Darren Gitelman, without all the GUI elements and\n');
        fprintf('other SPM-centric features.\n');
        fprintf('This version requires PALM (for file i/o).\n');
        fprintf('\n');
        fprintf('_____________________________________\n');
        fprintf('Anderson M. Winkler\n');
        fprintf('Hospital Israelita Albert Einstein\n');
        fprintf('May/2017\n');
        fprintf('http://brainder.org\n');
        return;
    end
end

% Take inputs:
narginchk(5,5);
fileI      = varargin{1};
TR         = varargin{2};
refslice   = varargin{3};
sliceorder = varargin{4};
fileO      = varargin{5};
if ischar(TR),         TR         = eval(TR);         end
if ischar(refslice),   refslice   = eval(refslice);   end
if ischar(sliceorder), sliceorder = eval(sliceorder); end

% Read input file, prepare output:
I = palm_miscread(fileI,false);
O = I;
O.filename = fileO;

% Sizes for later:
sz  = size(I.data);
nX  = sz(1); nY = sz(2); nZ = sz(3); nV = sz(4);
len = 2^(floor(log2(nV))+1);

% Prepare a slice order vector based on the NIFTI specification
if numel(sliceorder) == 1
    switch sliceorder
        case 1
            sliceorder = 1:nZ;
        case 2
            sliceorder = nZ:-1:1;
        case 3
            sliceorder = [1:2:nZ 2:2:nZ];
        case 4
            sliceorder = [nZ:-2:1 nZ-1:-2:1];
        case 5
            sliceorder = [2:2:nZ 1:2:nZ];
        case 6
            sliceorder = [nZ-1:-2:1 nZ:-2:1];
    end
end

% Timings:
TA        = TR - TR/nZ;
timing(1) = TA/(nZ-1);
timing(2) = TR - TA;

% Some general processing:
if nZ ~= numel(sliceorder)
    error('Number of slices in the file don''t match the vector with slice orderings.');
end
if ~ isequal(1:nZ,sort(sliceorder))
    if ~ all(sliceorder >= 0 & sliceorder <= TR)
        error('Input is neither slice indices nor slice times.');
    end
    unit = 'times';
else
    if ~ ismember(refslice,sliceorder)
        error('Reference slice should contain a slice index.');
    end
    unit = 'indices';
end

% Set up [time x voxels] matrix for holding image info
stack = zeros([len nX]);

% Compute shifting amount from reference slice and slice order
if isequal(unit,'times')
    shiftamount = (sliceorder-refslice)/TR;
else
    rslice      = find(sliceorder == refslice);
    [~,idx]     = sort(sliceorder);
    shiftamount = (idx-rslice) * timing(1)/TR;
end

% For loop to perform correction slice by slice
for z = 1:nZ
    slices = squeeze(I.data(:,:,z,:));
    phi    = zeros(1,len);
    
    % Phi represents a range of phases up to the Nyquist frequency
    % Shifted phi 1 to right.
    for f = 1:len/2
        phi(f+1) = -1*shiftamount(z)*2*pi/(len/f);
    end
    
    % Mirror phi about the center
    % 1 is added on both sides to reflect Matlab's 1 based indices
    phi(len/2+1+1:len) = -fliplr(phi(1+1:len/2));
    
    % Transform phi to the frequency domain and take the complex transpose:
    shifter = [cos(phi) + sin(phi)*sqrt(-1)].';
    shifter = shifter(:,ones(size(stack,2),1));
    
    % Loop over columns:
    for y = 1:nY
        
        % Extract columns from slices:
        stack(1:nV,:) = reshape(slices(:,y,:),[nX nV])';
        
        % Fill in continous function to avoid edge effects:
        for g = 1:size(stack,2)
            stack(nV+1:end,g) = linspace(stack(nV,g),stack(1,g),len-nV)';
        end
        
        % Shift the columns:
        stack = real(ifft(fft(stack,[],1).*shifter,[],1));
        
        % Re-insert shifted columns:
        slices(:,y,:) = reshape(stack(1:nV,:)',[nX 1 nV]);
    end
    
    % Write out the slice for all volumes:
    for v = 1:nV
        O.data(:,:,z,v) = slices(:,:,v);
    end
end

% Save to disk:
palm_miscwrite(O);

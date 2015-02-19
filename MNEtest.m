% Generate a test data set for MNE receptive field estimation
% Brad Theilman, Marvin Thielk
% 18 February 2015

clear
close all

%% Data Generation
binsize = 1;


stimuli = [];
totspikes = [];
njack = 4;

% We want 16 x 16 receptive fields
sz = 16;
Ndim = sz^2;

% Generate Gabor Patches
gb1 = gabor_fn(1, pi/4, 6, 3, 0, sz);
gb2 = gabor_fn(1, pi/2, 7, 3, 0, sz);

% Generate White Noise Stimulus
stimuli = randn(Ndim, 2^14);

% Generate J matrix from model
v1 = reshape(gb1, Ndim, 1);
v2 = reshape(gb2, Ndim, 1);

qmat = zeros(Ndim, Ndim);
qmat(:, 1) = v1;
qmat(:, 2) = v2;
lambdmat = zeros(Ndim, Ndim);
lambdmat(1, 1) = 1;
lambdmat(2, 2) = -1;
J_act = qmat*lambdmat * pinv(qmat);


% Compute Output spikes
p = 1 ./ (1+ exp(sum(stimuli'.*(stimuli'*J_act), 2) ));
spikes = rand(size(p)) < p;

stimmaster = stimuli;
spikesmaster = spikes;
Nsamples = length(stimuli);

for jack = 1:njack;
    stimuli = stimmaster;
    spikes = spikesmaster;
    if njack ~= 1
        Ntest = floor(Nsamples/njack);
        teststim = stimuli(:, 1+(jack-1)*Ntest:jack*Ntest);
        testresp = spikes(1+(jack-1)*Ntest:jack*Ntest);
        stimuli(:, 1+(jack-1)*Ntest:jack*Ntest) = [];
        spikes(1+(jack-1)*Ntest:jack*Ntest) = [];
    end

    %% MNE Estimation
    pfinal = MNEfit(stimuli', spikes, teststim', testresp, 2, 0);
    
    % Sort the eigenvalues
    avgfr = pfinal(1);
    sta = pfinal(2:(Ndim+1));
    J = reshape(pfinal((Ndim+2):end), [Ndim, Ndim]);
    Jmine = J;
    [evec, eval] = eig(J);
    
    [vals, inds] = sort(diag(eval));
    evecssort = evec(:, inds);
    sz1 = 16;
    sz2 = 16;

    figure();
    colormap(gray);
    imagesc(reshape(evecssort(:, 1), [sz1 sz2]));
    title(strcat('First Eigenvector: ', num2str(vals(1))));

    
    figure();
    colormap(gray);
    imagesc(reshape(evecssort(:, end), [sz1 sz2]));
    title(strcat('Last Eigenvector: ', num2str(vals(end))));
    
    
    
    figure();
    colormap(gray);
    imagesc(reshape(sta, [sz1 sz2]));
    title('STA')
    
    figure();
    
    plot(vals, '.');
    title('Eigenvalues')
    
end



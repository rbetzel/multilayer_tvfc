% add full directory to path
addpath(genpath(pwd));

% read parameters from config file
config = loadjson('config.json.multi');

% define variables
winsz = config.winsz;
gamma = config.gamma;
omega = config.omega;
ts = h5read(config.ts,'/timeseries');

% number of time points/nodes
[P,N] = size(ts);

% calculate number of windows
T = floor(P/winsz);

% generate time-varying fc matrices
A = zeros(N,N,T);
for i = 1:T
    idx = (1:winsz) + (i - 1)*winsz;
    A(:,:,i) = corr(ts(idx,:));
end

% generate and cluster multilayer network
B = spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
for s=1:T
    idx = (1:N) + (s - 1)*N;
    b = (A(:,:,s) - gamma).*~eye(N);
    B(idx,idx) = A(:,:,s);
end

% symmetrize matrix (due to small imprecisions)
B = (B + B')/2;

% add interlayer coupling
B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);

% cluster
[S,Q] = genlouvain(B,[],false);

% reshape (output 1) -- Q3 = how do I make sure that output 1 and output 2
% are saved somewhere?
S = reshape(S,N,T);

% calculate flexibility (output 2)
flx = nanmean(S(:,1:end - 1) ~= S(:,2:end),2);

% write output as csv files
dlmwrite(config.outnamecomm,'S');
dlmwrite(config.outnameflex,'flx');
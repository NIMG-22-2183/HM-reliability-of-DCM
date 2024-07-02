% This code test free energies associated with RMEL with and without compound 
% symmetry structure for the random effect as a proxy of reliability. 
% Reliability of dynamic causal modelling of resting state magnetoencephalography 

addpath('C:\spm12')
addpath(genpath('C:'))
spm('defaults','eeg')
%=====================Baseline vs two-weeks data================================

load('Data.mat')
a              = [1:14];
T              = D(a,1);
R              = D(a,2);
for i = 1:length(a)   
    F(i,1)     =  T{i, 1}.F;
    F(i,2)     =  R{i, 1}.F;
end
N                       = 14;
fr                      = abs(F); 
YY1                     = (fr(:,1)-mean(fr(:,1)))*(fr(:,2)-mean(fr(:,2)))'; % sample covariance COV(X,Y)= E(X-mx)(Y-my))
[~,~,~,Free(:,1)]       = spm_ar_reml_HB(YY1,[],1,14);
[~,~,~,Free(:,2)]       = spm_ar_reml_HB(YY1,[],2,14); % m=2 uses compund symetry

FF = Free; 
FS_labels=10; FS_ticks=10; fs_ticks=10;
figure('color','white','units','centimeters','position',[4 4 10 10],'papersize',[10 10],'filename','E.pdf')
set(gca,'fontsize',fs_ticks)
bar(FF, 'k');
set(gca,'fontsize',FS_labels)
xticks([1:2]);
names = {'standard Cov','compound symmetry Cov'};
set(gca,'XTickLabel',names, 'fontsize',10); % 'FontWeight','bold'
% xtickangle(45+45)
xlabel('Model','fontsize',15)
ylabel('Free energy','fontsize',15)
box off
axis tight;

function [C,h,Ph,F] = spm_ar_reml_HB(YY,X,m,N)
% ReML estimation of covariance components from y*y'
% FORMAT [C,h,Ph,F] = spm_ar_reml(YY,X,m,N);
%
% YY  - (m x m) sample covariance matrix Y*Y'  {Y = (m x N) data matrix}
% X   - (m x p) design matrix
% m   - (1) order of AR(m) model
% N   - number of samples
%
%
% C   - (m x m) estimated errors = h(1)*Q{1} + h(2)*Q{2} + ...
% h   - (q x 1) ReML hyperparameters h: normalised AR coeficients
% Ph  - (q x q) conditional precision of h (unnormalised)
%
% F   - [-ve] free energy F = log evidence = p(Y|X,Q) = ReML objective
%
% Performs a Fisher-Scoring ascent on F to find ReML variance parameter
% estimates.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner & Karl Friston
% $Id: spm_ar_reml.m 5219 2013-01-29 17:07:07Z spm $

% assume a single sample if not specified 
% modifed for compound symetry test
%--------------------------------------------------------------------------

try
    N;
catch
    N  = 1;
end

% assume AR(1) if not specified
%--------------------------------------------------------------------------
try
    m;
catch
    m  = 1;
end

% ortho-normalise X
%--------------------------------------------------------------------------
if isempty(X)
    X = sparse(length(YY),1);
else
    X = orth(full(X));
end

% initialise h
%--------------------------------------------------------------------------
% m     = m + 1;
n     = length(YY);
dh    = zeros(m,1);
dFdh  = zeros(m,1);
dFdhh = zeros(m,m);
L     = zeros(m,m);


% initialise and specify hyperpriors
%--------------------------------------------------------------------------
hE  = sparse(1,1,1,m,1);
hP  = speye(m,m)/exp(32);
h   = hE;

% initialise precision components
%--------------------------------------------------------------------------
for i = 1:m 
    Q{i} = spdiags(ones(n,2),[(1 - i) (i - 1)],n,n); % AR Cov assume diag(Q)=2
end
if m == 2,
 Q{1,2} = Q{1,2}-Q{1,2};        % Set AR shape to zero and then we reconstruct it as below! 
 Q{1,2} = 1*(ones(n,n)-eye(n)); % we want all off dig to be equal to 1 askin to ar COV est
 Q{1,1} = 2*eye(n);
end 

% scale data
%--------------------------------------------------------------------------
Ys    = norm(YY,1)/N;
YY    = YY/Ys;

% ReML (EM/VB)
%--------------------------------------------------------------------------
for k = 1:640

    % compute current estimate of covariance
    %----------------------------------------------------------------------
    iC    = sparse(n,n);
%     for i = 1:m
%             iC = iC + Q{i}*h(i);
%     end
     
     iC = iC + Q{1}*(h(1));
     if m==2
     iC = iC + Q{2}*((h(2)));
     end
    
    C     = spm_inv(iC);

    % E-step: conditional covariance cov(B|y) {Cq}
    %======================================================================
    Cq    = spm_pinv(X'*iC*X);

    % M-step: ReML estimate of hyperparameters
    %======================================================================

    % Gradient dF/dh (first derivatives)
    %----------------------------------------------------------------------
    P     = C - X*Cq*X';
    U     = (iC*YY/N*iC - iC)*P;
    for i = 1:m

        % dF/dh = -trace(dF/diC*iC*Q{i}*iC)
        %------------------------------------------------------------------
        QP{i}     = Q{i}*P;
        dFdh(i)   = -trace(QP{i}*U)*N/2;

    end

    % Expected curvature E{dF/dhh} (second derivatives)
    %----------------------------------------------------------------------
    for i = 1:m
        for j = i:m

            % dF/dhh = -trace{P*Q{i}*P*Q{j}}
            %--------------------------------------------------------------
            dFdhh(i,j) = -trace(QP{i}*QP{j})*N/2;
            dFdhh(j,i) =  dFdhh(i,j);

        end
    end
    
    % add hyperpriors
    %----------------------------------------------------------------------
    e     = h     - hE;
    dFdh  = dFdh  - hP*e;
    dFdhh = dFdhh - hP;
    
    % update regulariser
    %----------------------------------------------------------------------
    if ~rem(k,8)
       L  = speye(m,m)*norm(dFdhh)/128/4;
    end
    
    % Fisher scoring: update dh = -inv(ddF/dhh)*dF/dh
    %----------------------------------------------------------------------
    Ph    = -dFdhh;
    dh    = -pinv(dFdhh - L)*dFdh;

    % preclude numerical overflow
    %----------------------------------------------------------------------
    h     = h + dh;
    
    % Convergence (1% change in log-evidence)
    %======================================================================
    dF    = dFdh'*dh;
    fprintf('%-30s: %i %30s%e\n','  ReML Iteration',k,'...',full(dF));
    if dF < 1e-1, break, end

end

% log evidence = ln p(y|X,Q) = ReML objective = F = trace(R'*iC*R*YY)/2 ...
%--------------------------------------------------------------------------
if nargout > 3
    R = P*iC;
    F = - trace(R*YY*R')/2 ...
        - e'*hP*e/2 ...
        - N*n*log(2*pi)/2 ...
        - N*spm_logdet(C)/2 ...
        + N*spm_logdet(Cq)/2 ...
        -   spm_logdet(Ph)/2 ...
        +   spm_logdet(hP)/2;
end

% rescale (NB - Q{1) = 2*I
%--------------------------------------------------------------------------
C    =  C*Ys;
h    = -h(2:m)/(h(1)*2);

end


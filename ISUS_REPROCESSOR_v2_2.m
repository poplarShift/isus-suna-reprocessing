function [out,spec] = ISUS_REPROCESSOR_v2_2(CHANNEL, CHANNEL_DF, T, S, calfile, progressbarflag)
% Reprocesses absorption spectra sampled by the Satlantic ISUS/SUNA using co-located CTD data.
%
%
% v2.2: added spec structure in output...
%
% achim@npolar.no
% Aug 2015

if nargin<6
    progressbarflag=0;
end

if length(T) == 1;
    T = T*ones(size(CHANNEL,1),1) ;
end

if length(S) == 1;
    S = S*ones(size(CHANNEL,1),1) ;
end


cal=load(calfile);

% wavelength
LL=cal.lambda;
% position vector 217-240 nm
idL = find(LL>217 & LL<240);
L=LL(idL);
% calibration temperature
Tcal=20;
% ref spectrum from DIW
IR=cal.ref_DIW;
% NO3 spectrum
ENO3 = cal.ENO3 ;
% A,B,C,D coefficients for temperature dependency of ESW on sample
% temperature (Sakamoto et.al. 2009)
coefA =  1.1500276 ;
coefB =  0.02840;
coefC = -0.3101349;
coefD =  0.001222 ;


NTR=nan(size(CHANNEL,1),1);
Apmean=NTR;
rmse=NTR;

prog=0;
if progressbarflag
    init_progressbar()
end
for k=1:length(NTR)
    if progressbarflag
        prog=progressbar(k,length(NTR),prog);
    end
    % --- CALCULATE QUANTITIES (SPECTRA):
    %   I:  intensity
    %   ID: dark intensity
    %   ASM: measured absorbance
    %   ESW: estimated seawater extinction (full waveband)
    %   ASE: estimated seawater absorption (full waveband)
    %   Ap:  absorption spectrum corrected for seawater absorbance (incl. temp, full waveband)
    
    % signal intensity
    I=CHANNEL(k,:)';
    % dark current
    ID=CHANNEL_DF( k ,:)';
    % measured absorbance
    ASM = -log10((I-ID)./(IR-ID)) ;
    
    % in situ temperature
    Tis=T(k);
    %     if isnan(Tis)
    %         Tis=2;
    %     end
    % extinction coeff. of seawater at meas.temperature, extended beyond
    % Sakamoto et al's range... normalized to 1 PSU !
    % -- calfile ESW, scaled to in-situ temp
    % ESW  = cal.ESW(:) .* (coefA+coefB*Tis) ./ (coefA+coefB*Tcal) .* exp (coefD*(Tis-Tcal) .* (LL-210)) ;
    % -- sakamoto's fit
    ESW = (coefA+coefB*Tis) .* exp( (coefC+coefD*Tis).*(LL-210)) / 35;
    % expected absorption spectrum due to seawater
    ASE = ESW * S(k) ;
    % absorption spectrum corrected for seawater absorbance (incl. temp)
    Ap = ASM - ASE ;
    
    % --- COEFFICIENT MATRIX
    XX = [ENO3 LL] ;
    
    % coeff matrix incl. constant offset
    LXX = [ones(size(LL)) XX] ;
    % coefficient matrix in fitting waveband (idL)
    X = XX(idL,:) ;
    NN = size(XX,2)+1 ;
    index_bg = setdiff(1:NN,2)' ;
    if sum(~isnan(Ap))>length(Ap)/2 % only fit if sufficiently many data points exist
        % linear regression
        p = robustfit(X,Ap(idL),'ols') ; % ordinary least squares...
        NTR(k,1) = p(2);
        Afit = LXX * p; % fitted absorption
        ANO3 = ENO3 * p(2) ; % NO3 absorption
        rmse(k,1) = sqrt(nanmean((Afit(idL)-Ap(idL)).^2)) ; % RMS error of fit
        pout(k,:)=p;
        Abg = LXX(:,index_bg) * p(index_bg); % background absorption
    else
        Afit = nan(size(Ap));
        Abg = Afit ;
        ANO3 = Afit ;
        NTR(k,1)=NaN;
        Apmean(k,1)=NaN;
        p=nan(1,NN);
    end
    
    % record some spectral info for later plotting
    spec.full.Ap  (k,:) = Ap  ;
    spec.full.ASE (k,:) = ASE ;
    spec.full.ASM (k,:) = ASM ;
    spec.full.ANO3(k,:) = ANO3;
    spec.full.Afit(k,:) = Afit;
    spec.full.Abg (k,:) = Abg ;
    spec.full.I   (k,:) = I   ;
    spec.full.ID  (k,:) = ID  ;
    
    spec.fit.Ap  (k,:) = Ap  (idL) ;
    spec.fit.ASE (k,:) = ASE (idL) ;
    spec.fit.ASM (k,:) = ASM (idL) ;
    spec.fit.ANO3(k,:) = ANO3(idL) ;
    spec.fit.Afit(k,:) = Afit(idL) ;
    spec.fit.Abg (k,:) = Abg (idL) ;
    spec.fit.I  (k,:) = I   (idL) ;
    spec.fit.ID (k,:) = ID  (idL) ;
    
    spec.p (k,:) = p ;
    
end

spec.full.IR (k,:) = IR      ;
spec.fit.IR  (k,:) = IR(idL) ;

spec.full.lambda = LL';
spec.fit.lambda = L';

% create header
HEADER=create_readme(mfilename,'') ;

out.HEADER=HEADER;
out.NTR=NTR;
out.T=T;
out.S=S;
out.rmse=rmse;

if progressbarflag
    fprintf('\n')
end

%%--------------------------------------------------------%%
%%                   TVP-VAR forecasting package                      %%
%%  This code is modified from Primiceri (2005), Nakajima (2011) and  Chan and Eisenstat (2017)
%%  Function files for MCMC estimation of Time-Varying Parameter VAR model
%%  with stochastic volatility
%%  setvar.m    sets variables or options
%%  ssmooth.m   implements simulation smoother
%%  svsamp.m    implements multi-move sampler for SV model
%%  ftsvar.m    computes time-series variance  
%%  fAt.m, fXt.m, fXh.m, finvm.m  functions for some matrices 
%%  prior_ols.m  prior function for subsample




%%--------------------------------------------------------%%

close all;clear all;

nsim=10000;
nburn =0.1*nsim;       
tau = 40;
y=xlsread('awm.xlsx');% load data
yearlab=1970.25:0.25:2017.75;
my = y(tau+1:end,:);  
yearlab=yearlab(tau+1:end);

asvar = {'p'; 'x'; 'i'};    % variable names
nlag = 2;                   % # of lags

setvar('data', my, asvar, nlag); % set data

setvar('fastimp', 1);
global m_my  m_nl m_ns m_nk m_fli m_flSb m_nimp m_flfi ...
       m_iseed m_dvb0 m_dVb0 m_dva0 m_dVa0 m_dvh0 m_dVh0 m_k
   
tic;

%%--- set default options ---%%
if isempty(m_fli) == 1
    m_fli = 1;
end
if isempty(m_flSb) == 1
    m_flSb = 0;
end
if isempty(m_nimp) == 1
    m_nimp = 12 + 1;
end
if isempty(m_flfi) == 1
    m_flfi = 1;
end
if isempty(m_iseed) == 1
    m_iseed = 1;
end

rand('state', m_iseed);
randn('state', m_iseed);


%%--- set variables ---%%

ns = m_ns;  % # of time periods
T0=find(yearlab==2005.5);
nk = m_nk;  % # of series
nl = m_nl;  % # of lags
nb = nk * (nk*nl + m_fli);  % # of coefficients in beta 
na = nk * (nk-1) / 2;       % # of parameters in a
Y1 = my( nl+1:end, : );
ns=size(Y1,1);

mSigb = eye(nb) * 0.01;
mSiga = eye(na) * 0.01;
mSigh = eye(nk) * 0.01;

vidb = 1 : nb;
if m_fli == 1
    vidi = (0 : nk-1) * (nk*nl+1) + 1;
	vidb(vidi) = [];
end
[v1, v2] = find(triu(reshape(1:nk^2, nk, nk)', 1));
vida = (v1-1)*nk + v2;
%%--- Prior ---%%
[Beta,VBOLS,sigma_u] = prior_ols(y(1:tau,:),nl);
%Find lower triangular matrix At
C0=chol(sigma_u);
C0=C0./repmat(diag(C0),1,nk);
C0=inv(C0)';
C0=reshape(C0,nk*nk,1);
sigma_OLS=diag(sigma_u);
if isempty(m_dvb0) == 1
  if m_flSb == 1
    m_dvb0 = 25;          % Sigma ~ IW(vb0, I*Vb0)
    m_dVb0 = 1e-4;
  else
    m_dvb0 = 40;          % sigb_i^2 ~ IG(va0/2, Va0/2) 
    m_dVb0 = 2*1e-4;
  end
elseif m_flSb == 0
    m_dvb0 = m_dvb0*2;
    m_dVb0 = m_dVb0*2;
end   
if isempty(m_dva0) == 1
  m_dva0 = 8;             % siga_i^2 ~ IG(va0/2, Va0/2)
  m_dVa0 = 2*1e-4;    
end
if isempty(m_dvh0) == 1
  m_dvh0 = 8;             % sigh_i^2 ~ IG(vh0/2, Vh0/2)
  m_dVh0 = 2*1e-4;    
end

vb0 = reshape(Beta,nb,1);       % b_1 ~ N(b0, Sb0)
mSb0 = VBOLS;
va0 = zeros(na, 1);
for i=1:na
    va0(i,1)=C0(vida(i));
end                      % a_1 ~ N(a0, Sa0)
mSa0 = eye(na) * 10;
vh0 = log(sigma_OLS);      % h_1 ~ N(h0, Sh0)
mSh0 = eye(nk) * 10;

mS0 = eye(nb) * m_dVb0;

yhat1  = zeros(ns-1-T0+1,nk+1);  
yhat2  = zeros(ns-2-T0+1,nk+1);
yhat3  = zeros(ns-3-T0+1,nk+1);
yhat4  = zeros(ns-4-T0+1,nk+1);
yfo1  = zeros(ns-1-T0+1,nk);  
yfo2  = zeros(ns-2-T0+1,nk);
yfo3  = zeros(ns-3-T0+1,nk);
yfo4  = zeros(ns-4-T0+1,nk);

L = eye(nk);
nK = floor(m_ns/30)-1;

%Starting the recursive forecasting exercise
for t = T0:ns-1
    disp([ num2str(t-ns) ' more loops to go... ' ] );
    Y =Y1(1:t,:);
    T = size(Y,1);
    if m_fli == 1
    vy = zeros(1, nk);
    else
    vy = mean(Y);
    end
    Y = Y - ones(T, 1) * vy;
    myh = zeros(T, nk);
    mya = zeros(T, nk);
    amX = zeros(nk, nb, T);
    amXh = zeros(nk, na, T);
    amG2 = zeros(nk, nk, T);
    
    for i = nl+1 : T
    amX(:, :, i) = fXt(Y(i-nl:i-1, :), m_fli);
    end

    mb = zeros(T, nb);
    ma = zeros(T, na);
    mh = zeros(T, nk);
    dnub = m_dvb0 + T - nl - 1;
    dnua = m_dva0 + T - nl - 1;
    dnuh = m_dvh0 + T - nl - 1;
    tempyhat1 = zeros(nsim,nk+1); 
    tempyhat2 = zeros(nsim,nk+1); 
    tempyhat3 = zeros(nsim,nk+1); 
    tempyhat4= zeros(nsim,nk+1); 
    yfore1=zeros(nsim,nk);
    yfore2=zeros(nsim,nk);
    yfore3=zeros(nsim,nk);
    yfore4=zeros(nsim,nk);
        
    %%------------- S A M P L I N G   S T A R T --------------%%
    for m_k = -nburn : nsim

  %%--- sampling beta ---%%

    for i = nl+1 : T
        mAinv = finvm(fAt(ma(i, :), nk)); 
        amG2(:, :, i) = mAinv * diag(exp(mh(i,:))) * mAinv';
        mai(i, :) = mAinv(vida)';
    end
  
    mb(nl+1:end, :) ...
     = ssmooth(Y(nl+1:end,:), amX(:,:,nl+1:end), ...
               amG2(:,:,nl+1:end), mSigb, vb0, mSb0)';
    
    
  %%--- sampling a ---%%
    
    for i = nl+1 : T
       myh(i, :) = Y(i, :) - mb(i, :) * amX(:, :, i)';
       amXh(:, :, i) = fXh(myh(i, :), nk, na);
       amG2(:, :, i) = diag(exp(mh(i, :)));
    end
  
    ma(nl+1:end, :) ...
     = ssmooth(myh(nl+1:end,:), amXh(:,:,nl+1:end), ...
               amG2(:,:,nl+1:end), mSiga, va0, mSa0)';
  
  %%--- sampling h ---%%

    for i = nl+1 : T
        mya(i, :) = myh(i, :) * fAt(ma(i, :), nk)';
    end
           
    for i = 1 : nk
        mh(nl+1:end, i) ...
         = svsamp(mya(nl+1:end,i), mh(nl+1:end,i), ...
                  mSigh(i,i), vh0(i), mSh0(i,i), nk);
    end


  %%--- sampling Sigma ---%%
  
    mdif = diff(mb(nl+1:end, :));
    if m_flSb == 1 % non- diagonal
      mSb = inv(mS0 + mdif'*mdif);
      mSb = (mSb + mSb')/2;
      [mL, p] = chol(mSb, 'lower');
      if p > 0
        mSb = diag(diag(mSb));
      end
      mSigb = inv(wishrnd(mSb, dnub));
      mSigb = (mSigb + mSigb')/2;
    else
      vSb = m_dVb0 + sum(mdif.^2);
      mSigb = diag(1 ./ gamrnd(dnub/2, 2./vSb));
    end
    
    vSa = m_dVa0 + sum(diff(ma(nl+1:end, :)).^2);
    mSiga = diag(1 ./ gamrnd(dnua/2, 2./vSa));
    
    vSh = m_dVh0 + sum(diff(mh(nl+1:end, :)).^2);
    mSigh = diag(1 ./ gamrnd(dnuh/2, 2./vSh));
    if m_k > 0%% compute various forecasts
        i = m_k;
        beta_new = zeros(3,nb);
        beta_t = mb(end,:)';
        for j=1:4
        beta_new(j,:) = (beta_t + sqrt(mSigb)*randn(nb,1))';
        beta_t = beta_new(j,:)';
        end
        tmpbeta1 = reshape(beta_new(1,:),nl*nk+1,nk)';
        mu1 = tmpbeta1(:,1);
        A1 = tmpbeta1(:,2:end);
        Xt = reshape(Y(end:-1:end-nl+1,:)', nk*nl,1);    
        z1 = mu1 + A1*Xt;
        h1 = mh(end,:)' + sqrt(mSigh)*randn(nk,1);
        a1 = ma(end,:)' + sqrt(mSiga)*randn(na,1);
        L(vida) = a1; 
        Phi0 = L(:,:);
        Sig1 = Phi0*diag(exp(h1))*Phi0';
        Ytp1 = z1 +chol(Sig1,'lower')*randn(nk,1);
        yfore1(i,:)=Ytp1;
        
        
        tmpbeta2 = reshape(beta_new(2,:),nl*nk+1,nk)';
        mu2 = tmpbeta2(:,1);
        A2 = tmpbeta2(:,2:end);
        Yt1 = [Ytp1'; Y(end:-1:end-nl+2,:)]';
        Xt = reshape(Yt1(:,1:nl),nk*nl,1);
        z2 = mu2 + A2*Xt ;
        h2 = h1 + sqrt(mSigh)*randn(nk,1);
        a2 = a1 + sqrt(mSiga)*randn(na,1);
        L(vida) = a2; 
        Phi0 = L(:,:);
        Sig2 = Phi0*diag(exp(h2))*Phi0';
        Ytp2 = z2 +chol(Sig2,'lower')*randn(nk,1);
        yfore2(i,:)=Ytp2;
        
        tmpbeta3 = reshape(beta_new(3,:),nl*nk+1,nk)';
        mu3 = tmpbeta3(:,1);
        A3 = tmpbeta3(:,2:end);
        Yt2 = [Ytp2';Ytp1' ; Y(end:-1:end-nl+3,:)]';
        Xt = reshape(Yt2(:,1:nl),nk*nl,1);
        z3 = mu3 + A3*Xt ;
        h3 = h2 + sqrt(mSigh)*randn(nk,1);
        a3 = a2 + sqrt(mSiga)*randn(na,1);
        L(vida) = a3; 
        Phi0 = L(:,:);
        Sig3 = Phi0*diag(exp(h3))*Phi0';
        Ytp3 = z3 +chol(Sig3,'lower')*randn(nk,1);
        yfore3(i,:)=Ytp3;
        
        
        tmpbeta4 = reshape(beta_new(4,:),nl*nk+1,nk)';
        mu4 = tmpbeta4(:,1);
        A4 = tmpbeta4(:,2:end);
        Yt3 = [Ytp3';Ytp2';Ytp1'; Y(end:-1:end-nl+4,:)]';
        Xt = reshape(Yt3(:,1:nl),nk*nl,1);
        z4 = mu4 + A4*Xt ;
        h4 = h3 + sqrt(mSigh)*randn(nk,1);
        a4 = a3 + sqrt(mSiga)*randn(na,1);
        L(vida) = a4; 
        Phi0 = L(:,:);
        Sig4 = Phi0*diag(exp(h4))*Phi0';
        Ytp4 = z4 +chol(Sig4,'lower')*randn(nk,1);
        yfore4(i,:)=Ytp4;
        
        tempyhat1(i,:) = [log(normpdf(Y1(t+1,:),Ytp1',sqrt(diag(Sig1))')) log(mvnpdf(Y1(t+1,:)',Ytp1,Sig1)) ];
        if t<=ns-2
            tempyhat2(i,:)= [log(normpdf(Y1(t+2,:),Ytp2',sqrt(diag(Sig2))')) log(mvnpdf(Y1(t+2,:)',Ytp2,Sig2))];
        end
        if t<=ns-3
            tempyhat3(i,:) = [log(normpdf(Y1(t+3,:),Ytp3',sqrt(diag(Sig3))')) log(mvnpdf(Y1(t+3,:)',Ytp3,Sig3)) ];
        end
        if t<=ns-4
            tempyhat4(i,:) =[log(normpdf(Y1(t+4,:),Ytp4',sqrt(diag(Sig4))')) log(mvnpdf(Y1(t+4,:)',Ytp4,Sig4)) ];
        end
               
        
       
    end 
    maxy1 = max(tempyhat1);
    maxy2 = max(tempyhat2);
    maxy3 = max(tempyhat3);  
    maxy4 = max(tempyhat4); 
    meany1= mean(yfore1);
    meany2= mean(yfore2);
    meany3= mean(yfore3);
    meany4= mean(yfore4);
    end
    
        yhat1(t-T0+1,:)  = log(mean(exp(tempyhat1-repmat(maxy1,nsim,1))))+maxy1; 
        yfo1(t-T0+1,:)=meany1;
if  t<=ns-2
        yhat2(t-T0+1,:) = log(mean(exp(tempyhat2-repmat(maxy2,nsim,1))))+maxy2; 
        yfo2(t-T0+1,:)=meany2;
end
if  t<=ns-3
        yhat3(t-T0+1,:) = log(mean(exp(tempyhat3-repmat(maxy3,nsim,1))))+maxy3;
        yfo3(t-T0+1,:)=meany3;
end
if  t<=ns-4
        yhat4(t-T0+1,:) = log(mean(exp(tempyhat4-repmat(maxy4,nsim,1))))+maxy4; 
        yfo4(t-T0+1,:)=meany4;
end

end
results = zeros(4,4);

results(1,:) = sum(yhat1(4:end,:));
results(2,:) = sum(yhat2(3:end,:));
results(3,:) = sum(yhat3(2:end,:));
results(4,:) = sum(yhat4(1:end,:));
RMSFE=zeros(4,3);
RMSFE(1,:)=sqrt(mean((Y1(T0+1:end,:) - yfo1(:,:)).^2))';
RMSFE(2,:)=sqrt(mean((Y1(T0+2:end,:) - yfo2(:,:)).^2))';
RMSFE(3,:)=sqrt(mean((Y1(T0+3:end,:) - yfo3(:,:)).^2))';
RMSFE(4,:)=sqrt(mean((Y1(T0+4:end,:) - yfo4(:,:)).^2))';
clc;
fprintf('Sum of log predictive likelihoods for TVPVAR model lag %d', nlag)
fprintf('\n'); 
fprintf('                 | inflation, GDP growth, interest rate, joint\n'); 
fprintf('one-step-ahead   | %.1f,       %.1f,       %.1f,         %.1f\n', results(1,1), results(1,2), results(1,3),results(1,4)); 
fprintf('two-step-ahead   | %.1f,       %.1f,       %.1f,         %.1f\n', results(2,1), results(2,2), results(2,3),results(2,4)); 
fprintf('three-step-ahead | %.1f,       %.1f,       %.1f,         %.1f\n', results(3,1), results(3,2), results(3,3),results(3,4)); 
fprintf('four-step-ahead  | %.1f,       %.1f,       %.1f,         %.1f\n', results(4,1), results(4,2), results(4,3),results(4,4)); 
fprintf('\n'); 

fprintf('RMSFE for TVPVAR model ')
fprintf('\n'); 
fprintf('                  | inflation, GDP growth, interest rate \n'); 
fprintf('1-quarter-ahead   | %.2f,       %.2f,       %.2f\n', RMSFE(1,1), RMSFE(1,2), RMSFE(1,3)); 
fprintf('2-quarter-ahead   | %.2f,       %.2f,       %.2f\n', RMSFE(2,1), RMSFE(2,2), RMSFE(2,3)); 
fprintf('3-quarter-ahead   | %.2f,       %.2f,       %.2f\n', RMSFE(3,1), RMSFE(3,2), RMSFE(3,3)); 
fprintf('4-quarter-ahead   | %.2f,       %.2f,       %.2f\n', RMSFE(4,1), RMSFE(4,2), RMSFE(4,3)); 




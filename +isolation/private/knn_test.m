function [fnErr,fpErr,isolation,details] = knn_test(analogAp,noiseEvents,options,active)
% Exact blockwise Joshua metrics. active = [IS FN FP], independently enabled.
% FN=Nfn/(Nfn+Ns), FP=Nfp/Ns (paper equations 12/13); K is never reduced.
% Ties follow [AP;noise] input order: AP precedes noise at equal distance.
% FN needs the distance to AP neighbor h=(K+1)/2 and whether at least h
% noise neighbors are STRICTLY closer. Early rejection is therefore exact.
if nargin<3, options=struct; end
if nargin<4, active=[true true true]; end
cfg=isolation_config(options); active=logical(active);
ns=size(analogAp,1); nn=size(noiseEvents,1);
fnErr=NaN; fpErr=NaN; isolation=NaN;
details=struct('status','insufficient_events','spikes',ns,'noise',nn, ...
    'characteristicDistance',NaN,'K',cfg.K,'Completed',false(1,3), ...
    'TargetAffinity',[],'FNCount',NaN,'FPCount',NaN,'NoiseDistancesEvaluated',0, ...
    'FNDefinition','Nfn/(Nfn+Ns)','MetricVersion','Joshua2007-full-v1');
if ns<2 || nn<1, return; end
assert(size(analogAp,2)==size(noiseEvents,2),'isolation:FeatureWidth','Feature widths differ.');
assert(isreal(analogAp) && isreal(noiseEvents) && all(isfinite(analogAp),'all') && all(isfinite(noiseEvents),'all'), ...
    'isolation:InvalidFeatures','Waveforms must be real and finite.');
a=double(analogAp); b=double(noiseEvents); na=sum(a.^2,2); nb=sum(b.^2,2);
block=cfg.DistanceBlockSize; colBlock=cfg.DistanceColumnBlockSize;
canKNN=(ns+nn-1)>=cfg.K; active(2:3)=active(2:3) & canKNN;
h=(cfg.K+1)/2; d0=NaN;
if active(1)
    distanceSum=0;
    for first=1:block:ns
        rows=first:min(first+block-1,ns);
        for start=1:colBlock:ns
            cols=start:min(start+colBlock-1,ns);
            d=distances(a(rows,:),a(cols,:),na(rows),na(cols));
            d(rows(:)==cols)=0;
            distanceSum=distanceSum+sum(sqrt(d),'all');
        end
    end
    % Preserve the established Ns^2 normalization, including zero self terms.
    d0=cfg.DistanceScale*distanceSum/(ns*ns); details.characteristicDistance=d0;
    if ~(isfinite(d0) && d0>0), active(1)=false; details.status='zero_spike_spread'; end
end
probabilities=nan(ns,1); fp=0; nearestAP=[];
if active(2), nearestAP=inf(h,nn); end
for first=1:block:ns
    rows=first:min(first+block-1,ns); nrows=numel(rows);
    reference=inf(nrows,1); ws=zeros(nrows,1); wn=ws;
    nearestNoise=inf(nrows,h);
    for start=1:colBlock:nn
        cols=start:min(start+colBlock-1,nn);
        d2=distances(a(rows,:),b(cols,:),na(rows),nb(cols));
        if active(1)
            d=sqrt(d2); next=min(reference,min(d,[],2));
            shrink=exp((next-reference)/d0);
            ws=ws.*shrink; wn=wn.*shrink+sum(exp(-(d-next)/d0),2); reference=next;
        end
        if active(3), nearestNoise=mink([nearestNoise,d2],h,2); end
        if active(2), nearestAP(:,cols)=mink([nearestAP(:,cols);d2],h,1); end %#ok<AGROW> fixed h-by-nn allocation
    end
    if active(1) || active(3)
        apVotes=zeros(nrows,1);
        for start=1:colBlock:ns
            cols=start:min(start+colBlock-1,ns);
            d2=distances(a(rows,:),a(cols,:),na(rows),na(cols));
            d2(rows(:)==cols)=Inf;
            if active(1)
                d=sqrt(d2); next=min(reference,min(d,[],2)); shrink=exp((next-reference)/d0);
                ws=ws.*shrink+sum(exp(-(d-next)/d0),2); wn=wn.*shrink; reference=next;
            end
            if active(3), apVotes=apVotes+sum(isfinite(d2) & d2<=nearestNoise(:,end),2); end
        end
        if active(1), probabilities(rows)=ws./(ws+wn); end
        if active(3), fp=fp+sum(apVotes<h); end
    end
end
if active(1)
    isolation=mean(probabilities); details.TargetAffinity=probabilities; details.Completed(1)=true;
end
if active(3), fpErr=fp/ns; details.FPCount=fp; details.Completed(3)=true; end
if active(2)
    radius=nearestAP(end,:)'; clear nearestAP
    fn=0;
    for first=1:block:nn
        queries=(first:min(first+block-1,nn))'; count=zeros(numel(queries),1);
        % The first 64 references usually resolve distant background events.
        starts=[1,65:colBlock:nn];
        for j=1:numel(starts)
            start=starts(j);
            if j==1, last=min(64,nn); else, last=min(start+colBlock-1,nn); end
            cols=start:last;
            d2=distances(b(queries,:),b(cols,:),nb(queries),nb(cols));
            d2(queries==cols)=Inf;
            details.NoiseDistancesEvaluated=details.NoiseDistancesEvaluated+numel(d2);
            count=count+sum(d2<radius(queries),2);
            keep=count<h; queries=queries(keep); count=count(keep);
            if isempty(queries), break; end
        end
        fn=fn+numel(queries);
    end
    fnErr=fn/(fn+ns); details.FNCount=fn; details.Completed(2)=true;
end
if ~strcmp(details.status,'zero_spike_spread'), details.status='ok'; end
if ~canKNN && ~strcmp(details.status,'zero_spike_spread'), details.status='insufficient_knn_population'; end
end

function d2=distances(x,y,nx,ny)
% Full waveform Euclidean distances: no PCA or feature rescaling.
d2=max(0,nx+ny'-2*x*y');
end

function [T50,best_ms,details]=historical_t50(counts,rate_factor)
%HISTORICAL_T50 Literal complex-cell estimator from the archived phase_code.m.
% counts: phase x cell x 2000, 1 ms histogram bins covering 0.5--2.5 s.
% The historical scale is 16.66 (not exactly 1000/60).
% Return bin numbers (301 means no crossing), WITHOUT the plot's +4 ms.
assert(ndims(counts)==3 && size(counts,3)==2000);
kernel=exp(-((1:21)-11).^2/(2*5^2));
kernel=kernel/sum(kernel);
T50=301*ones(size(counts,1),size(counts,2));
details.bg_hz=zeros(size(T50));details.steady_hz=zeros(size(T50));
for phase=1:size(counts,1)
    for cellnum=1:size(counts,2)
        rate=conv(reshape(double(counts(phase,cellnum,:)),1,[])*rate_factor,kernel,'same');
        bg=mean(rate(1:1000));ss=mean(rate(1200:2000));
        details.bg_hz(phase,cellnum)=bg;
        details.steady_hz(phase,cellnum)=ss;
        threshold=bg+0.5*(ss-bg);
        for t=1:300
            if min(rate(1000+t:1000+t+25))>threshold
                T50(phase,cellnum)=t;
                break;
            end
        end
    end
end
best_ms=min(T50,[],1);
end

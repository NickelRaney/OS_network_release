function validate_result(res,param)
%VALIDATE_RESULT Detect corrupt output without changing accepted spikes.
names={'VE','VI','HE','HI','L6','amb','lgn','spike'};
for k=1:numel(names)
    assert(all(isfinite(res.(names{k})),'all'),'Nonfinite %s output.',names{k});
end
assert(size(res.spike,2)==16090);
for k=1:size(res.spike,2)
    n=res.spike(1,k);
    assert(n>=0 && n==floor(n) && n<size(res.spike,1));
    t=res.spike(2:n+1,k);
    assert(all(t>0 & t<=param.duration) && all(diff(t)>0));
    assert(all(abs(t*1000/param.dt-round(t*1000/param.dt))<1e-8));
    assert(all(res.spike(n+2:end,k)==0));
end
end

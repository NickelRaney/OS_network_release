function [out,lgn]=lgnout(param,kernel,input,onset)
%LGNOUT Original local helper, with an explicit optional onset argument.
% LGN timestamps are milliseconds. The visual drive is sampled every 1 ms.
assert(param.dt==0.1,'OSNetwork:TimeStep','Legacy LGN helper requires dt=0.1 ms.');
c=param.c;
IB=param.IB;
init=0;
lgn=zeros(110,2000);
if nargin<4
    out=LGN3(param,kernel,input);
else
    out=LGN4(param,kernel,input,onset);
end
for i=1:110
    res=NLIF(param,max(c*kron(out(i,:),ones(1,10))+IB,0),init);
    spike=res.spike;
    lgn(i,1:spike(1)+1)=spike(1:spike(1)+1);
end
end

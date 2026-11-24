function [out]=LGN4(param,kernel,input,tot)
 %sudden global input
 %tot means turn on time


pos=param.LGN_pos; % the LGN information
duration=param.duration*1000;
sker=kernel.sker';

nLGN=size(pos,1);
out=[];

sp_on=zeros(1,200);

speed=abs(param.speed);

for t=1:duration
    sp_on(t)=input(1,ceil(t*speed*1000):ceil(t*speed*1000)+600)*sker;
end
sp_on=sp_on/1000;


maxshift=round(max(pos(:,1))/speed);



for i=1:nLGN
    shift=round(pos(i,1)/speed);
    sp=pos(i,2)*[zeros(1,tot*1000),sp_on(maxshift-shift+1:end),sp_on(end)*ones(1,maxshift-shift)];
    out(i,:)=conv_full(sp,kernel,pos(i,2));
end

end





function a=conv_full(b,kernel,i)
    if i==1
        ker=kernel.tker_on;
    else
        ker=kernel.tker_off;
    end
    b=[b(1)*ones(1,length(ker)),b,b(end)*ones(1,length(ker))];
    b=conv(b,ker);
    a=b(length(ker)+1:end-length(ker));

end

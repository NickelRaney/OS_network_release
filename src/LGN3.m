function [out]=LGN3(param,kernel,input)
% In fact we only need to calculate the output of one single on cell and
% off cell. All other cell will get the same output with some delay.

% The input will be the same input as LGN2, plus delay time and its cell type.
% Note that the 'shift' variable is given here now.

% We delete the saccade function here.


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

%3.3333(number of HC)*0.25(size of HC)*sqrt(2)~1.2
% tail=round(1.2/speed);
tail=round(max(pos(:,1))*1000+100);




for i=1:nLGN
    shift=round(pos(i,1)/speed);
    sp=pos(i,2)*[sp_on(1)*ones(1,shift),sp_on,sp_on(end)*ones(1,tail-shift)];
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

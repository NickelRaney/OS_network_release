function [res]=NLIF(param,It,init)

tau_lk=param.tau_leak;
lambda=param.lambda;
s_ex=param.s_ex;
duration=param.duration*1000; %(ms)
v_th=param.v_th;
dt=param.dt;




esi_e=exprnd(1/lambda,ceil(lambda*(duration*1.4)),1); %external spike interval
while min(sum(esi_e))<duration
    esi_e=exprnd(1/lambda,ceil(lambda*(duration*1.4)),1); %external spike interval
end


ex_e=zeros(ceil(duration/dt),1);
%convert esi to the external effect at each time point.



t=esi_e(1);
count=1;
while t<duration
    ind=ceil(t/dt);
    ex_e(ind)=ex_e(ind)+sign(rand(1)-0.5)*exp((t-ind*dt)*tau_lk);
    count=count+1;
    t=t+esi_e(count);
end

v=zeros(ceil(duration/dt),1);
spike=zeros(1,ceil(duration/1000*30));


volt=init;

for step=1:duration/dt
    time = step*dt;
    v(step)=volt;
    volt=volt+(-tau_lk*volt+It(step))*dt+ex_e(step)*s_ex;

    if volt>v_th
        volt=0;
        spike(1)=spike(1)+1;
        spike(spike(1)+1)=time;
    end
end
res.vtra=v;
res.spike=spike;
res.es=ex_e;
end

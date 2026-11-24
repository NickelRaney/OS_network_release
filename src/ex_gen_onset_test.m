function [res]=ex_gen_onset_test(param,cat_e,cat_i,LVwire,LVwire_i,lgn,onset_time)
dt=param.dt;
duration=param.duration*1000;
ne=size(cat_e,1);
ni=size(cat_i,1);
tau_ampa_r = 0.5;
tau_ampa_d = 2;
step=duration/dt;

q_ampa_r=exp(-dt/tau_ampa_r);
q_ampa_d=exp(-dt/tau_ampa_d);

nlgn=110;
lgn_ampa_r_e=zeros(step,nlgn);
lgn_ampa_d_e=zeros(step,nlgn);

lgn_ampa_r_i=zeros(step,nlgn);
lgn_ampa_d_i=zeros(step,nlgn);

i_decay=([ 1.2800    1.1400    1.0700    1.0378    1.0112, 1.00, 0.99, 0.982, 0.975, 0.975, 0.977, 0.980, 0.985, 0.988, 0.990, 0.993, 0.995, 0.998, 0.999, 1.000]-1)*15/28+1;
e_decay = ([1.25, 1.28, 1.28, 1.278, 1.275, 1.268, 1.257, 1.240, 1.221, 1.200, 1.177, 1.153, 1.128, 1.104, 1.080, 1.058, 1.038, 1.022, 1.008, 1.000]-1)*11/25+1;

% Preserve the audited onset arrays; steady input gains do not replace them.
if isfield(param,'onset_gain_e_initial') && param.onset_gain_e_initial~=1.11
    e_decay=(e_decay-1)*(param.onset_gain_e_initial-1)/0.11+1;
end
if isfield(param,'onset_gain_i_initial') && param.onset_gain_i_initial~=1.15
    i_decay=(i_decay-1)*(param.onset_gain_i_initial-1)/0.15+1;
end

for i=1:nlgn
   flag=1;
    cofe=e_decay(1);
    cofi=i_decay(1);

     for j=1:lgn(i,1)
        t=lgn(i,j+1);
        ind=round(t/dt);
        if t>duration
            break;
        end

        if t>onset_time*1000
            if flag==1
            nstart=j;
            flag=0;
            end

            if j-nstart+1<=20
                cofe=e_decay(j-nstart+1);
            else
                cofe=1;
            end

            if j-nstart+1<=20
                cofi=i_decay(j-nstart+1);
            else
                cofi=1;
            end
        end





        lgn_ampa_r_e(ind,i)=lgn_ampa_r_e(ind,i)+cofe*1/(tau_ampa_d-tau_ampa_r);
        lgn_ampa_d_e(ind,i)=lgn_ampa_d_e(ind,i)+cofe*1/(tau_ampa_d-tau_ampa_r);

        lgn_ampa_r_i(ind,i)=lgn_ampa_r_i(ind,i)+cofi*1/(tau_ampa_d-tau_ampa_r);
        lgn_ampa_d_i(ind,i)=lgn_ampa_d_i(ind,i)+cofi*1/(tau_ampa_d-tau_ampa_r);
     end
end

for i=2:size(lgn_ampa_r_e,1)
    lgn_ampa_r_e(i,:)=lgn_ampa_r_e(i,:)+lgn_ampa_r_e(i-1,:)*q_ampa_r;
    lgn_ampa_d_e(i,:)=lgn_ampa_d_e(i,:)+lgn_ampa_d_e(i-1,:)*q_ampa_d;

    lgn_ampa_r_i(i,:)=lgn_ampa_r_i(i,:)+lgn_ampa_r_i(i-1,:)*q_ampa_r;
    lgn_ampa_d_i(i,:)=lgn_ampa_d_i(i,:)+lgn_ampa_d_i(i-1,:)*q_ampa_d;
end

lgn_net_e=lgn_ampa_d_e-lgn_ampa_r_e;
lgn_net_i=lgn_ampa_d_i-lgn_ampa_r_i;

I = []; J = []; V = [];
for i = 1:ne
    if cat_e(i,2) ~= 0
        lgn_indices = LVwire(i, 1:cat_e(i,2));
        I = [I, lgn_indices];
        J = [J, repmat(i, 1, length(lgn_indices))];
        V = [V, ones(1, length(lgn_indices))];
    end
end
conn_matrix_e = sparse(I, J, V, nlgn, ne);

I = []; J = []; V = [];
for i = 1:ni

        lgn_indices = LVwire_i(i, 1:cat_i(i,1));
        I = [I, lgn_indices];
        J = [J, repmat(i, 1, length(lgn_indices))];
        V = [V, ones(1, length(lgn_indices))];

end
conn_matrix_i = sparse(I, J, V, nlgn, ni);

res.exe=lgn_net_e;
res.exi=lgn_net_i;
res.conn_matrix_e=conn_matrix_e;
res.conn_matrix_i=conn_matrix_i;
end

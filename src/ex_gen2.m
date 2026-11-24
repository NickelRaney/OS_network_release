function [res]=ex_gen2(param,cat_e,cat_i,LVwire,LVwire_i,lgn)
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

for i=1:nlgn

    cofe=1.11;
    if isfield(param,'lgn_gain_e_initial'), cofe=param.lgn_gain_e_initial; end
    cofi=1.15;
    if isfield(param,'lgn_gain_i_initial'), cofi=param.lgn_gain_i_initial; end

    flag_e=0;
    flag1_e=0;

    flag_i=0;
    flag1_i=0;

    for j=1:lgn(i,1)
        t=lgn(i,j+1);
        ind=round(t/dt);
        if t>duration
            break;
        end

        if j>=12 && lgn(i,j+1)-lgn(i,j-10)<240 && flag_e==0
            cofe=1.12;
            flag_e=1;
        end

        if j>=20 && lgn(i,j+1)-lgn(i,j-18)<400 && flag1_e==0
            cofe=1;
            flag1_e=1;

        end

        if j>=7 && lgn(i,j+1)-lgn(i,j-5)<140 && flag_i==0
            cofi=1.05;
            flag_i=1;
        end

        if j>=12 && lgn(i,j+1)-lgn(i,j-10)<240 && flag1_i==0
            cofi=1;
            flag1_i=1;

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

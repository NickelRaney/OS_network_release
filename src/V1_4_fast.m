function [res]=V1_4_fast(param,init,ex)
%V1_4_FAST Core MATLAB LIF network simulator for macaque V1 layer 4Ca.
% Optional fields param.f1_ee_scale and param.f1_ei_scale control F1
% depression scaling (defaults: 2.25 and 1.5 for Disordered-LGN).
ee_scale=2.25; ei_scale=1.5;
if isfield(param,'f1_ee_scale'), ee_scale=param.f1_ee_scale; end
if isfield(param,'f1_ei_scale'), ei_scale=param.f1_ei_scale; end


e_decay = ([1.2 1.2 1.2 1.20 1.19 1.16 1.12 1.07 1 0.96 0.94 0.92]-0.92)*ee_scale/3+0.92;
i_decay = ([1.1 1.02 0.98 0.96 0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95]-0.95)*ei_scale/1.5+0.95;
ei_decay=[1.05 1 0.98];
see_dep_recover_time=133;
sei_dep_recover_time=50;

% -------- unpack params (unchanged) --------
cmee=param.cmee;
cmei=param.cmei; % Note: logic implies this is I->E based on usage below
cmie=param.cmie; % Note: logic implies this is E->I
cmii=param.cmii;
map_ei=param.map_ei;
L6tfm=param.L6tfm;
seetfm=param.seetfm;
l6cof=param.l6cof;
conn_matrix_e=ex.conn_matrix_e;
conn_matrix_i=ex.conn_matrix_i;
see = param.see;
sie = param.sie*ones(13689,1);
sei = param.sei*ones(2401,1);
sei_list = param.sei_list;
selgn=param.selgn;
silgn=param.silgn;
segabab = param.segabab;
seifac = param.seifac;
sii = param.sii;
sel6 = param.sel6;
sil6 = param.sil6;
l6facmax=param.l6facmax;
l6facier=param.l6facier;
ndep_e=param.ndep_e;
ndep_i=param.ndep_i;
pdep_e=param.pdep_e;
pdep_i=param.pdep_i;
sincre_e=param.sincre_e;
sbasic_e=param.sie;
supb_e=sbasic_e;
slob_e=sbasic_e-5*sincre_e;
sincre_i=param.sincre_i;
sbasic_i=param.sei;
supb_i=sbasic_i;
slob_i=param.slob_i;

% -------- Decay Tables --------


tau_ampa_r = 0.5; tau_ampa_d = 2;
tau_nmda_r = 2;   tau_nmda_d = 80;
tau_gaba_r = 0.5; tau_gaba_d = 5;
tau_gaba_b_r=3;   tau_gaba_b_d=100;
tau_thrd_e = 5.05;
dt = param.dt;
ref_e = param.tau_re;
ref_i = param.tau_ri;
duration = param.duration*1000;
onset=param.onset*10000;
ne=size(cmee,1);
ni=size(cmei,1);
exe=ex.exe;
exi=ex.exi;

% Pre-compute sparse synaptic weight matrices
% E->E and E->I: linear updates via sparse matrix multiplication
cme_EE = sparse(cmee * see);
cme_EI = sparse(cmie .* sie);

% I->I: linear update
cmi_II = sparse(cmii * sii);

% I->E: nonlinear update max(W * boost - seifac, 0).
% Transpose cmei (Ni x Ne) to cmi_IE_T (Ne x Ni) for fast column slicing.
raw_cmi_IE = sparse(cmei .* sei);
cmi_IE_T   = raw_cmi_IE';

% Kernel & Mask (unchanged)
kernel=exp(-((1:21)'-11).^2/240)*exp(-((1:21)-11).^2/240);
kernel=kernel/sum(sum(kernel));
mask.weight_low = conv2(m_large(param.mask,10),kernel,'valid');
mask.weight_high= conv2(m_large(~param.mask,10),kernel,'valid');
mask.sign=param.mask;

% Init State
if isempty(init)
    ve=zeros(1,ne); vi=zeros(1,ni);
    hee=zeros(4,ne); hie=zeros(4,ni);
    hei=zeros(4,ne); hii=zeros(2,ni);
    hee_sum=zeros(1,ne); hie_sum=zeros(1,ni);
    hei_sum=zeros(1,ne); hii_sum=zeros(1,ni);
else
    ve=init.ve; vi=init.vi;
    hee=0.03*rand(4,ne); hie=0.03*rand(4,ni);
    hei=0.05*rand(4,ne); hii=0.05*rand(2,ni);
    hee_sum=0.8*(hee(2,:)-hee(1,:))+0.2*(hee(4,:)-hee(3,:));
    hie_sum=0.8*(hie(2,:)-hie(1,:))+0.2*(hie(4,:)-hie(3,:));
    hei_sum=hei(2,:)-hei(1,:)+hei(4,:)-hei(3,:);
    hii_sum=hii(2,:)-hii(1,:);
end
refc_e=zeros(1,ne); refc_i=zeros(1,ni);
trhd_e=zeros(1,ne);
spike_e=zeros(1000,ne); spike_i=zeros(1000,ni);

% Recording
VE = zeros(ceil(duration/dt/10)+1,ne);
VI = zeros(ceil(duration/dt/10)+1,ni);
HE = zeros(ceil(duration/dt/10)+1,ne+ni);
HI = zeros(ceil(duration/dt/10)+1,ne+ni);
amb = zeros(ceil(duration/dt/10)+1,ne);
lgn = zeros(ceil(duration/dt/10)+1,ne);
L6 = zeros(ceil(duration/dt/10)+1,ne+ni);
spikecount_short=zeros(70,ne);
spikecount_gabab=zeros(300,ne);
spikecount_dep_e=zeros(duration/100+1,13689);
spikecount_dep_i=zeros(duration/100+1,2401);
scind=1; gbind=1; dep_ind=0;

% Noise
mean_e=0.0124*0.8; theta=0.4*0.2; sigma_e=0.1*sqrt(0.8*0.2)*mean_e;
mean_i=0.0124*0.75; sigma_i=0.1*sqrt(0.8*0.2)*mean_i;
exe_amb=rand(1,ne)*2*mean_e; exi_amb=rand(1,ni)*2*mean_i;
see_rec=ones(1,13689); sie_rec=ones(1,13689); sei_rec=ones(1,2401);
see_check=see_rec; sie_check=sie_rec; sei_check=sei_rec;

% Constants
dec_ar=exp(-dt/tau_ampa_r); dec_ad=exp(-dt/tau_ampa_d);
dec_nr=exp(-dt/tau_nmda_r); dec_nd=exp(-dt/tau_nmda_d);
dec_gr=exp(-dt/tau_gaba_r); dec_gd=exp(-dt/tau_gaba_d);
dec_gbr=exp(-dt/tau_gaba_b_r); dec_gbd=exp(-dt/tau_gaba_b_d);
dec_thr=exp(-dt/tau_thrd_e);
inv_amp = 1/(tau_ampa_d - tau_ampa_r);
amp_to_nmda = (tau_ampa_d - tau_ampa_r)/(tau_nmda_d - tau_nmda_r);
inv_gab = 1/(tau_gaba_d - tau_gaba_r);

% ===================== MAIN LOOP =====================
for step=2:duration/dt
    t=dt*step;

    % 1. Ambient & Refractory
    exe_amb = exe_amb - dt*theta*(exe_amb-mean_e) + sigma_e*sqrt(dt)*randn(1,ne);
    exi_amb = exi_amb - dt*theta*(exi_amb-mean_i) + sigma_i*sqrt(dt)*randn(1,ni);

    rind_e = (refc_e < 1e-7); rind_i = (refc_i < 1e-7);
    refc_e(~rind_e)=refc_e(~rind_e)-dt; refc_i(~rind_i)=refc_i(~rind_i)-dt;

    % 2. Membrane Updates
    exe_row = exe(step,:) * conn_matrix_e * selgn + exe_amb;
    exi_row = exi(step,:) * conn_matrix_i * silgn + exi_amb;

    ve = ve - (ve/20 + (exe_row+hee_sum).*(ve - 14/3) + hei_sum.*(ve + 2/3)) .* rind_e * dt;
    vi = vi - (vi/15 + (exi_row+hie_sum).*(vi - 14/3) + hii_sum.*(vi + 2/3)) .* rind_i * dt;

    % 3. Decay
    hee(1,:)=hee(1,:).*dec_ar; hee(2,:)=hee(2,:).*dec_ad;
    hee(3,:)=hee(3,:).*dec_nr; hee(4,:)=hee(4,:).*dec_nd;
    hie(1,:)=hie(1,:).*dec_ar; hie(2,:)=hie(2,:).*dec_ad;
    hie(3,:)=hie(3,:).*dec_nr; hie(4,:)=hie(4,:).*dec_nd;
    hei(1,:)=hei(1,:).*dec_gr; hei(2,:)=hei(2,:).*dec_gd;
    hei(3,:)=hei(3,:).*dec_gbr; hei(4,:)=hei(4,:).*dec_gbd;
    hii(1,:)=hii(1,:).*dec_gr; hii(2,:)=hii(2,:).*dec_gd;
    trhd_e = trhd_e .* dec_thr;

    % 4. Spikes
    sind_e = find(ve > 1+trhd_e);
    sind_i = find(vi > 1);
    spikecount_e = numel(sind_e);
    spikecount_i = numel(sind_i);
    spikecount_short(scind,:)=0;
    spikecount_gabab(gbind,:)=0;

    % === Excitatory Spikes (Vectorized Linear) ===
    if spikecount_e>0
        ve(sind_e)=0;
        refc_e(sind_e)=refc_e(sind_e)+ref_e(sind_e);
        boosts_EE = sparse(1, ne);
        boosts_EI = sparse(1, ne);

        for ii=1:spikecount_e
            ind=sind_e(ii);
            % --- Boost Calculation (Standard) ---
            spc=spike_e(2:spike_e(1,ind)+1,ind);
            prevspc_200= round((spike_e(1,ind)+sum(sign(spc - t + 199.99)))/2);
            prevspc = round((spike_e(1,ind)+sum(sign(spc - t + 99.99)))/2);

            if spike_e(1,ind)<ndep_e || rand(1)<pdep_e^(prevspc-ndep_e+1)
                t_pre=spike_e(spike_e(1,ind)+1,ind);
                if see_check(ind)>0.92
                    if t-t_pre<50, see_rec(ind)=see_rec(ind)+1;
                    elseif t-t_pre>see_dep_recover_time, see_rec(ind)=see_rec(ind)-1; end
                    see_rec(ind)=min(max(see_rec(ind),1),12);
                    see_boost=e_decay(see_rec(ind)); sie_boost=i_decay(see_rec(ind));
                elseif see_check(ind)<0.92
                    see_boost=seetfm(prevspc_200+1); sie_boost=seetfm(prevspc_200+1)*sie(ind)/0.0176;
                else
                    if seetfm(prevspc_200+1)<0.92
                        see_boost=seetfm(prevspc_200+1); sie_boost=seetfm(prevspc_200+1)*sie(ind)/0.0176;
                    elseif t-t_pre>see_dep_recover_time
                        see_rec(ind)=max(see_rec(ind)-1,1);
                        see_boost=e_decay(see_rec(ind)); sie_boost=i_decay(see_rec(ind));
                    else
                        see_boost=0.92; sie_boost=0.92;
                    end
                end
                see_check(ind)=see_boost; sie_check(ind)=sie_boost;

                boosts_EE(ind) = boosts_EE(ind) + see_boost;
                boosts_EI(ind) = boosts_EI(ind) + sie_boost;

                trhd_e(ind)=0.5;
                spike_e(1,ind)=spike_e(1,ind)+1;
                spike_e(spike_e(1,ind)+1,ind)=t;
                spikecount_short(scind,ind)=1;
            end
        end

        if nnz(boosts_EE) > 0 || nnz(boosts_EI) > 0
            C1 = full(boosts_EE * cme_EE) * inv_amp;
            C2 = full(boosts_EI * cme_EI) * inv_amp;
            hee(1,:)=hee(1,:)+C1; hee(2,:)=hee(2,:)+C1;
            hee(3,:)=hee(3,:)+C1*amp_to_nmda; hee(4,:)=hee(4,:)+C1*amp_to_nmda;
            hie(1,:)=hie(1,:)+C2; hie(2,:)=hie(2,:)+C2;
            hie(3,:)=hie(3,:)+C2*amp_to_nmda; hie(4,:)=hie(4,:)+C2*amp_to_nmda;
        end
    end

    % === Inhibitory Spikes (Robust I-to-E Logic) ===
    if spikecount_i>0
        vi(sind_i)=0;
        refc_i(sind_i)=refc_i(sind_i)+ref_i(sind_i);

        cI_accumulator_E = sparse(1, ne);
        boosts_II = sparse(1, ni);

        for ii=1:spikecount_i
            ind=sind_i(ii);
            prevspc = round((spike_i(1,ind)+sum(sign(spike_i(2:spike_i(1,ind)+1,ind)-t+99.99)))/2);
            if spike_i(1,ind)<ndep_i || rand(1)<pdep_i^(prevspc-ndep_i+1)
                t_pre=spike_i(spike_i(1,ind)+1,ind);

                % --- Boost Calculation ---
                if sei_check(ind)>0.98
                    if t-t_pre<25, sei_rec(ind)=sei_rec(ind)+1;
                    elseif t-t_pre>sei_dep_recover_time, sei_rec(ind)=sei_rec(ind)-1; end
                    sei_rec(ind)=min(max(sei_rec(ind),1),3); sei_boost=ei_decay(sei_rec(ind));
                elseif sei_check(ind)<0.98
                    sei_boost=sei(ind)/0.051;
                else
                    if sei(ind)/0.051<0.98, sei_boost=sei(ind)/0.051;
                    elseif t-t_pre>sei_dep_recover_time
                        sei_rec(ind)=max(sei_rec(ind)-1,1); sei_boost=ei_decay(sei_rec(ind));
                    else
                        sei_boost=0.98;
                    end
                end
                sei_check(ind)=sei_boost;

                % I-to-E nonlinear synaptic update: max(W * sei_boost - seifac, 0)
                w_col = cmi_IE_T(:, ind);

                if nnz(w_col) > 0
                    [target_idxs, ~, w_vals] = find(w_col);

                    if isscalar(seifac)
                        eff_seifac = seifac;
                    else
                        eff_seifac = seifac(target_idxs);
                        if size(eff_seifac, 2) > 1, eff_seifac = eff_seifac'; end
                    end

                    vals = w_vals * sei_boost - eff_seifac;
                    keep_mask = vals > 0;

                    if any(keep_mask)
                        cI_accumulator_E = cI_accumulator_E + ...
                            sparse(1, target_idxs(keep_mask), vals(keep_mask), 1, ne);
                    end
                end

                % --- I to I Update (Linear) ---
                boosts_II(ind) = boosts_II(ind) + 1;

                spike_i(1,ind)=spike_i(1,ind)+1;
                spike_i(spike_i(1,ind)+1,ind)=t;
                spikecount_gabab(gbind,:)=spikecount_gabab(gbind,:)+cmei(ind,:);
            end
        end

        if nnz(boosts_II) > 0 || nnz(cI_accumulator_E) > 0
            C1 = full(cI_accumulator_E) * inv_gab;
            C2 = full(boosts_II * cmi_II) * inv_gab;
            hei(1,:)=hei(1,:)+C1; hei(2,:)=hei(2,:)+C1;
            hii(1,:)=hii(1,:)+C2; hii(2,:)=hii(2,:)+C2;
        end
    end

    % Circular Pointers & GABAB
    scind=scind+1; if scind>70, scind=1; end
    gbind=gbind+1; if gbind>300, gbind=1; end
    gbsum=sum(spikecount_gabab,1);
    if max(gbsum)>265
        gb_sp=(gbsum>265);
        cb = gb_sp*segabab * (1/(tau_gaba_b_d - tau_gaba_b_r));
        spikecount_gabab(:,gb_sp)=0;
        hei(3,:)=hei(3,:)+cb; hei(4,:)=hei(4,:)+cb;
    end

    % L6 feedback & facilitation
    l6uptind=mod(step-2,10);
    if l6uptind==0
        l6step=round((step-2)/10);
        [L6_spike_e,L6_spike_i]=fL6(spikecount_short,mask,kernel,map_ei,L6tfm,l6cof);
        L6(l6step+1,:)=sum([L6_spike_e,L6_spike_i],1);
        if l6step+1>=100
            fac=1 + l6facmax ./ (exp(-(sum(L6(l6step-98:l6step+1,1:13689),1)/3.5 - 40)/4.5) + 1);
            fac(fac<1.05)=-0.2625+1.25*fac(fac<1.05);
            fac=max(fac,1);
            L6_spike_e=L6_spike_e.*fac;
            L6_spike_i=L6_spike_i.*fac(map_ei)*l6facier;
        end
    end

    c1=sel6.*L6_spike_e(l6uptind+1,:) * inv_amp;
    c2=sil6.*L6_spike_i(l6uptind+1,:) * inv_amp;
    hee(1,:)=hee(1,:)+c1; hee(2,:)=hee(2,:)+c1;
    hee(3,:)=hee(3,:)+c1*amp_to_nmda; hee(4,:)=hee(4,:)+c1*amp_to_nmda;
    hie(1,:)=hie(1,:)+c2; hie(2,:)=hie(2,:)+c2;
    hie(3,:)=hie(3,:)+c2*amp_to_nmda; hie(4,:)=hie(4,:)+c2*amp_to_nmda;

    hee_sum=0.8*(hee(2,:)-hee(1,:))+0.2*(hee(4,:)-hee(3,:));
    hie_sum=0.8*(hie(2,:)-hie(1,:))+0.2*(hie(4,:)-hie(3,:));
    hei_sum=hei(2,:)-hei(1,:)+hei(4,:)-hei(3,:);
    hii_sum=hii(2,:)-hii(1,:);

    if mod(step,1000)==0
        dep_ind=dep_ind+1;
        spikecount_dep_e(dep_ind,:)=spike_e(1,:);
        spikecount_dep_i(dep_ind,:)=spike_i(1,:);
        if dep_ind>=3
            diff_e=spikecount_dep_e(dep_ind,:)-spikecount_dep_e(dep_ind-2,:);
            sie=sie+sign(sbasic_e-(diff_e'-6)*sincre_e-sie)*sincre_e;
            sie=min(sie,supb_e); sie=max(sie,slob_e);
            diff_i=spikecount_dep_i(dep_ind,:)-spikecount_dep_i(dep_ind-2,:);
            sei=sei_list(diff_i+1);
        end
    end

    if mod(step,10)==0
        idx = step/10;
        VE(idx,:)=ve(:); VI(idx,:)=vi(:);
        amb(idx,:)=exe_amb; lgn(idx,:)=exe_row-exe_amb;
        HE(idx,:)=[hee_sum,hie_sum]; HI(idx,:)=[hei_sum,hii_sum];
    end
end

res.VE = VE; res.VI = VI; res.HE = HE; res.HI = HI;
res.L6 = L6; res.amb = amb; res.lgn = lgn;
res.spike=[spike_e,spike_i]; res.spike(2:end,:)=res.spike(2:end,:)/1000;
end

function [L6_spike_e,L6_spike_i]=fL6(spike,mask,kernel,map_ei,L6tfm,l6cof)
fr=reshape(sum(spike,1),117,117)/1.4*200;
fr_low = fr.*sign(mask.sign);
fr_high= fr.*sign(~mask.sign);
l6 = 0.89*conv2(m_large(fr_low,10),kernel,'valid')./mask.weight_low ...
   + 0.11*conv2(m_large(fr_high,10),kernel,'valid')./mask.weight_high;
l6=reshape(l6,1,13689);
l6e=L6tfm(round(l6*10)+1);
l6i=l6e(map_ei);
fre=l6e.*l6cof/10000;
fri=l6i*100/10000;
L6_spike_e = rand(10,13689) < fre;
L6_spike_i = rand(10,2401)  < fri;
end

function out = m_large(m,n)
out = [m(end-(n-1):end,end-(n-1):end), m(end-(n-1):end,:), m(end-(n-1):end,1:n); ...
       m(:,end-(n-1):end),               m(:,:),             m(:,1:n); ...
       m(1:n,end-(n-1):end),             m(1:n,:),           m(1:n,1:n)];
end

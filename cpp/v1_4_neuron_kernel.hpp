// Exact double-precision neuron kernels. Compile with /fp:strict /arch:AVX2
// (or -ffp-contract=off -fno-fast-math -mavx2). Components use SoA layout:
// state[component * n + neuron]. No reduction or operation is reassociated.
#pragma once
#include <cstddef>
#include <immintrin.h>
#include <vector>

namespace v1fast {
using Size = std::size_t;

// Appends in the same ascending-neuron order as the scalar reference.
// Threshold is optional (I neurons use 1), and callers clear spiking first.
inline void detect(Size n, double* voltage, double* refractory,
    const double* refractoryPeriod, const double* threshold,
    std::vector<Size>& spiking) {
    const __m256d one=_mm256_set1_pd(1.0);
    Size j=0;
    for(;j+4<=n;j+=4) {
        const __m256d limit=threshold?_mm256_add_pd(one,_mm256_loadu_pd(threshold+j)):one;
        const int mask=_mm256_movemask_pd(_mm256_cmp_pd(_mm256_loadu_pd(voltage+j),limit,_CMP_GT_OQ));
        if(mask) for(Size lane=0;lane<4;++lane) if(mask&(1<<lane)) {
            const Size index=j+lane;
            spiking.push_back(index);
            voltage[index]=0;
            refractory[index]=refractory[index]+refractoryPeriod[index];
        }
    }
    for(;j<n;++j) if(voltage[j]>(threshold?1+threshold[j]:1)) {
        spiking.push_back(j);
        voltage[j]=0;
        refractory[j]=refractory[j]+refractoryPeriod[j];
    }
}

inline void advance(Size n, double* voltage, double* refractory,
    double* ambient, double* external, const double* lgnAlreadyScaled,
    const double* normal, const double* heSum, const double* hiSum,
    double* threshold, double* excitation, double* inhibition,
    Size inhibitionComponents, double dt, double theta, double mean,
    double noise, double membraneTau, double vE, double vI,
    const double* decay) {
    const double dtTheta=dt*theta;
    const __m256d dT=_mm256_set1_pd(dt), dTheta=_mm256_set1_pd(dtTheta);
    const __m256d dMean=_mm256_set1_pd(mean), dNoise=_mm256_set1_pd(noise);
    const __m256d dTau=_mm256_set1_pd(membraneTau), dVE=_mm256_set1_pd(vE);
    const __m256d dVI=_mm256_set1_pd(vI), epsilon=_mm256_set1_pd(1e-7);
    const __m256d one=_mm256_set1_pd(1.0);
    __m256d dec[9];
    for(Size k=0;k<9;++k) dec[k]=_mm256_set1_pd(decay[k]);
    Size j=0;
    for(;j+4<=n;j+=4) {
        const __m256d oldAmbient=_mm256_loadu_pd(ambient+j);
        const __m256d am=_mm256_add_pd(
            _mm256_sub_pd(oldAmbient,_mm256_mul_pd(dTheta,_mm256_sub_pd(oldAmbient,dMean))),
            _mm256_mul_pd(dNoise,_mm256_loadu_pd(normal+j)));
        _mm256_storeu_pd(ambient+j,am);
        const __m256d ex=_mm256_add_pd(_mm256_loadu_pd(lgnAlreadyScaled+j),am);
        _mm256_storeu_pd(external+j,ex);
        const __m256d ref=_mm256_loadu_pd(refractory+j);
        const __m256d active=_mm256_cmp_pd(ref,epsilon,_CMP_LT_OQ);
        _mm256_storeu_pd(refractory+j,_mm256_blendv_pd(_mm256_sub_pd(ref,dT),ref,active));
        const __m256d v=_mm256_loadu_pd(voltage+j);
        const __m256d exc=_mm256_mul_pd(_mm256_add_pd(ex,_mm256_loadu_pd(heSum+j)),_mm256_sub_pd(v,dVE));
        const __m256d inh=_mm256_mul_pd(_mm256_loadu_pd(hiSum+j),_mm256_add_pd(v,dVI));
        const __m256d change=_mm256_mul_pd(_mm256_mul_pd(
            _mm256_add_pd(_mm256_add_pd(_mm256_div_pd(v,dTau),exc),inh),
            _mm256_and_pd(active,one)),dT);
        _mm256_storeu_pd(voltage+j,_mm256_sub_pd(v,change));
        for(Size k=0;k<4;++k)
            _mm256_storeu_pd(excitation+k*n+j,_mm256_mul_pd(_mm256_loadu_pd(excitation+k*n+j),dec[k]));
        for(Size k=0;k<inhibitionComponents;++k)
            _mm256_storeu_pd(inhibition+k*n+j,_mm256_mul_pd(_mm256_loadu_pd(inhibition+k*n+j),dec[4+k]));
        if(threshold)
            _mm256_storeu_pd(threshold+j,_mm256_mul_pd(_mm256_loadu_pd(threshold+j),dec[8]));
    }
    for(;j<n;++j) {
        ambient[j]=ambient[j]-dtTheta*(ambient[j]-mean)+noise*normal[j];
        external[j]=lgnAlreadyScaled[j]+ambient[j];
        const bool active=refractory[j]<1e-7;
        if(!active) refractory[j]=refractory[j]-dt;
        voltage[j]=voltage[j]-(voltage[j]/membraneTau+(external[j]+heSum[j])*(voltage[j]-vE)+hiSum[j]*(voltage[j]+vI))*static_cast<double>(active)*dt;
        for(Size k=0;k<4;++k) excitation[k*n+j]=excitation[k*n+j]*decay[k];
        for(Size k=0;k<inhibitionComponents;++k) inhibition[k*n+j]=inhibition[k*n+j]*decay[4+k];
        if(threshold) threshold[j]=threshold[j]*decay[8];
    }
}

inline void addExcitation(Size n, double* excitation, const double* input,
    double invamp, double ampnmda) {
    const __m256d inv=_mm256_set1_pd(invamp), nmda=_mm256_set1_pd(ampnmda);
    Size j=0;
    for(;j+4<=n;j+=4) {
        const __m256d a=_mm256_mul_pd(_mm256_loadu_pd(input+j),inv);
        const __m256d b=_mm256_mul_pd(a,nmda);
        for(Size k=0;k<4;++k)
            _mm256_storeu_pd(excitation+k*n+j,_mm256_add_pd(_mm256_loadu_pd(excitation+k*n+j),k<2?a:b));
    }
    for(;j<n;++j) {
        const double a=input[j]*invamp;
        excitation[j]+=a; excitation[n+j]+=a;
        excitation[2*n+j]+=a*ampnmda; excitation[3*n+j]+=a*ampnmda;
    }
}

inline void addInhibition(Size n, double* inhibition, const double* input,
    double invgab) {
    const __m256d inv=_mm256_set1_pd(invgab);
    Size j=0;
    for(;j+4<=n;j+=4) {
        const __m256d a=_mm256_mul_pd(_mm256_loadu_pd(input+j),inv);
        _mm256_storeu_pd(inhibition+j,_mm256_add_pd(_mm256_loadu_pd(inhibition+j),a));
        _mm256_storeu_pd(inhibition+n+j,_mm256_add_pd(_mm256_loadu_pd(inhibition+n+j),a));
    }
    for(;j<n;++j) {
        const double a=input[j]*invgab;
        inhibition[j]+=a; inhibition[n+j]+=a;
    }
}

// l6Inputs is contiguous per neuron. scaleCount==1 broadcasts scale[0].
inline void finishL6(Size n, double* excitation, const double* inhibition,
    Size inhibitionComponents, double* heSum, double* hiSum,
    const double* l6Inputs, const double* scale, Size scaleCount,
    double invamp, double ampnmda) {
    const __m256d inv=_mm256_set1_pd(invamp), nmda=_mm256_set1_pd(ampnmda);
    const __m256d c08=_mm256_set1_pd(0.8), c02=_mm256_set1_pd(0.2);
    const __m256d scalarScale=_mm256_set1_pd(scale[0]);
    Size j=0;
    for(;j+4<=n;j+=4) {
        const __m256d s=scaleCount==1?scalarScale:_mm256_loadu_pd(scale+j);
        const __m256d a=_mm256_mul_pd(_mm256_mul_pd(s,_mm256_loadu_pd(l6Inputs+j)),inv);
        const __m256d b=_mm256_mul_pd(a,nmda);
        __m256d h[4];
        for(Size k=0;k<4;++k) {
            h[k]=_mm256_add_pd(_mm256_loadu_pd(excitation+k*n+j),k<2?a:b);
            _mm256_storeu_pd(excitation+k*n+j,h[k]);
        }
        _mm256_storeu_pd(heSum+j,_mm256_add_pd(_mm256_mul_pd(c08,_mm256_sub_pd(h[1],h[0])),_mm256_mul_pd(c02,_mm256_sub_pd(h[3],h[2]))));
        __m256d hi=_mm256_sub_pd(_mm256_loadu_pd(inhibition+n+j),_mm256_loadu_pd(inhibition+j));
        if(inhibitionComponents==4)
            hi=_mm256_sub_pd(_mm256_add_pd(hi,_mm256_loadu_pd(inhibition+3*n+j)),_mm256_loadu_pd(inhibition+2*n+j));
        _mm256_storeu_pd(hiSum+j,hi);
    }
    for(;j<n;++j) {
        const double a=scale[scaleCount==1?0:j]*l6Inputs[j]*invamp;
        excitation[j]+=a; excitation[n+j]+=a;
        excitation[2*n+j]+=a*ampnmda; excitation[3*n+j]+=a*ampnmda;
        heSum[j]=0.8*(excitation[n+j]-excitation[j])+0.2*(excitation[3*n+j]-excitation[2*n+j]);
        hiSum[j]=inhibition[n+j]-inhibition[j];
        if(inhibitionComponents==4) hiSum[j]=hiSum[j]+inhibition[3*n+j]-inhibition[2*n+j];
    }
}
} // namespace v1fast

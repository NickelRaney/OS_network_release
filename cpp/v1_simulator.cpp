// Standalone V1_4 model. No MATLAB headers, libraries, or runtime calls.
#include "v1_simulator.hpp"
#include "v1_rng.hpp"
#include "v1_l6.hpp"
#include "v1_4_neuron_kernel.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace v1standalone {
namespace {
using Index=std::size_t;
using Vec=std::vector<double>;
using Node=mat5::Node;
void require(bool ok,const char* message) {if(!ok) throw std::runtime_error(message);}
const Node* field(const Node* s,const char* name) {return &s->field(name);}
const double* doubles(const Node* a) {require(a->bytes.empty(),"Expected double numeric data.");return a->values.data();}
double value(const Node* s,const char* name) {return s->field(name).scalar();}
Vec numeric(const Node& a) {Vec v(a.numel());for(Index j=0;j<v.size();++j)v[j]=a.doubleAt(j);return v;}
Vec vectorField(const Node* s,const char* name,Index expected) {
    const Node& a=s->field(name);require(a.numel()>=expected,"Parameter vector too short.");
    Vec v(expected);for(Index j=0;j<expected;++j)v[j]=a.doubleAt(j);return v;
}
struct Broadcast {
    Vec storage;const double* p;Index n;
    Broadcast(const Node* s,const char* name,Index expected):storage(numeric(s->field(name))) {
        p=storage.data();n=storage.size();require(n==1 || n==expected,"Invalid broadcast parameter shape.");
    }
    double operator[](Index j) const {return p[n==1?0:j];}
};
struct Sparse {
    std::vector<std::uint16_t> row;
    std::vector<std::uint32_t> col;
    Vec weights;
    bool uniform=true;
    double common=0;
    Index rows=0,cols=0;
    Sparse()=default;
    explicit Sparse(const Node* a) {
        require(a->dims.size()==2 && !a->jc.empty(),"Expected a sparse matrix.");
        rows=a->dims[0];cols=a->dims[1];require(rows<=65536,"Too many graph rows.");
        row.reserve(a->ir.size());for(auto index:a->ir){require(index<rows,"Sparse row index out of range.");row.push_back(static_cast<std::uint16_t>(index));}
        col=a->jc;weights=a->values;compact();
    }
    void compact() {
        if(!weights.empty()) common=weights[0];
        for(double w:weights)if(w!=common){uniform=false;break;}
        if(uniform){Vec empty;weights.swap(empty);}
    }
    static Sparse transpose(const Node& a,double scale,bool binary=false) {
        require(a.dims.size()==2,"Connectivity must be a matrix.");
        Sparse g;g.rows=a.dims[1];g.cols=a.dims[0];
        require(g.rows<=65536,"Too many graph targets.");
        std::vector<std::uint32_t> degree(g.cols,0);
        auto visit=[&](auto&& f) {
            if(!a.jc.empty()) {
                for(Index target=0;target<g.rows;++target)
                    for(Index k=a.jc[target];k<a.jc[target+1];++k)f(a.ir[k],target,a.values[k]);
            } else if(!a.bytes.empty()) {
                for(Index target=0;target<g.rows;++target)for(Index source=0;source<g.cols;++source)
                    if(a.bytes[target*g.cols+source])f(source,target,1.0);
            } else {
                for(Index target=0;target<g.rows;++target)for(Index source=0;source<g.cols;++source) {
                    const double w=a.values[target*g.cols+source];if(w!=0)f(source,target,w);
                }
            }
        };
        visit([&](Index source,Index,double w){
            require(std::isfinite(w) && (!binary || w==1),"GABAB requires binary connectivity.");
            if(w*scale!=0)++degree[source];
        });
        g.col.resize(g.cols+1);std::uint64_t total=0;
        for(Index j=0;j<g.cols;++j){g.col[j]=static_cast<std::uint32_t>(total);total+=degree[j];}
        require(total<=UINT32_MAX,"Too many graph edges.");g.col[g.cols]=static_cast<std::uint32_t>(total);
        g.row.resize(static_cast<Index>(total));g.weights.resize(static_cast<Index>(total));auto cursor=g.col;
        visit([&](Index source,Index target,double w){if(w*scale!=0){const auto k=cursor[source]++;g.row[k]=static_cast<std::uint16_t>(target);g.weights[k]=w*scale;}});
        g.compact();return g;
    }
    double weight(Index k) const {return uniform?common:weights[k];}
    void scatter(Index source,double boost,Vec& accumulator) const {
        if(uniform){const double amount=boost*common;for(Index k=col[source];k<col[source+1];++k)accumulator[row[k]]+=amount;}
        else for(Index k=col[source];k<col[source+1];++k)accumulator[row[k]]+=boost*weights[k];
    }
};
// SIMD lanes cover independent time samples. Each lane keeps the original
// ascending source order; no reduction, FMA, or change to the numerical sum.
struct ExternalBlock {
    static constexpr Index width=16;
    Vec data;
    Index start=0,count=0,neurons;
    explicit ExternalBlock(Index n):data(width*n),neurons(n) {}
    const double* at(Index sample) const {return data.data()+(sample-start)*neurons;}
    void fill(Index sample,Index total,const Sparse& graph,const double* input,Index ld,double scale) {
        start=sample;count=std::min(width,total-sample);
        const __m256d scaled=_mm256_set1_pd(scale);
        for(Index j=0;j<neurons;++j) {
            if(count==width) {
                __m256d a=_mm256_setzero_pd(),b=a,c=a,d=a;
                for(Index k=graph.col[j];k<graph.col[j+1];++k) {
                    const double* x=input+sample+ld*graph.row[k];
                    const __m256d w=_mm256_set1_pd(graph.weight(k));
                    a=_mm256_add_pd(a,_mm256_mul_pd(_mm256_loadu_pd(x),w));
                    b=_mm256_add_pd(b,_mm256_mul_pd(_mm256_loadu_pd(x+4),w));
                    c=_mm256_add_pd(c,_mm256_mul_pd(_mm256_loadu_pd(x+8),w));
                    d=_mm256_add_pd(d,_mm256_mul_pd(_mm256_loadu_pd(x+12),w));
                }
                double values[width];
                _mm256_storeu_pd(values,_mm256_mul_pd(a,scaled));
                _mm256_storeu_pd(values+4,_mm256_mul_pd(b,scaled));
                _mm256_storeu_pd(values+8,_mm256_mul_pd(c,scaled));
                _mm256_storeu_pd(values+12,_mm256_mul_pd(d,scaled));
                for(Index t=0;t<width;++t) data[t*neurons+j]=values[t];
            } else {
                for(Index t=0;t<count;++t) {
                    double sum=0;
                    for(Index k=graph.col[j];k<graph.col[j+1];++k)
                        sum=sum+input[sample+t+ld*graph.row[k]]*graph.weight(k);
                    data[t*neurons+j]=sum*scale;
                }
            }
        }
    }
};
struct History {
    std::vector<Vec> times;
    std::vector<Index> p100, p200;
    explicit History(Index n): times(n), p100(n,0), p200(n,0) {
        for(auto& v: times) v.reserve(128);
    }
    Index count(Index i, double t, double offset, std::vector<Index>& pointers) {
        auto& p=pointers[i]; const Vec& v=times[i];
        // Preserve (spc-t)+offset, including equality; dt>0 means unique times.
        while (p<v.size() && ((v[p]-t)+offset)<0) ++p;
        return v.size()-p;
    }
};
struct Recording {
    Result& result;
    // Time-major staging changes data movement only; all arithmetic is unchanged.
    struct Block {
        static constexpr Index capacity=32;
        double* destination=nullptr;
        Index rows=0,cols=0,start=0,used=0;
        Vec data;
        void init(double* dst,Index nr,Index nc) {destination=dst;rows=nr;cols=nc;data.resize(capacity*nc);}
        double* begin(Index r) {if(!used) start=r;return data.data()+used*cols;}
        void end() {++used;if(used==capacity) flush();}
        void flush() {
            if(!used) return;
            // 16-column tiles keep both source rows and destination lines hot.
            for(Index tile=0;tile<cols;tile+=16) for(Index j=tile;j<std::min(tile+16,cols);++j)
                for(Index r=0;r<used;++r) destination[start+r+rows*j]=data[r*cols+j];
            used=0;
        }
    } ve,vi,he,hi,amb,lgn,l6;
    Index rows;
    Recording(Result& output,Index nr, Index ne, Index ni):result(output),rows(nr) {
        ve.init(add("VE",ne),nr,ne);vi.init(add("VI",ni),nr,ni);
        he.init(add("HE",ne+ni),nr,ne+ni);hi.init(add("HI",ne+ni),nr,ne+ni);
        l6.init(add("L6",ne+ni),nr,ne+ni);amb.init(add("amb",ne),nr,ne);lgn.init(add("lgn",ne),nr,ne);
    }
    void finish() {ve.flush();vi.flush();he.flush();hi.flush();l6.flush();amb.flush();lgn.flush();}
    double* add(const char* name, Index columns) {
        Matrix& a=result.arrays[name];a.rows=rows;a.columns=columns;a.data.resize(rows*columns,0);return a.data.data();
    }
    void spikes(const History& e, const History& i) {
        Index re=1000,ri=1000;
        for(const auto& v:e.times) re=std::max(re,v.size()+1);
        for(const auto& v:i.times) ri=std::max(ri,v.size()+1);
        re=std::max(re,ri);
        Matrix& a=result.arrays["spike"];a.rows=re;a.columns=e.times.size()+i.times.size();a.data.resize(re*a.columns,0);
        double* dst=a.data.data();Index col=0;
        for(const History* h:{&e,&i}) for(const Vec& v:h->times) {
            dst[col*re]=static_cast<double>(v.size());
            for(Index k=0;k<v.size();++k) dst[col*re+k+1]=v[k]/1000;
            ++col;
        }
        for(const Vec& v:e.times)result.eSpikes+=v.size();
        for(const Vec& v:i.times)result.iSpikes+=v.size();
    }
};

} // namespace
Result simulate(const Node& parameters,const Node& initial,const Node& external,const Options& options) {
    constexpr Index ne=13689,ni=2401;
    Result result;
    const Node* c=&parameters;const Node* ex=&external;
    const auto start=std::chrono::steady_clock::now();auto tick=start;
    auto mark=[&](int section){if(options.profile){const auto now=std::chrono::steady_clock::now();result.sections[section]+=std::chrono::duration<double>(now-tick).count();tick=now;}};
    const double dt=value(c,"dt"),duration=(options.durationSeconds>0?options.durationSeconds:value(c,"duration"))*1000;
    require(dt>0 && std::isfinite(dt) && std::isfinite(duration) && duration>=100,"Use positive dt and duration at least 0.1 s.");
    require(std::abs(duration/100-std::round(duration/100))<1e-10,"Duration must be a multiple of 0.1 s.");
    const Index steps=static_cast<Index>(std::floor(duration/dt));
    const Index nr=static_cast<Index>(std::ceil(duration/dt/10))+1;
    result.steps=steps-1;result.durationSeconds=duration/1000;result.dtMilliseconds=dt;
    Sparse ee=Sparse::transpose(c->field("cmee"),value(c,"see"));
    Sparse ei=Sparse::transpose(c->field("cmie"),value(c,"sie"));
    Sparse ie=Sparse::transpose(c->field("cmei"),value(c,"sei"));
    Sparse ii=Sparse::transpose(c->field("cmii"),value(c,"sii"));
    Sparse gb=Sparse::transpose(c->field("cmei"),1,true);
    require(ee.rows==ne && ee.cols==ne && ei.rows==ni && ei.cols==ne && ie.rows==ne && ie.cols==ni && ii.rows==ni && ii.cols==ni,"Expected 13689 E and 2401 I neurons.");
    Sparse lgnE(field(ex,"conn_matrix_e")),lgnI(field(ex,"conn_matrix_i"));
    const Node& exe=ex->field("exe");const Node& exi=ex->field("exi");
    const double* xe=doubles(&exe);const double* xi=doubles(&exi);
    const Index lde=exe.dims[0],ldi=exi.dims[0];
    require(lde>=steps && ldi>=steps && exe.dims[1]==lgnE.rows && exi.dims[1]==lgnI.rows && lgnE.cols==ne && lgnI.cols==ni,"External inputs do not cover the requested simulation.");
    const double selgn=value(c,"selgn"),silgn=value(c,"silgn"),segabab=value(c,"segabab");
    const double meanE=0.0124*0.8,meanI=0.0124*0.75,theta=0.4*0.2;
    const double noiseE=(0.1*std::sqrt(0.8*0.2)*meanE)*std::sqrt(dt),noiseI=(0.1*std::sqrt(0.8*0.2)*meanI)*std::sqrt(dt);
    const double ndepE=value(c,"ndep_e"),ndepI=value(c,"ndep_i");
    const double sbasic=value(c,"sie"),sincre=value(c,"sincre_e"),slob=sbasic-5*sincre,initialSei=value(c,"sei");
    double decay[9];const double taus[]={0.5,2,2,80,0.5,5,3,100,5.05};
    for(Index j=0;j<9;++j)decay[j]=std::exp(-dt/taus[j]);
    const double invamp=1/(2-0.5),ampnmda=(2-0.5)/(80-2),invgab=1/(5-0.5),invgabab=1.0/(100-3),vE=14.0/3,vI=2.0/3;
    auto optionalScale=[&](double overrideVal,const char* fieldName,double defaultVal){
        double s=defaultVal;
        if(overrideVal>0) s=overrideVal;
        else { auto it=c->fields.find(fieldName); if(it!=c->fields.end()) s=it->second.scalar(); }
        require(std::isfinite(s) && s>0,"F1 depression scales must be positive finite numbers.");
        return s;
    };
    const double eeScale=optionalScale(options.f1EeScale,"f1_ee_scale",2.25);
    const double eiScale=optionalScale(options.f1EiScale,"f1_ei_scale",1.5);
    result.f1EeScale=eeScale;result.f1EiScale=eiScale;
    double eDecay[]={1.2,1.2,1.2,1.20,1.19,1.16,1.12,1.07,1,0.96,0.94,0.92};
    double iDecay[]={1.1,1.02,0.98,0.96,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95};
    for(Index j=0;j<12;++j){eDecay[j]=(eDecay[j]-0.92)*eeScale/3+0.92;iDecay[j]=(iDecay[j]-0.95)*eiScale/1.5+0.95;}
    const double eiDecay[]={1.05,1,0.98};
    const Vec seetfm=numeric(c->field("seetfm")),seiList=numeric(c->field("sei_list"));const Index ntfm=seetfm.size(),nSei=seiList.size();
    Vec probE(steps+1),probI(steps+1);
    for(Index j=0;j<=steps;++j){probE[j]=std::pow(value(c,"pdep_e"),static_cast<double>(j)-ndepE+1);probI[j]=std::pow(value(c,"pdep_i"),static_cast<double>(j)-ndepI+1);}
    Broadcast seifac(c,"seifac",ne),sel6(c,"sel6",ne),sil6(c,"sil6",ni);
    Vec refE=vectorField(c,"tau_re",ne),refI=vectorField(c,"tau_ri",ni);
    Rng rng(options.seed);
    Vec ve(ne,0),vi(ni,0),hee(4*ne,0),hie(4*ni,0),hei(4*ne,0),hii(2*ni,0);
    const bool hasInit=!options.emptyInitialization && !initial.fields.empty();
    if(hasInit){
        ve=vectorField(&initial,"ve",ne);vi=vectorField(&initial,"vi",ni);
        auto fillState=[&](Vec& state,Index n,Index components,double scale){for(Index j=0;j<n;++j)for(Index k=0;k<components;++k)state[k*n+j]=scale*rng.uniform();};
        fillState(hee,ne,4,0.03);fillState(hie,ni,4,0.03);fillState(hei,ne,4,0.05);fillState(hii,ni,2,0.05);
    }
    Vec heSum(ne),hiESum(ni),hiSum(ne),hiISum(ni),ambE(ne),ambI(ni);
    for(Index j=0;j<ne;++j){heSum[j]=0.8*(hee[ne+j]-hee[j])+0.2*(hee[3*ne+j]-hee[2*ne+j]);hiSum[j]=hei[ne+j]-hei[j]+hei[3*ne+j]-hei[2*ne+j];ambE[j]=rng.uniform()*2*meanE;}
    for(Index j=0;j<ni;++j){hiESum[j]=0.8*(hie[ni+j]-hie[j])+0.2*(hie[3*ni+j]-hie[2*ni+j]);hiISum[j]=hii[ni+j]-hii[j];ambI[j]=rng.uniform()*2*meanI;}
    Vec refCountE(ne,0),refCountI(ni,0),threshold(ne,0),seeCheck(ne,1),seiCheck(ni,1),sie(ne,sbasic),sei(ni,initialSei);
    std::vector<Index> seeRec(ne,1),seiRec(ni,1),spikingE,spikingI;spikingE.reserve(ne);spikingI.reserve(ni);
    History histE(ne),histI(ni);
    Vec externalE(ne),externalI(ni),ce(ne,0),ci(ni,0),gabaSum(ne,0);
    std::vector<std::uint16_t> gabaHistory(300*ne,0);ExternalBlock inputE(ne),inputI(ni);
    Vec l6DataE(10*ne),l6DataI(10*ni),l6Raw(ne+ni),l6Fac(ne),fre(ne),fri(ni),uniformL6(10*(ne+ni));
    std::vector<Index> mapEI(ni);const Vec map=numeric(c->field("map_ei"));require(map.size()==ni,"Invalid E-to-I mapping size.");
    for(Index j=0;j<ni;++j){require(map[j]>=1 && map[j]<=ne && map[j]==std::floor(map[j]),"Invalid E-to-I map.");mapEI[j]=static_cast<Index>(map[j])-1;}
    L6Config l6Config;l6Config.mask=numeric(c->field("mask"));l6Config.transform=numeric(c->field("L6tfm"));l6Config.coefficient=numeric(c->field("l6cof"));
    l6Config.facilitationMaximum=numeric(c->field("l6facmax"));l6Config.mapEI=mapEI;l6Config.inhibitoryFacilitation=value(c,"l6facier");
    L6 l6(std::move(l6Config));const double facIER=l6.inhibitoryFacilitation();
    std::vector<unsigned char> shortHistory(70*ne,0);Vec shortSum(ne,0),l6Sum(ne,0),depE(3*ne,0),depI(3*ni,0);
    std::vector<std::uint8_t> l6Window(100*ne,0);Index shortIndex=0,gabaIndex=0,l6Count=0,depCount=0;
    Vec normals(ne+ni),releases(ne+ni);Recording rec(result,nr,ne,ni);
    std::uint64_t eventEdges=0,randomReleases=0;
    mark(0);
    for(Index step=2;step<=steps;++step) {
        const double t=dt*static_cast<double>(step);
        if(step-1>=inputE.start+inputE.count) {
            inputE.fill(step-1,steps,lgnE,xe,lde,selgn);
            inputI.fill(step-1,steps,lgnI,xi,ldi,silgn);
        }
        mark(11);
        rng.fillNormal(normals.data(),normals.size());
        const double* normalE=normals.data();const double* normalI=normalE+ne;
        mark(1);
        v1fast::advance(ne,ve.data(),refCountE.data(),ambE.data(),externalE.data(),inputE.at(step-1),normalE,
            heSum.data(),hiSum.data(),threshold.data(),hee.data(),hei.data(),4,dt,theta,meanE,noiseE,20,vE,vI,decay);
        v1fast::advance(ni,vi.data(),refCountI.data(),ambI.data(),externalI.data(),inputI.at(step-1),normalI,
            hiESum.data(),hiISum.data(),nullptr,hie.data(),hii.data(),2,dt,theta,meanI,noiseI,15,vE,vI,decay);
        mark(2);
        spikingE.clear();spikingI.clear(); Index needRandom=0;
        v1fast::detect(ne,ve.data(),refCountE.data(),refE.data(),threshold.data(),spikingE);
        v1fast::detect(ni,vi.data(),refCountI.data(),refI.data(),nullptr,spikingI);
        for(Index j:spikingE) if(!(histE.times[j].size()<ndepE)) ++needRandom;
        for(Index j:spikingI) if(!(histI.times[j].size()<ndepI)) ++needRandom;
        // rand calls can be batched: eligibility is known before processing,
        // and distinct neurons cannot affect each other's history this step.
        const double* random=releases.data();Index randomIndex=0;
        if(needRandom)rng.fillUniform(releases.data(),needRandom);
        randomReleases+=needRandom;
        for(Index j=0;j<ne;++j) {
            const Index si=shortIndex*ne+j,gi=gabaIndex*ne+j;
            shortSum[j]-=shortHistory[si]; shortHistory[si]=0;
            gabaSum[j]-=gabaHistory[gi]; gabaHistory[gi]=0;
        }
        mark(3);
        std::fill(ce.begin(),ce.end(),0); std::fill(ci.begin(),ci.end(),0);
        bool anyE=false;
        for(Index ind:spikingE) {
            const Index n=histE.times[ind].size();
            const Index prev200=histE.count(ind,t,199.99,histE.p200);
            const Index prev=histE.count(ind,t,99.99,histE.p100);
            const bool accepted=n<ndepE || random[randomIndex++]<probE[prev];
            if(!accepted) continue;
            const double tpre=n?histE.times[ind].back():0;
            double boostE=0,boostI=0;
            auto tfm=[&]() { require(prev200<ntfm,"seetfm index exceeds reference lookup table."); return seetfm[prev200]; };
            if(seeCheck[ind]>0.92) {
                if(t-tpre<50) seeRec[ind]=std::min<Index>(seeRec[ind]+1,12);
                else if(t-tpre>133) seeRec[ind]=seeRec[ind]>1?seeRec[ind]-1:1;
                boostE=eDecay[seeRec[ind]-1];boostI=iDecay[seeRec[ind]-1];
            } else if(seeCheck[ind]<0.92) {
                boostE=tfm();boostI=tfm()*sie[ind]/0.0176;
            } else {
                if(tfm()<0.92) { boostE=tfm();boostI=tfm()*sie[ind]/0.0176; }
                else if(t-tpre>133) {
                    seeRec[ind]=seeRec[ind]>1?seeRec[ind]-1:1;
                    boostE=eDecay[seeRec[ind]-1];boostI=iDecay[seeRec[ind]-1];
                } else { boostE=0.92;boostI=0.92; }
            }
            seeCheck[ind]=boostE;
            ee.scatter(ind,boostE,ce);ei.scatter(ind,boostI,ci);
            eventEdges+=ee.col[ind+1]-ee.col[ind]+ei.col[ind+1]-ei.col[ind];
            anyE=anyE || boostE!=0 || boostI!=0;
            threshold[ind]=0.5;histE.times[ind].push_back(t);
            shortHistory[shortIndex*ne+ind]=1;shortSum[ind]+=1;
        }
        if(anyE) {
            v1fast::addExcitation(ne,hee.data(),ce.data(),invamp,ampnmda);
            v1fast::addExcitation(ni,hie.data(),ci.data(),invamp,ampnmda);
        }
        std::fill(ce.begin(),ce.end(),0); std::fill(ci.begin(),ci.end(),0); bool anyI=false;
        for(Index ind:spikingI) {
            const Index n=histI.times[ind].size();
            const Index prev=histI.count(ind,t,99.99,histI.p100);
            const bool accepted=n<ndepI || random[randomIndex++]<probI[prev];
            if(!accepted) continue;
            const double tpre=n?histI.times[ind].back():0;double boost=0;
            if(seiCheck[ind]>0.98) {
                if(t-tpre<25) seiRec[ind]=std::min<Index>(seiRec[ind]+1,3);
                else if(t-tpre>50) seiRec[ind]=seiRec[ind]>1?seiRec[ind]-1:1;
                boost=eiDecay[seiRec[ind]-1];
            } else if(seiCheck[ind]<0.98) { boost=sei[ind]/0.051; }
            else {
                if(sei[ind]/0.051<0.98) boost=sei[ind]/0.051;
                else if(t-tpre>50) { seiRec[ind]=seiRec[ind]>1?seiRec[ind]-1:1;boost=eiDecay[seiRec[ind]-1]; }
                else boost=0.98;
            }
            seiCheck[ind]=boost;
            for(Index k=ie.col[ind];k<ie.col[ind+1];++k) {
                const Index target=ie.row[k];const double w=ie.weight(k)*boost-seifac[target];
                if(w>0) ce[target]=ce[target]+w;
            }
            ii.scatter(ind,1,ci);
            eventEdges+=ie.col[ind+1]-ie.col[ind]+ii.col[ind+1]-ii.col[ind];
            for(Index k=gb.col[ind];k<gb.col[ind+1];++k) {
                const Index target=gb.row[k];
                ++gabaHistory[gabaIndex*ne+target];gabaSum[target]+=1;
            }
            histI.times[ind].push_back(t);anyI=true;
        }
        require(randomIndex==needRandom,"Internal random stream accounting mismatch.");
        if(anyI) {
            v1fast::addInhibition(ne,hei.data(),ce.data(),invgab);
            v1fast::addInhibition(ni,hii.data(),ci.data(),invgab);
        }
        mark(4);
        shortIndex=(shortIndex+1)%70;gabaIndex=(gabaIndex+1)%300;
        bool gabaTriggered=false;
        for(double total:gabaSum) if(total>265) {gabaTriggered=true;break;}
        if(gabaTriggered) for(Index j=0;j<ne;++j) {
            const bool triggered=gabaSum[j]>265;
            const double cb=static_cast<double>(triggered)*segabab*invgabab;
            hei[2*ne+j]+=cb;hei[3*ne+j]+=cb;
            if(triggered) {
                for(Index k=0;k<300;++k) gabaHistory[k*ne+j]=0;
                gabaSum[j]=0;
            }
        }
        mark(5);
        const Index l6offset=(step-2)%10;
        if(l6offset==0) {
            const Index slot=l6Count%100;
            if(l6Count>=100) for(Index j=0;j<ne;++j) l6Sum[j]-=l6Window[slot*ne+j];
            l6.probabilities(shortSum.data(),fre.data(),fri.data());
            rng.fillUniform(uniformL6.data(),uniformL6.size());
            const double* ep=fre.data();const double* ip=fri.data();const double* rawp=uniformL6.data();
                const bool enabled=l6Count+1>=100;
                for(Index j=0;j<ne;++j) {
                    double samples[10];unsigned count=0;
                    for(Index sampleIndex=0;sampleIndex<10;++sampleIndex) {samples[sampleIndex]=rawp[j*10+sampleIndex]<ep[j];count+=static_cast<unsigned>(samples[sampleIndex]);}
                    l6Raw[j]=count;
                    const Index total=static_cast<Index>(l6Sum[j])+count;
                    if(total>1000) throw std::runtime_error("L6 facilitation count exceeds exact lookup bounds.");
                    l6Fac[j]=enabled?l6.facilitation(j,static_cast<unsigned>(total)):1;
                    for(Index sampleIndex=0;sampleIndex<10;++sampleIndex) l6DataE[sampleIndex*ne+j]=samples[sampleIndex]*l6Fac[j];
                }
                for(Index j=0;j<ni;++j) {
                    const Index mapped=mapEI[j];
                    unsigned count=0;
                    for(Index sampleIndex=0;sampleIndex<10;++sampleIndex) {
                        const double sample=rawp[(ne+j)*10+sampleIndex]<ip[j];count+=static_cast<unsigned>(sample);
                        l6DataI[sampleIndex*ni+j]=enabled?(sample*l6Fac[mapped])*facIER:sample;
                    }
                    l6Raw[ne+j]=count;
                }
            rawp=l6Raw.data();
            std::copy(rawp,rawp+ne+ni,rec.l6.begin(l6Count));rec.l6.end();
            for(Index j=0;j<ne;++j) {l6Window[slot*ne+j]=static_cast<std::uint8_t>(rawp[j]);l6Sum[j]+=rawp[j];}
            ++l6Count;
        }
        mark(6);
        v1fast::finishL6(ne,hee.data(),hei.data(),4,heSum.data(),hiSum.data(),l6DataE.data()+ne*l6offset,
            sel6.p,sel6.n,invamp,ampnmda);
        v1fast::finishL6(ni,hie.data(),hii.data(),2,hiESum.data(),hiISum.data(),l6DataI.data()+ni*l6offset,
            sil6.p,sil6.n,invamp,ampnmda);
        mark(7);
        if(step%1000==0) {
            const Index slot=depCount%3,old=(depCount+1)%3;
            for(Index j=0;j<ne;++j) {
                const double count=static_cast<double>(histE.times[j].size());
                depE[slot*ne+j]=count;
                if(depCount>=2) {
                    const double diff=count-depE[old*ne+j];
                    const double delta=sbasic-(diff-6)*sincre-sie[j];
                    const double sign=(delta>0)-(delta<0);
                    sie[j]=sie[j]+sign*sincre;
                    sie[j]=std::max(std::min(sie[j],sbasic),slob);
                }
            }
            for(Index j=0;j<ni;++j) {
                const double count=static_cast<double>(histI.times[j].size());depI[slot*ni+j]=count;
                if(depCount>=2) {
                    const Index diff=static_cast<Index>(count-depI[old*ni+j]);
                    require(diff<nSei,"sei_list index exceeds reference lookup table.");sei[j]=seiList[diff];
                }
            }
            ++depCount;
        }
        mark(8);
        if(step%10==0) {
            const Index r=step/10-1;
            double* outVe=rec.ve.begin(r);double* outVi=rec.vi.begin(r);
            double* outHe=rec.he.begin(r);double* outHi=rec.hi.begin(r);
            double* outAmb=rec.amb.begin(r);double* outLgn=rec.lgn.begin(r);
            for(Index j=0;j<ne;++j) {
                outVe[j]=ve[j];outHe[j]=heSum[j];outHi[j]=hiSum[j];
                outAmb[j]=ambE[j];outLgn[j]=externalE[j]-ambE[j];
            }
            for(Index j=0;j<ni;++j) {
                outVi[j]=vi[j];outHe[ne+j]=hiESum[j];outHi[ne+j]=hiISum[j];
            }
            rec.ve.end();rec.vi.end();rec.he.end();rec.hi.end();rec.amb.end();rec.lgn.end();
        }
        mark(9);
    }
    rec.finish();
    rec.spikes(histE,histI);
    mark(10);
    result.eventEdges=eventEdges;result.releaseDraws=randomReleases;
    result.simulationSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    return result;
}
std::vector<mat5::FieldView> outputFields(const Result& result) {
    std::vector<mat5::FieldView> fields;
    for(const char* name:{"VE","VI","HE","HI","L6","amb","lgn","spike"}) {
        const Matrix& a=result.arrays.at(name);fields.push_back({name,{a.rows,a.columns},a.data.data(),a.data.size()});
    }
    return fields;
}
std::string summaryJson(const Result& result,const Options& options,double loadSeconds,double saveSeconds) {
    std::ostringstream s;s<<std::setprecision(17);
    s<<"{\n  \"seed\": "<<options.seed<<",\n  \"rng\": \""<<Rng::algorithm<<"\",\n  \"f1_ee_scale\": "<<result.f1EeScale
     <<",\n  \"f1_ei_scale\": "<<result.f1EiScale<<",\n  \"duration_seconds\": "<<result.durationSeconds
     <<",\n  \"dt_ms\": "<<result.dtMilliseconds<<",\n  \"steps\": "<<result.steps
     <<",\n  \"load_seconds\": "<<loadSeconds<<",\n  \"simulation_seconds\": "<<result.simulationSeconds<<",\n  \"save_seconds\": "<<saveSeconds
     <<",\n  \"e_spikes\": "<<result.eSpikes<<",\n  \"i_spikes\": "<<result.iSpikes
     <<",\n  \"e_rate_hz\": "<<result.eSpikes/(13689*result.durationSeconds)<<",\n  \"i_rate_hz\": "<<result.iSpikes/(2401*result.durationSeconds)
     <<",\n  \"event_edges\": "<<result.eventEdges<<",\n  \"release_draws\": "<<result.releaseDraws;
    const auto& spikes=result.arrays.at("spike");
    auto cellCV=[&](Index begin,Index end) {
        double sum=0,sum2=0;for(Index j=begin;j<end;++j){const double count=spikes.data[j*spikes.rows];sum+=count;sum2+=count*count;}
        const double count=static_cast<double>(end-begin),mean=sum/count;
        return mean>0?std::sqrt(std::max(0.0,sum2/count-mean*mean))/mean:0;
    };
    s<<",\n  \"e_cellrate_cv\": "<<cellCV(0,13689)<<",\n  \"i_cellrate_cv\": "<<cellCV(13689,16090);
    s<<",\n  \"fields\": {";bool first=true;
    for(const auto& entry:result.arrays){
        const auto& a=entry.second;bool finite=true;double minimum=std::numeric_limits<double>::infinity(),maximum=-minimum;long double sum=0,sum2=0;
        for(double v:a.data){finite=finite&&std::isfinite(v);minimum=std::min(minimum,v);maximum=std::max(maximum,v);sum+=v;sum2+=static_cast<long double>(v)*v;}
        if(!finite)throw std::runtime_error("Simulation produced non-finite output in "+entry.first);
        if(!first)s<<",";first=false;
        s<<"\n    \""<<entry.first<<"\": {\"rows\": "<<a.rows<<", \"columns\": "<<a.columns<<", \"finite\": true, \"mean\": "<<static_cast<double>(sum/a.data.size())
         <<", \"mean_square\": "<<static_cast<double>(sum2/a.data.size())<<", \"min\": "<<minimum<<", \"max\": "<<maximum<<"}";
    }
    s<<"\n  },\n  \"section_seconds\": {";
    const char* names[]={"setup","normal_rng","membrane_lgn_decay","detection_release_rng_expire","events","gabab_reset","l6","l6_integration_sums","homeostasis","recording","finish","lgn_projection"};
    for(Index k=0;k<12;++k){if(k)s<<",";s<<"\n    \""<<names[k]<<"\": "<<result.sections[k];}
    s<<"\n  }\n}\n";return s.str();
}
} // namespace v1standalone

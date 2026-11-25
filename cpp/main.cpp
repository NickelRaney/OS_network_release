#include "v1_simulator.hpp"
#include "v1_rng.hpp"
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _MSC_VER
#include <intrin.h>
#endif

namespace {
using Clock=std::chrono::steady_clock;
double elapsed(Clock::time_point start){return std::chrono::duration<double>(Clock::now()-start).count();}
void usage(){
    std::cout<<"V1_4 standalone simulator (C++17, double, AVX2)\n"
    <<"Usage: v1_4_standalone --input parameters.mat [options]\n"
    <<"  --output result.mat   Save the same eight result fields as a MAT struct res\n"
    <<"  --no-save             Compute all outputs in memory without writing a MAT\n"
    <<"  --summary report.json Save timings and result statistics\n"
    <<"  --seed INTEGER        Native RNG seed (default 20261005)\n"
    <<"  --duration SECONDS    Override saved duration (multiples of 0.1 s)\n"
    <<"  --f1-ee SCALE         Override F1 E->E scale (default: param field or 2.25)\n"
    <<"  --f1-ei SCALE         Override F1 E->I scale (default: param field or 1.5)\n"
    <<"  --external VARNAME    External input struct name in MAT file (default ex_vert_long)\n"
    <<"  --repeat N            Repeat from the same seed; save last result only\n"
    <<"  --empty-init          Start with zero voltages/synaptic states\n"
    <<"  --profile             Include component timings\n"
    <<"  --force               Allow replacement of explicitly selected output files\n";
}
bool avx2Available(){
#ifdef _MSC_VER
    int info[4];__cpuid(info,1);
    if(!(info[2]&(1<<27)) || !(info[2]&(1<<28)))return false;
    if((_xgetbv(0)&6)!=6)return false;
    __cpuidex(info,7,0);return (info[1]&(1<<5))!=0;
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
    return __builtin_cpu_supports("avx2");
#else
    return false;
#endif
}
std::uint64_t integer(const std::string& text){
    if(text.empty() || text[0]=='-')throw std::runtime_error("Expected a nonnegative integer.");
    std::size_t consumed=0;auto value=std::stoull(text,&consumed);
    if(consumed!=text.size())throw std::runtime_error("Invalid integer: "+text);
    return value;
}
bool samePath(const std::string& a,const std::string& b){
    std::error_code error;
    if(std::filesystem::equivalent(a,b,error) && !error)return true;
    const auto first=std::filesystem::weakly_canonical(std::filesystem::absolute(a));
    const auto second=std::filesystem::weakly_canonical(std::filesystem::absolute(b));
#ifdef _WIN32
    return _wcsicmp(first.c_str(),second.c_str())==0;
#else
    return first==second;
#endif
}
}
int main(int argc,char** argv){
    try {
        const auto programStart=Clock::now();
        std::string input,output="v1_4_result.mat",summary,externalVar="ex_vert_long";
        bool save=true,force=false;std::size_t repeats=1;
        v1standalone::Options options;
        for(int k=1;k<argc;++k){
            const std::string arg=argv[k];
            auto next=[&](){if(k+1>=argc)throw std::runtime_error("Missing value for "+arg);return std::string(argv[++k]);};
            if(arg=="--input")input=next();else if(arg=="--output")output=next();
            else if(arg=="--summary")summary=next();else if(arg=="--seed")options.seed=integer(next());
            else if(arg=="--duration"){const auto text=next();std::size_t used;options.durationSeconds=std::stod(text,&used);if(used!=text.size() || !(options.durationSeconds>0))throw std::runtime_error("Invalid duration.");}
            else if(arg=="--f1-ee"){const auto text=next();std::size_t used;options.f1EeScale=std::stod(text,&used);if(used!=text.size() || !(options.f1EeScale>0))throw std::runtime_error("Invalid --f1-ee scale.");}
            else if(arg=="--f1-ei"){const auto text=next();std::size_t used;options.f1EiScale=std::stod(text,&used);if(used!=text.size() || !(options.f1EiScale>0))throw std::runtime_error("Invalid --f1-ei scale.");}
            else if(arg=="--external")externalVar=next();
            else if(arg=="--repeat")repeats=static_cast<std::size_t>(integer(next()));
            else if(arg=="--no-save")save=false;else if(arg=="--empty-init")options.emptyInitialization=true;
            else if(arg=="--profile")options.profile=true;else if(arg=="--force")force=true;
            else if(arg=="--help" || arg=="-h"){usage();return 0;}
            else throw std::runtime_error("Unknown option: "+arg);
        }
        if(input.empty()){usage();return 2;}
        if(repeats<1 || repeats>1000)throw std::runtime_error("Repeat count must be 1..1000.");
        if(!avx2Available())throw std::runtime_error("This build requires an AVX2-capable x86-64 CPU and OS.");
        if(save && !summary.empty() && samePath(output,summary))throw std::runtime_error("MAT output and JSON summary must use different paths.");
        for(const auto& path:std::vector<std::string>{save?output:"",summary}) {
            if(path.empty())continue;
            if(std::filesystem::exists(path) && !force)throw std::runtime_error("Output already exists; use another name or --force: "+path);
            if(samePath(path,input))throw std::runtime_error("Output must not replace the input parameter file.");
        }
        std::cout<<std::fixed<<std::setprecision(6)<<"Loading "<<input<<"\n";
        const auto loadStart=Clock::now();
        auto variables=mat5::read(input,{"param_v1","init_v1",externalVar});
        const double loadSeconds=elapsed(loadStart);
        const auto& param=variables.at("param_v1");const auto& initial=variables.at("init_v1");const auto& external=variables.at(externalVar);
        std::cout<<"Loaded in "<<loadSeconds<<" s; seed "<<options.seed<<"; "<<v1standalone::Rng::algorithm<<"\n";
        std::vector<double> timings;double saveSeconds=0,summarySeconds=0;std::string report;
        for(std::size_t trial=0;trial<repeats;++trial){
            auto result=v1standalone::simulate(param,initial,external,options);
            timings.push_back(result.simulationSeconds);
            std::cout<<"Trial "<<(trial+1)<<": simulation "<<result.simulationSeconds<<" s; E spikes "<<result.eSpikes<<"; I spikes "<<result.iSpikes<<"\n";
            if(trial+1==repeats){
                if(save){const auto timer=Clock::now();mat5::writeStruct(output,"res",v1standalone::outputFields(result));saveSeconds=elapsed(timer);std::cout<<"Saved "<<output<<" in "<<saveSeconds<<" s\n";}
                const auto timer=Clock::now();report=v1standalone::summaryJson(result,options,loadSeconds,saveSeconds);summarySeconds=elapsed(timer);
            }
        }
        auto sorted=timings;std::sort(sorted.begin(),sorted.end());const auto n=sorted.size();const double median=n%2?sorted[n/2]:(sorted[n/2-1]+sorted[n/2])/2;
        const double totalSeconds=elapsed(programStart);
        std::ostringstream extra;extra<<std::setprecision(17)<<"\n  \"summary_seconds\": "<<summarySeconds<<",\n  \"total_seconds\": "<<totalSeconds
            <<",\n  \"median_simulation_seconds\": "<<median<<",\n  \"trial_seconds\": [";
        for(std::size_t j=0;j<timings.size();++j){if(j)extra<<",";extra<<timings[j];}extra<<"],";
        report.insert(1,extra.str());
        if(!summary.empty()){std::ofstream stream(summary);if(!stream)throw std::runtime_error("Cannot open summary: "+summary);stream<<report;if(!stream)throw std::runtime_error("Summary write failed.");}
        std::cout<<"Median simulation "<<median<<" s; result validation/statistics "<<summarySeconds<<" s; total "<<totalSeconds<<" s\n";
        return 0;
    }catch(const std::exception& error){std::cerr<<"Error: "<<error.what()<<"\n";return 1;}
}

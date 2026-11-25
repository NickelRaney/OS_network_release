#pragma once
#include "mat5.hpp"
#include <array>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace v1standalone {
struct Options {
    std::uint64_t seed=20261005;
    double durationSeconds=0; // zero uses the saved parameter
    double f1EeScale=0;       // zero uses param.f1_ee_scale if present, else 2.25
    double f1EiScale=0;       // zero uses param.f1_ei_scale if present, else 1.5
    bool emptyInitialization=false;
    bool profile=false;
};
struct Matrix {
    std::size_t rows=0,columns=0;
    std::vector<double> data;
};
struct Result {
    std::map<std::string,Matrix> arrays;
    double simulationSeconds=0,durationSeconds=0,dtMilliseconds=0;
    double f1EeScale=2.25,f1EiScale=1.5;
    std::array<double,12> sections{};
    std::uint64_t eSpikes=0,iSpikes=0,eventEdges=0,releaseDraws=0;
    std::size_t steps=0;
};
Result simulate(const mat5::Node& parameters,const mat5::Node& initial,
                const mat5::Node& external,const Options& options);
std::vector<mat5::FieldView> outputFields(const Result& result);
std::string summaryJson(const Result& result,const Options& options,
                        double loadSeconds,double saveSeconds);
} // namespace v1standalone

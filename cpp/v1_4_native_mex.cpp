// MATLAB input/output adapter for the independent C++ engine. No MATLAB
// callback, random-number API, MAT file, or simulation formula is used here.
#include "mex.h"
#include "v1_simulator.hpp"
#include "v1_rng.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#ifdef _MSC_VER
#include <intrin.h>
#endif

namespace {
using Clock = std::chrono::steady_clock;
using Node = mat5::Node;
constexpr std::size_t ne = 13689, ni = 2401;
const char* const parameterFields[] = {
    "dt", "duration", "cmee", "cmie", "cmei", "cmii", "see", "sie", "sei", "sii",
    "selgn", "silgn", "segabab", "ndep_e", "ndep_i", "sincre_e", "seetfm", "sei_list",
    "pdep_e", "pdep_i", "seifac", "sel6", "sil6", "tau_re", "tau_ri", "map_ei",
    "mask", "L6tfm", "l6cof", "l6facmax", "l6facier"
};
const char* resultFields[] = {"VE", "VI", "HE", "HI", "L6", "amb", "lgn", "spike"};
const char* sectionFields[] = {
    "setup", "normal_rng", "membrane_lgn_decay", "detection_release_rng_expire",
    "events", "gabab_reset", "l6", "l6_integration_sums", "homeostasis",
    "recording", "finish", "lgn_projection"
};

void require(bool condition, const std::string& message) {
    if (!condition) throw std::invalid_argument(message);
}
double elapsed(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}
struct MxDestroy { void operator()(mxArray* value) const noexcept { if (value) mxDestroyArray(value); } };
using Mx = std::unique_ptr<mxArray, MxDestroy>;
Mx owned(mxArray* value) {
    if (!value) throw std::bad_alloc();
    return Mx(value);
}
void set(mxArray* object, const char* name, Mx value) {
    mxSetField(object, 0, name, value.release());
}
void scalar(mxArray* object, const char* name, double value) {
    set(object, name, owned(mxCreateDoubleScalar(value)));
}
void integer(mxArray* object, const char* name, std::uint64_t value) {
    auto out = owned(mxCreateUninitNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL));
    *static_cast<std::uint64_t*>(mxGetData(out.get())) = value;
    set(object, name, std::move(out));
}
bool avx2Available() {
#ifdef _MSC_VER
    int info[4];
    __cpuid(info, 0);
    if (info[0] < 7) return false;
    __cpuid(info, 1);
    if (!(info[2] & (1 << 27)) || !(info[2] & (1 << 28))) return false;
    if ((_xgetbv(0) & 6) != 6) return false;
    __cpuidex(info, 7, 0);
    return (info[1] & (1 << 5)) != 0;
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
    return __builtin_cpu_supports("avx2");
#else
    return false;
#endif
}
void scalarStruct(const mxArray* value, const std::string& name) {
    require(value && mxIsStruct(value) && mxGetNumberOfElements(value) == 1,
            name + " must be a scalar struct.");
}
const mxArray* field(const mxArray* object, const char* name, const char* owner) {
    const mxArray* value = mxGetField(object, 0, name);
    require(value != nullptr, std::string(owner) + "." + name + " is required.");
    return value;
}
void realMatrix(const mxArray* value, const std::string& name) {
    require((mxIsNumeric(value) || mxIsLogical(value)) && !mxIsComplex(value) &&
                mxGetNumberOfDimensions(value) == 2,
            name + " must be a real numeric or logical matrix.");
}
void shape(const mxArray* value, std::size_t rows, std::size_t columns, const std::string& name) {
    realMatrix(value, name);
    require(mxGetM(value) == rows && mxGetN(value) == columns,
            name + " has the wrong matrix size.");
}
template<class T> void copyNumeric(const mxArray* value, std::vector<double>& out) {
    const T* source = static_cast<const T*>(mxGetData(value));
    for (std::size_t i = 0; i < out.size(); ++i) out[i] = static_cast<double>(source[i]);
}
void numericValues(const mxArray* value, std::vector<double>& out, const std::string& name) {
    switch (mxGetClassID(value)) {
    case mxDOUBLE_CLASS:
        if (!out.empty()) std::memcpy(out.data(), mxGetData(value), out.size() * sizeof(double));
        break;
    case mxSINGLE_CLASS: copyNumeric<float>(value, out); break;
    case mxINT8_CLASS: copyNumeric<std::int8_t>(value, out); break;
    case mxUINT8_CLASS: copyNumeric<std::uint8_t>(value, out); break;
    case mxINT16_CLASS: copyNumeric<std::int16_t>(value, out); break;
    case mxUINT16_CLASS: copyNumeric<std::uint16_t>(value, out); break;
    case mxINT32_CLASS: copyNumeric<std::int32_t>(value, out); break;
    case mxUINT32_CLASS: copyNumeric<std::uint32_t>(value, out); break;
    case mxINT64_CLASS: copyNumeric<std::int64_t>(value, out); break;
    case mxUINT64_CLASS: copyNumeric<std::uint64_t>(value, out); break;
    case mxLOGICAL_CLASS: copyNumeric<mxLogical>(value, out); break;
    default: throw std::invalid_argument(name + " uses an unsupported numeric class.");
    }
    for (double v : out)
        if (!std::isfinite(v)) throw std::invalid_argument(name + " must contain only finite values.");
}
Node numericNode(const mxArray* value, const std::string& name,
                 bool allowSparse = false, bool retainLogical = true) {
    realMatrix(value, name);
    Node node;
    node.dims = {mxGetM(value), mxGetN(value)};
    if (mxIsSparse(value)) {
        require(allowSparse, name + " must be a full matrix.");
        node.kind = Node::Sparse;
        const mwIndex* jc = mxGetJc(value);
        const mwIndex* ir = mxGetIr(value);
        const std::size_t count = jc[node.dims[1]];
        require(count <= UINT32_MAX && count <= mxGetNzmax(value), name + " has too many sparse entries.");
        node.jc.resize(node.dims[1] + 1);
        node.ir.resize(count);
        require(jc[0] == 0, name + " has invalid sparse column pointers.");
        for (std::size_t j = 0; j < node.jc.size(); ++j) {
            if (jc[j] > count || (j > 0 && jc[j - 1] > jc[j]))
                throw std::invalid_argument(name + " has invalid sparse column pointers.");
            node.jc[j] = static_cast<std::uint32_t>(jc[j]);
        }
        for (std::size_t j = 0; j < count; ++j) {
            if (ir[j] >= node.dims[0] || ir[j] > UINT32_MAX)
                throw std::invalid_argument(name + " has an invalid sparse row index.");
            node.ir[j] = static_cast<std::uint32_t>(ir[j]);
        }
        // The engine's CSC graph representation uses real values even for
        // logical sparse matrices; dense logical connectivity stays compact.
        node.values.resize(count);
        numericValues(value, node.values, name);
    } else if (mxIsLogical(value) && retainLogical) {
        node.logical = true;
        const std::size_t count = mxGetNumberOfElements(value);
        node.bytes.resize(count);
        static_assert(sizeof(mxLogical) == sizeof(std::uint8_t), "Unexpected MATLAB logical size");
        if (count) std::memcpy(node.bytes.data(), mxGetData(value), count);
    } else {
        node.values.resize(mxGetNumberOfElements(value));
        numericValues(value, node.values, name);
    }
    return node;
}
Node structure() {
    Node node;
    node.kind = Node::Struct;
    node.dims = {1, 1};
    return node;
}
v1standalone::Options readOptions(const mxArray* value) {
    v1standalone::Options options;
    if (!value) return options;
    scalarStruct(value, "options");
    for (int i = 0; i < mxGetNumberOfFields(value); ++i) {
        const std::string name = mxGetFieldNameByNumber(value, i);
        require(name == "seed" || name == "profile", "Unknown option: " + name);
    }
    if (const mxArray* seed = mxGetField(value, 0, "seed")) {
        require(!mxIsComplex(seed) && !mxIsSparse(seed) && mxGetNumberOfElements(seed) == 1,
                "options.seed must be a real uint64 or double scalar.");
        if (mxIsUint64(seed)) {
            options.seed = *static_cast<const std::uint64_t*>(mxGetData(seed));
        } else {
            require(mxIsDouble(seed), "options.seed must be uint64 or double.");
            const double number = mxGetScalar(seed);
            require(std::isfinite(number) && number >= 0 && number <= 9007199254740992.0 && number == std::floor(number),
                    "A double seed must be an integer from 0 through flintmax; use uint64 for larger seeds.");
            options.seed = static_cast<std::uint64_t>(number);
        }
    }
    if (const mxArray* profile = mxGetField(value, 0, "profile")) {
        require(mxIsLogicalScalar(profile), "options.profile must be a logical scalar.");
        options.profile = mxIsLogicalScalarTrue(profile);
    }
    return options;
}
struct Inputs { Node parameters, initial, external; };
Inputs readInputs(const mxArray* parameters, const mxArray* initial, const mxArray* external) {
    scalarStruct(parameters, "param");
    scalarStruct(external, "ex");
    shape(field(parameters, "cmee", "param"), ne, ne, "param.cmee");
    shape(field(parameters, "cmie", "param"), ne, ni, "param.cmie");
    shape(field(parameters, "cmei", "param"), ni, ne, "param.cmei");
    shape(field(parameters, "cmii", "param"), ni, ni, "param.cmii");
    shape(field(parameters, "mask", "param"), 117, 117, "param.mask");
    Inputs input;
    input.parameters = structure();
    input.external = structure();
    for (const char* name : parameterFields) {
        const bool graph = std::strcmp(name, "cmee") == 0 || std::strcmp(name, "cmie") == 0 ||
                           std::strcmp(name, "cmei") == 0 || std::strcmp(name, "cmii") == 0;
        input.parameters.fields.emplace(name,
            numericNode(field(parameters, name, "param"), std::string("param.") + name, graph));
    }
    for (const char* optName : {"f1_ee_scale", "f1_ei_scale"}) {
        if (const mxArray* optVal = mxGetField(parameters, 0, optName)) {
            input.parameters.fields.emplace(optName,
                numericNode(optVal, std::string("param.") + optName, false));
        }
    }
    const double dt = input.parameters.field("dt").scalar();
    const double duration = input.parameters.field("duration").scalar() * 1000;
    require(dt > 0 && std::isfinite(duration) && duration >= 100, "Use positive param.dt and param.duration >= 0.1 seconds.");
    const double stepCount = std::floor(duration / dt);
    require(std::isfinite(stepCount) && stepCount >= 2 &&
                stepCount < static_cast<double>(std::numeric_limits<std::size_t>::max() / (sizeof(double) * (ne + ni))),
            "The simulation time-step count is outside the supported range.");
    if (!mxIsEmpty(initial)) {
        scalarStruct(initial, "init");
        input.initial = structure();
        for (const char* name : {"ve", "vi"}) {
            input.initial.fields.emplace(name,
                numericNode(field(initial, name, "init"), std::string("init.") + name));
        }
    }
    for (const char* name : {"exe", "exi", "conn_matrix_e", "conn_matrix_i"}) {
        const bool graph = std::strncmp(name, "conn_", 5) == 0;
        const mxArray* value = field(external, name, "ex");
        require(!graph || mxIsSparse(value), std::string("ex.") + name + " must be sparse.");
        input.external.fields.emplace(name, numericNode(value, std::string("ex.") + name, graph, false));
    }
    return input;
}
Mx convertResult(v1standalone::Result& result) {
    auto output = owned(mxCreateStructMatrix(1, 1, 8, resultFields));
    for (const char* name : resultFields) {
        auto& value = result.arrays.at(name);
        require(value.columns == 0 || value.rows <= std::numeric_limits<std::size_t>::max() / value.columns,
                "Native result dimensions overflow.");
        require(value.rows * value.columns == value.data.size() &&
                    value.data.size() <= std::numeric_limits<std::size_t>::max() / sizeof(double),
                "Native result shape or byte count is invalid.");
        auto array = owned(mxCreateUninitNumericMatrix(value.rows, value.columns, mxDOUBLE_CLASS, mxREAL));
        if (!value.data.empty()) std::memcpy(mxGetData(array.get()), value.data.data(), value.data.size() * sizeof(double));
        set(output.get(), name, std::move(array));
        // Free each C++ allocation immediately; do not retain two complete
        // ~2 GB output structs during conversion of the full saved case.
        std::vector<double>().swap(value.data);
    }
    return output;
}
Mx convertStats(const v1standalone::Result& result, const v1standalone::Options& options,
                double inputSeconds, double outputSeconds) {
    const char* names[] = {
        "simulation_seconds", "input_conversion_seconds", "output_conversion_seconds", "total_seconds",
        "seed", "rng", "e_spikes", "i_spikes", "duration_seconds", "dt_ms", "steps", "section_seconds"
    };
    auto stats = owned(mxCreateStructMatrix(1, 1, 12, names));
    scalar(stats.get(), "simulation_seconds", result.simulationSeconds);
    scalar(stats.get(), "input_conversion_seconds", inputSeconds);
    scalar(stats.get(), "output_conversion_seconds", outputSeconds);
    scalar(stats.get(), "total_seconds", 0);
    integer(stats.get(), "seed", options.seed);
    set(stats.get(), "rng", owned(mxCreateString(v1standalone::Rng::algorithm)));
    integer(stats.get(), "e_spikes", result.eSpikes);
    integer(stats.get(), "i_spikes", result.iSpikes);
    scalar(stats.get(), "duration_seconds", result.durationSeconds);
    scalar(stats.get(), "dt_ms", result.dtMilliseconds);
    scalar(stats.get(), "steps", static_cast<double>(result.steps));
    auto sections = owned(mxCreateStructMatrix(1, 1, options.profile ? 12 : 0, options.profile ? sectionFields : nullptr));
    if (options.profile)
        for (std::size_t j = 0; j < 12; ++j) scalar(sections.get(), sectionFields[j], result.sections[j]);
    set(stats.get(), "section_seconds", std::move(sections));
    return stats;
}
void run(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    const auto begin = Clock::now();
    require(nrhs == 3 || nrhs == 4, "Usage: [res, stats] = v1_4_native_mex(param, init, ex, options).");
    require(nlhs <= 2, "At most two outputs are supported: res and stats.");
    require(avx2Available(), "This MEX requires an AVX2-capable x86-64 CPU and operating system.");
    const auto options = readOptions(nrhs == 4 ? prhs[3] : nullptr);
    v1standalone::Result result;
    double inputSeconds;
    {
        auto input = readInputs(prhs[0], prhs[1], prhs[2]);
        inputSeconds = elapsed(begin);
        result = v1standalone::simulate(input.parameters, input.initial, input.external, options);
        // Release input copies before allocating MATLAB result arrays.
    }
    const auto outputStart = Clock::now();
    auto output = convertResult(result);
    const double outputSeconds = elapsed(outputStart);
    Mx stats;
    if (nlhs > 1) {
        stats = convertStats(result, options, inputSeconds, outputSeconds);
        *mxGetDoubles(mxGetField(stats.get(), 0, "total_seconds")) = elapsed(begin);
    }
    if (nlhs > 0) plhs[0] = output.release();
    if (nlhs > 1) plhs[1] = stats.release();
}
} // namespace

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    // MATLAB's error routine does not unwind C++ stack frames. Catch into a
    // fixed buffer, leave the catch scope, then report after every RAII owner
    // and caught exception has been destroyed.
    char error[2048] = {};
    try {
        run(nlhs, plhs, nrhs, prhs);
        return;
    } catch (const std::exception& caught) {
        std::strncpy(error, caught.what(), sizeof(error) - 1);
    } catch (...) {
        std::strncpy(error, "Unknown C++ error in the V1 native adapter.", sizeof(error) - 1);
    }
    mexErrMsgIdAndTxt("V1Native:Error", "%s", error);
}

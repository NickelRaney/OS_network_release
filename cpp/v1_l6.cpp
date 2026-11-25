#include "v1_l6.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#if defined(__AVX2__) || defined(_M_AVX2)
#include <immintrin.h>
#define V1_L6_AVX2 1
#endif

namespace v1standalone {
namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::invalid_argument(message);
}
void finiteVector(const std::vector<double>& values, const char* message) {
    for (double v : values) require(std::isfinite(v), message);
}
}

L6::L6(L6Config config)
    : config_(std::move(config)), normLow_(excitatoryCount), normHigh_(excitatoryCount),
      low_(paddedSide * paddedSide), high_(paddedSide * paddedSide),
      tmpLow_(side * paddedSide), tmpHigh_(side * paddedSide),
      filteredLow_(excitatoryCount), filteredHigh_(excitatoryCount),
      indices_(excitatoryCount) {
    require(config_.mask.size() == excitatoryCount, "L6 mask must contain 117 x 117 elements.");
    require(!config_.transform.empty(), "L6 transform must not be empty.");
    require(config_.coefficient.size() == 1 || config_.coefficient.size() == excitatoryCount,
            "L6 coefficient must be scalar or contain 13689 elements.");
    require(config_.facilitationMaximum.size() == 1 ||
            config_.facilitationMaximum.size() == excitatoryCount,
            "L6 facilitation maximum must be scalar or contain 13689 elements.");
    require(config_.mapEI.size() == inhibitoryCount, "L6 mapEI must contain 2401 elements.");
    finiteVector(config_.transform, "L6 transform must be finite.");
    finiteVector(config_.coefficient, "L6 coefficient must be finite.");
    finiteVector(config_.facilitationMaximum, "L6 facilitation maximum must be finite.");
    require(std::isfinite(config_.inhibitoryFacilitation), "L6 inhibitory facilitation must be finite.");
    for (double mask : config_.mask)
        require(mask == 0.0 || mask == 1.0, "L6 mask must be binary.");
    for (std::size_t index : config_.mapEI)
        require(index < excitatoryCount, "L6 mapEI contains an invalid excitatory index.");

    // Match the original outer product and sum(sum(kernel)) normalization.
    std::array<double, width> gaussian{};
    for (std::size_t k = 0; k < width; ++k) {
        const double x = static_cast<double>(k) - 10.0;
        gaussian[k] = std::exp(-(x * x) / 240.0);
    }
    double mass = 0;
    for (std::size_t col = 0; col < width; ++col) {
        double columnMass = 0;
        for (std::size_t row = 0; row < width; ++row) {
            const double v = gaussian[row] * gaussian[col];
            kernel_[row + width * col] = v;
            columnMass += v;
        }
        mass += columnMass;
    }
    for (double& value : kernel_) value /= mass;
    for (std::size_t k = 0; k < width; ++k) {
        vertical_[k] = kernel_[k + width * 10];
        horizontal_[k] = kernel_[10 + width * k] / kernel_[10 + width * 10];
    }
    for (std::size_t j = 0; j < paddedSide; ++j)
        padIndex_[j] = (j + side - 10) % side;
    for (std::size_t col = 0; col < paddedSide; ++col)
        for (std::size_t row = 0; row < paddedSide; ++row) {
            const double m = config_.mask[padIndex_[row] + side * padIndex_[col]];
            low_[row + paddedSide * col] = m;
            high_[row + paddedSide * col] = 1.0 - m;
        }
    convolveDirect(low_, normLow_);
    convolveDirect(high_, normHigh_);
    separableSupported_ = true;
    double maxNormalization = 0;
    for (std::size_t j = 0; j < excitatoryCount; ++j) {
        require(normLow_[j] > 0 && normHigh_[j] > 0,
                "L6 mask produces a zero normalization denominator.");
        if (normLow_[j] < 1e-8 || normLow_[j] > 2 ||
            normHigh_[j] < 1e-8 || normHigh_[j] > 2) separableSupported_ = false;
        maxNormalization = std::max(maxNormalization,
                                    10.0 * (0.89 / normLow_[j] + 0.11 / normHigh_[j]));
    }
    double originalMass = 0, factorVertical = 0, factorHorizontal = 0, difference = 0;
    for (std::size_t k = 0; k < width; ++k) {
        factorVertical += std::abs(vertical_[k]);
        factorHorizontal += std::abs(horizontal_[k]);
    }
    for (std::size_t col = 0; col < width; ++col)
        for (std::size_t row = 0; row < width; ++row) {
            const double v = kernel_[row + width * col];
            originalMass += std::abs(v);
            difference += std::abs(v - vertical_[row] * horizontal_[col]);
        }
    constexpr double eps = std::numeric_limits<double>::epsilon();
    const double factorMass = factorVertical * factorHorizontal;
    // 4096 eps bounds both the direct 441-term and separable 21+21-term sums.
    // Include measured non-separability plus reconstruction/bound rounding.
    const double kernelDifference = difference + 8 * eps * (originalMass + factorMass);
    const double convolutionError = kernelDifference + 4096 * eps * (originalMass + factorMass);
    zErrorScale_ = convolutionError * maxNormalization * (1 + 16 * eps);

    scalarFac_ = std::all_of(config_.facilitationMaximum.begin(), config_.facilitationMaximum.end(),
                            [&](double v) { return v == config_.facilitationMaximum.front(); });
    const std::size_t tables = scalarFac_ ? 1 : excitatoryCount;
    facTable_.resize(tables * 1001);
    for (std::size_t j = 0; j < tables; ++j)
        for (unsigned count = 0; count <= 1000; ++count) {
            double fac = 1.0 + config_.facilitationMaximum[j] /
                              (std::exp(-(static_cast<double>(count) / 3.5 - 40.0) / 4.5) + 1.0);
            if (fac < 1.05) fac = -0.2625 + 1.25 * fac;
            facTable_[j * 1001 + count] = std::max(fac, 1.0);
        }
}

void L6::convolveDirect(const std::vector<double>& input, std::vector<double>& out) const {
    std::fill(out.begin(), out.end(), 0.0);
    // Each output independently accumulates kernel columns then rows. SIMD
    // spans outputs, so it does not change the order of a cell's additions.
    for (std::size_t col = 0; col < side; ++col)
        for (std::size_t kc = 0; kc < width; ++kc)
            for (std::size_t kr = 0; kr < width; ++kr) {
                const double weight = kernel_[kr + width * kc];
                const double* source = input.data() + (col + width - 1 - kc) * paddedSide + width - 1 - kr;
                double* target = out.data() + col * side;
                std::size_t row = 0;
#if V1_L6_AVX2
                const __m256d w = _mm256_set1_pd(weight);
                for (; row + 4 <= side; row += 4) {
                    const __m256d product = _mm256_mul_pd(_mm256_loadu_pd(source + row), w);
                    _mm256_storeu_pd(target + row, _mm256_add_pd(_mm256_loadu_pd(target + row), product));
                }
#endif
                for (; row < side; ++row) target[row] += source[row] * weight;
            }
}

void L6::convolveSeparable() {
    std::fill(tmpLow_.begin(), tmpLow_.end(), 0.0);
    std::fill(tmpHigh_.begin(), tmpHigh_.end(), 0.0);
    for (std::size_t col = 0; col < paddedSide; ++col)
        for (std::size_t k = 0; k < width; ++k) {
            const double* lo = low_.data() + col * paddedSide + width - 1 - k;
            const double* hi = high_.data() + col * paddedSide + width - 1 - k;
            double* outLo = tmpLow_.data() + col * side;
            double* outHi = tmpHigh_.data() + col * side;
            std::size_t row = 0;
#if V1_L6_AVX2
            const __m256d w = _mm256_set1_pd(vertical_[k]);
            for (; row + 4 <= side; row += 4) {
                _mm256_storeu_pd(outLo + row, _mm256_add_pd(_mm256_loadu_pd(outLo + row),
                    _mm256_mul_pd(_mm256_loadu_pd(lo + row), w)));
                _mm256_storeu_pd(outHi + row, _mm256_add_pd(_mm256_loadu_pd(outHi + row),
                    _mm256_mul_pd(_mm256_loadu_pd(hi + row), w)));
            }
#endif
            for (; row < side; ++row) {
                outLo[row] += lo[row] * vertical_[k];
                outHi[row] += hi[row] * vertical_[k];
            }
        }
    std::fill(filteredLow_.begin(), filteredLow_.end(), 0.0);
    std::fill(filteredHigh_.begin(), filteredHigh_.end(), 0.0);
    for (std::size_t col = 0; col < side; ++col)
        for (std::size_t k = 0; k < width; ++k) {
            const double* lo = tmpLow_.data() + (col + width - 1 - k) * side;
            const double* hi = tmpHigh_.data() + (col + width - 1 - k) * side;
            double* outLo = filteredLow_.data() + col * side;
            double* outHi = filteredHigh_.data() + col * side;
            std::size_t row = 0;
#if V1_L6_AVX2
            const __m256d w = _mm256_set1_pd(horizontal_[k]);
            for (; row + 4 <= side; row += 4) {
                _mm256_storeu_pd(outLo + row, _mm256_add_pd(_mm256_loadu_pd(outLo + row),
                    _mm256_mul_pd(_mm256_loadu_pd(lo + row), w)));
                _mm256_storeu_pd(outHi + row, _mm256_add_pd(_mm256_loadu_pd(outHi + row),
                    _mm256_mul_pd(_mm256_loadu_pd(hi + row), w)));
            }
#endif
            for (; row < side; ++row) {
                outLo[row] += lo[row] * horizontal_[k];
                outHi[row] += hi[row] * horizontal_[k];
            }
        }
}

void L6::probabilities(const double* shortCounts, double* pE, double* pI, bool forceDirect) {
    require(shortCounts && pE && pI, "L6 probability buffers must not be null.");
    ++evaluations_;
    bool supported = separableSupported_ && !forceDirect;
    double maxRate = 0;
    for (std::size_t j = 0; j < excitatoryCount; ++j) {
        const double count = shortCounts[j];
        require(std::isfinite(count) && count >= 0, "L6 short history counts must be finite and nonnegative.");
        supported = supported && count <= 70 && count == std::floor(count);
        maxRate = std::max(maxRate, count / 1.4 * 200.0);
    }
    for (std::size_t col = 0; col < paddedSide; ++col)
        for (std::size_t row = 0; row < paddedSide; ++row) {
            const std::size_t source = padIndex_[row] + side * padIndex_[col];
            const double rate = shortCounts[source] / 1.4 * 200.0;
            const std::size_t target = row + paddedSide * col;
            low_[target] = rate * config_.mask[source];
            high_[target] = rate * (1.0 - config_.mask[source]);
        }
    bool fallback = !supported;
    if (supported) {
        convolveSeparable();
        double maxZ = 1;
        // Use pE temporarily for the pre-lookup rates, avoiding another buffer.
        for (std::size_t j = 0; j < excitatoryCount; ++j) {
            pE[j] = (0.89 * filteredLow_[j] / normLow_[j] +
                     0.11 * filteredHigh_[j] / normHigh_[j]) * 10.0;
            maxZ = std::max(maxZ, pE[j]);
        }
        const double radius = zErrorScale_ * maxRate +
                              512 * std::numeric_limits<double>::epsilon() * maxZ;
        for (std::size_t j = 0; j < excitatoryCount; ++j) {
            const double rounded = std::round(pE[j]);
            if (std::abs(pE[j] - rounded) >= 0.5 - radius) { fallback = true; break; }
        }
    }
    if (fallback) {
        ++fallbackCalls_;
        convolveDirect(low_, filteredLow_);
        convolveDirect(high_, filteredHigh_);
        for (std::size_t j = 0; j < excitatoryCount; ++j)
            pE[j] = (0.89 * filteredLow_[j] / normLow_[j] +
                     0.11 * filteredHigh_[j] / normHigh_[j]) * 10.0;
    }
    for (std::size_t j = 0; j < excitatoryCount; ++j) {
        const double index = std::round(pE[j]);
        if (!std::isfinite(index) || index < 0 || index >= static_cast<double>(config_.transform.size()))
            throw std::out_of_range("L6 transform index outside available table at excitatory neuron " +
                                    std::to_string(j) + ": " + std::to_string(index));
        indices_[j] = static_cast<std::size_t>(index);
        pE[j] = config_.transform[indices_[j]] *
                config_.coefficient[config_.coefficient.size() == 1 ? 0 : j] / 10000.0;
    }
    for (std::size_t j = 0; j < inhibitoryCount; ++j)
        pI[j] = config_.transform[indices_[config_.mapEI[j]]] * 100.0 / 10000.0;
}

double L6::facilitation(std::size_t eNeuron, unsigned windowCount) const {
    if (eNeuron >= excitatoryCount || windowCount > 1000)
        throw std::out_of_range("L6 facilitation history must be an integer in [0, 1000].");
    return facTable_[(scalarFac_ ? 0 : eNeuron) * 1001 + windowCount];
}

} // namespace v1standalone

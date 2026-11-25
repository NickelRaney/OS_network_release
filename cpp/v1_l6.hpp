#pragma once

#include <array>
#include <cstddef>
#include <vector>

namespace v1standalone {

// All matrices follow MATLAB column-major order; mapEI uses zero-based indices.
struct L6Config {
    std::vector<double> mask;
    std::vector<double> transform;
    std::vector<double> coefficient;
    std::vector<double> facilitationMaximum;
    std::vector<std::size_t> mapEI;
    double inhibitoryFacilitation = 1.0;
};

// Deterministic L6 equations only. The caller owns sampling and spike history.
class L6 {
public:
    static constexpr std::size_t side = 117;
    static constexpr std::size_t excitatoryCount = side * side;
    static constexpr std::size_t inhibitoryCount = 49 * 49;

    explicit L6(L6Config config);
    void probabilities(const double* shortCounts, double* pE, double* pI,
                       bool forceDirect = false);
    double facilitation(std::size_t eNeuron, unsigned windowCount) const;
    double inhibitoryFacilitation() const noexcept {
        return config_.inhibitoryFacilitation;
    }

    const std::vector<std::size_t>& lastIndices() const noexcept { return indices_; }
    std::size_t fallbackCalls() const noexcept { return fallbackCalls_; }
    std::size_t evaluations() const noexcept { return evaluations_; }

private:
    static constexpr std::size_t width = 21;
    static constexpr std::size_t paddedSide = side + width - 1;
    L6Config config_;
    std::array<double, width * width> kernel_{};
    std::array<double, width> vertical_{}, horizontal_{};
    std::vector<double> normLow_, normHigh_, low_, high_, tmpLow_, tmpHigh_;
    std::vector<double> filteredLow_, filteredHigh_, facTable_;
    std::vector<std::size_t> indices_;
    std::array<std::size_t, paddedSide> padIndex_{};
    bool separableSupported_ = false;
    bool scalarFac_ = true;
    double zErrorScale_ = 0.0;
    std::size_t fallbackCalls_ = 0, evaluations_ = 0;

    void convolveDirect(const std::vector<double>& input, std::vector<double>& out) const;
    void convolveSeparable();
};

} // namespace v1standalone

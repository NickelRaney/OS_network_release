#pragma once
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#if defined(_MSC_VER)
#define V1RNG_INLINE __forceinline
#else
#define V1RNG_INLINE inline __attribute__((always_inline))
#endif

// MATLAB-independent random stream. See RNG_README.md and vendor/
// for primary-source attribution and licenses. This does not reproduce MATLAB.
namespace v1standalone {
namespace tables {
#include "vendor/ziggurat_constants.h"
}

class Rng {
    static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559,
                  "Rng requires IEEE-754 binary64 double");
    std::uint64_t s_[4];
    static V1RNG_INLINE std::uint64_t rotl(std::uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }
    static std::uint64_t splitmix(std::uint64_t& state) {
        std::uint64_t z = (state += UINT64_C(0x9e3779b97f4a7c15));
        z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
        return z ^ (z >> 31);
    }
public:
    static constexpr const char* algorithm = "xoshiro256**-1.0/splitmix64/numpy-ziggurat-2.3.3";
    explicit Rng(std::uint64_t seed) {
        for (auto& word : s_) word = splitmix(seed);
    }
    // Blackman/Vigna xoshiro256** 1.0; public-domain reference in vendor/.
    V1RNG_INLINE std::uint64_t next() {
        const auto s0 = s_[0], s1 = s_[1], s2 = s_[2], s3 = s_[3];
        const auto result = rotl(s1 * 5, 7) * 9;
        const auto t = s1 << 17;
        const auto t2 = s2 ^ s0, t3 = s3 ^ s1;
        s_[0] = s0 ^ t3; s_[1] = s1 ^ t2;
        s_[2] = t2 ^ t; s_[3] = rotl(t3, 45);
        return result;
    }
    V1RNG_INLINE double uniform() {
        // All 53 mantissa bits, exactly [0,1); never returns 1.
        return static_cast<double>(static_cast<std::int64_t>(next() >> 11)) * 0x1.0p-53;
    }
    // Adapted from NumPy 2.3.3 random_standard_normal, BSD-3-Clause.
    // next_uint64/next_double and npy_log1p use local C++ equivalents; the
    // sign branch uses a bitwise equivalent XOR. Acceptance is unchanged.
    V1RNG_INLINE double normal() {
        using namespace tables;
        for (;;) {
            std::uint64_t r = next();
            const int idx = static_cast<int>(r & 0xff);
            r >>= 8;
            const std::uint64_t rabs = (r >> 1) & UINT64_C(0x000fffffffffffff);
            double x = static_cast<double>(static_cast<std::int64_t>(rabs)) * wi_double[idx];
            // The reference's random sign branch can be a 50% misprediction.
            // Set the IEEE-754 sign bit instead; negative zero is retained.
            std::uint64_t bits;
            std::memcpy(&bits, &x, sizeof bits);
            bits ^= (r & 1) << 63;
            std::memcpy(&x, &bits, sizeof x);
            if (rabs < ki_double[idx]) return x;
            if (idx == 0) {
                for (;;) {
                    const double xx = -ziggurat_nor_inv_r * std::log1p(-uniform());
                    const double yy = -std::log1p(-uniform());
                    if (yy + yy > xx * xx)
                        return ((rabs >> 8) & 1) ? -(ziggurat_nor_r + xx)
                                                : ziggurat_nor_r + xx;
                }
            } else if (((fi_double[idx - 1] - fi_double[idx]) * uniform() +
                         fi_double[idx]) < std::exp(-0.5 * x * x)) {
                return x;
            }
        }
    }
    void fillNormal(double* out, std::size_t n) {
        for (std::size_t i = 0; i < n; ++i) out[i] = normal();
    }
    void fillUniform(double* out, std::size_t n) {
        auto s0 = s_[0], s1 = s_[1], s2 = s_[2], s3 = s_[3];
        for (std::size_t i = 0; i < n; ++i) {
            const auto result = rotl(s1 * 5, 7) * 9;
            const auto t = s1 << 17;
            s2 ^= s0; s3 ^= s1;
            s1 ^= s2; s0 ^= s3;
            s2 ^= t; s3 = rotl(s3, 45);
            out[i] = static_cast<double>(static_cast<std::int64_t>(result >> 11)) * 0x1.0p-53;
        }
        s_[0] = s0; s_[1] = s1; s_[2] = s2; s_[3] = s3;
    }
};

} // namespace v1standalone
#undef V1RNG_INLINE


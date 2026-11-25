#pragma once
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <map>
#include <set>
#include <string>
#include <vector>

// Small MATLAB level-5 reader/writer. No MATLAB libraries are linked.
// Arrays use MATLAB column-major order. Scalar structs, real numeric arrays,
// logical arrays, and CSC sparse matrices are supported (including zlib data).
namespace mat5 {
struct Node {
    enum Kind { Numeric, Struct, Sparse };
    Kind kind = Numeric;
    std::vector<std::size_t> dims;
    std::vector<double> values;
    std::vector<std::uint8_t> bytes; // Only logical arrays use this storage.
    std::map<std::string, Node> fields;
    std::vector<std::uint32_t> ir, jc;
    bool logical = false;

    const Node& field(const std::string& name) const;
    Node& field(const std::string& name);
    std::size_t numel() const;
    double doubleAt(std::size_t index) const;
    double scalar() const;
};
using Variables = std::map<std::string, Node>;

// Empty wanted means read every variable. Requested missing names are errors.
// Unrequested compressed variables are skipped after inspecting their header.
Variables read(const std::filesystem::path& path,
               const std::set<std::string>& wanted = {});

struct FieldView {
    std::string name;
    std::vector<std::size_t> dims;
    const double* data = nullptr;
    std::size_t count = 0;
};
// Writes one scalar struct of real double arrays, directly from caller memory.
// The complete struct must fit MATLAB's 2 GiB per-variable compatibility
// limit (INT32_MAX payload bytes); larger results error before file creation.
void writeStruct(const std::filesystem::path& path,
                 const std::string& variableName,
                 const std::vector<FieldView>& fields);
} // namespace mat5

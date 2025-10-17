#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <iomanip>

struct TestResult {
    std::string name; // the name of the test
    bool passed;   // whether is the test passed
    std::string message;
};

inline std::vector<TestResult> __results;

inline void __record(const std::string& name, bool passed, const std::string& msg = "") {
    __results.push_back({name, passed, msg});
}

#define TEST_CASE(name) void name()

#define RUN_TEST(name) do { \
    try { \
        name(); \
        __record(#name, true); \
    } catch (const std::exception &e) { \
        __record(#name, false, e.what()); \
    } catch (...) { \
        __record(#name, false, "Unknown exception"); \
    } \
} while(0) // 

#define ASSERT_NEAR(val1, val2, tol) do { \
    double _a = (val1), _b = (val2); \
    if (std::fabs(_a - _b) > (tol)) { \
        std::ostringstream oss; \
        oss << std::setprecision(10) \
            << "ASSERT_NEAR failed: |" << #val1 << " - " << #val2 << "| = " \
            << std::fabs(_a - _b) << " > " << tol; \
        throw std::runtime_error(oss.str()); \
    } \
} while(0) // 

#define ASSERT_EQ(val1, val2) do { \
    if ((val1) != (val2)) { \
        std::ostringstream oss; \
        oss << "ASSERT_EQ failed: " << #val1 << " != " << #val2; \
        throw std::runtime_error(oss.str()); \
    } \
} while(0)

inline void printTestSummary() {
    int passed = 0, failed = 0;
    for (auto &r : __results) {
        if (r.passed) ++passed;
        else ++failed;
    }

    std::cout << "\n=============================\n";
    std::cout << "🧪  TEST SUMMARY\n";
    std::cout << "-----------------------------\n";
    for (auto &r : __results) {
        if (r.passed)
            std::cout << "✅  " << r.name << "\n";
        else
            std::cout << "❌  " << r.name << " → " << r.message << "\n";
    }
    std::cout << "-----------------------------\n";
    std::cout << "TOTAL: " << (passed + failed)
              << " | PASSED: " << passed
              << " | FAILED: " << failed << "\n";
    std::cout << "=============================\n\n";

    if (failed > 0) std::exit(1);
}

#include "tests.hpp"

int main() {
    try {
        run_all_tests();
    } catch (const std::exception& e) {
        std::cerr << "❌ Test failed: " << e.what() << "\n";
        return 1;
    }
    return 0;
}

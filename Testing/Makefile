
CXX       = g++
CXXFLAGS  = -std=c++17 -O2 -Wall -Wextra -pedantic
INCLUDES  = -I.                # project's root 

GEN_SRCS  = generated_force.cpp generated_gravity.cpp
TEST_SRCS = $(wildcard test*.cpp)
BIN_DIR   = build
TEST_BINS = $(patsubst %.cpp,$(BIN_DIR)/%,$(TEST_SRCS))

.PHONY: all test clean
all:
	clang++ -std=c++17 -O2 -Wall -Wextra -pedantic \
		test_main.cpp generated_force.cpp generated_gravity.cpp -o tests
	./tests


# All tests 
test: $(TEST_BINS)
	@set -e; for t in $(TEST_BINS); do echo "▶ $$t"; "$$t"; echo ""; done
	@echo "✅ All tests completed successfully."

# --- Her test dosyasını ayrı ayrı derle ---
$(BIN_DIR)/%: %.cpp $(GEN_SRCS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

clean:
	rm -rf $(BIN_DIR) *.o

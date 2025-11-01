CXX := g++ -std=c++17
CXXFLAGS := -std=c++17 -O3 -Wall -Wextra -march=native \
            -I. -IFunctors -Igenerated_files \
            -DMIE_SAFE

TARGET := benchmark
SRCS := main.cpp $(wildcard generated_files/*.cpp)
OBJS := $(SRCS:.cpp=.o)

TEST_SRC := testing.cpp
TEST_TARGET := testing

# Varsayılan hedef: önce test, sonra build, sonra run
all: testcheck $(TARGET)
	@echo "🚀 Running ./$(TARGET)..."
	@./$(TARGET)

# --- test executable ---
$(TEST_TARGET): $(TEST_SRC)
	$(CXX) $(CXXFLAGS) $(TEST_SRC) $(wildcard generated_files/*.cpp) -o $(TEST_TARGET)

# --- test kontrolü ---
testcheck: $(TEST_TARGET)
	@echo "Running pre-check tests..."
	@./$(TEST_TARGET) | tee test_output.txt
	@val=$$(grep "max|ΔF|" test_output.txt | grep -Eo "[0-9]+\.[0-9]+e[+-]?[0-9]+"); \
	thr=1e-12; \
	echo "→ Detected max|ΔF| = $$val (threshold $$thr)"; \
	awk -v v=$$val -v t=$$thr 'BEGIN {if (v+0 <= t) exit 0; else exit 1}' && \
		echo "Tests passed. Building benchmark..." || \
		( echo "Tests failed! benchmark build aborted."; exit 1 )


DELTA := $(shell ./testing | awk '/^\[MAXDELTA\]/{print $$2}')
# --- main target ---
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJS) $(TEST_TARGET) test_output.txt

.PHONY: all clean testcheck

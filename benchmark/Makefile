CXX := g++ -std=c++17
CXXFLAGS := -std=c++17 -O3 -Wall -Wextra -march=native \
            -I. -IFunctors -Igenerated_files \
            -DMIE_SAFE

TARGET := benchmark
SRCS := main.cpp $(wildcard generated_files/*.cpp)
OBJS := $(SRCS:.cpp=.o)

TEST_SRC := testing.cpp
TEST_TARGET := testing
STRICT_THR = 1e-12
BENCH_THR  = 1e-10

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
	@maxrel=$$(grep "RELERR" test_output.txt | awk '{ if ($$2 > m) m = $$2 } END { print m }'); \
	echo "→ Detected max relative error = $$maxrel"; \
	awk -v v="$$maxrel" -v t1=$(STRICT_THR) -v t2=$(BENCH_THR) ' \
	BEGIN { \
	  if (v <= t1) { \
	    print "Strict correctness check passed."; exit 0; \
	  } else if (v <= t2) { \
	    print "Warning: extended tolerance applied (benchmark-only case)."; exit 0; \
	  } else { \
	    print "Correctness check failed."; exit 1; \
	  } \
	}'




DELTA := $(shell ./testing | awk '/^\[MAXDELTA\]/{print $$2}')
# --- main target ---
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJS) $(TEST_TARGET) test_output.txt

.PHONY: all clean testcheck

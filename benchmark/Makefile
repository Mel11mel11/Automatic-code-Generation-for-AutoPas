CXX := g++ -std=c++17
CXXFLAGS := -std=c++17 -O3 -Wall -Wextra -march=native \
            -I. -IFunctors -Igenerated_files \
            -DMIE_SAFE            # ← 

TARGET := benchmark
SRCS := main.cpp $(wildcard generated_files/*.cpp)
OBJS := $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJS)

.PHONY: all clean

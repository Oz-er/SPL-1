# Compiler
CXX := g++

# Build modes (switchable)
BUILD ?= release

# Flags
CXXFLAGS_DEBUG := -Wall -Wextra -std=c++17 -O0 -g
CXXFLAGS_RELEASE := -Wall -Wextra -std=c++17 -O3

# Choose flags based on build type
ifeq ($(BUILD),debug)
    CXXFLAGS := $(CXXFLAGS_DEBUG)
else
    CXXFLAGS := $(CXXFLAGS_RELEASE)
endif

# Target executable
TARGET := DomainAdaption

# Source files
SRCS := Main.cpp TCA.cpp BDA.cpp CORAL.cpp

# Object files (auto-generated)
OBJS := $(SRCS:.cpp=.o)

# Default rule
all: $(TARGET)

# Linking step
$(TARGET): $(OBJS)
	@echo "🔗 Linking..."
	$(CXX) $(OBJS) -o $(TARGET)

# Compilation step (pattern rule)
%.o: %.cpp
	@echo "⚙️ Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build files
clean:
	@echo "🧹 Cleaning..."
	rm -f $(OBJS) $(TARGET)

# Rebuild everything
rebuild: clean all

# Run program
run: all
	./$(TARGET)
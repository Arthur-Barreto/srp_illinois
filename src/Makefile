# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -Wextra -O3

# Source and object files
SRC_DIR = src
OBJ_DIR = obj
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# Targets
TARGETS = read_csv.cpp

# Default target
all: $(TARGETS)

# Create object directory if it doesn't exist
$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)

# Compile global_search without OpenMP
read_csv: $(OBJ_DIR)/read_csv.o
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $^ -o $@

# Compile object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -c $< -o $@

# Clean up build artifacts
clean:
	rm -rf $(OBJ_DIR) $(TARGETS)
# Helpful references:
# https://makefiletutorial.com
# https://stackoverflow.com/a/30602701
# https://stackoverflow.com/a/27838026 

SRC_DIR := src
TESTS_SRC_DIR := tests
VERIFIERS_SRC_DIR := verifiers
VECTORS_SRC_DIR := vectors

BUILD_DIR := build
OBJ_DIR := $(BUILD_DIR)/obj
TESTS_OBJ_DIR := $(OBJ_DIR)/tests
TESTS_BIN_DIR := $(BUILD_DIR)/tests
TESTS_OUT_DIR := $(TESTS_BIN_DIR)/out
VECTORS_DIR := $(BUILD_DIR)/vectors


# Variables for general file 
SRC := $(wildcard $(SRC_DIR)/*.c)
OBJ := $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC))

# TESTS_NAMES := $(foreach token,$(wildcard $(TESTS_DIR)/test_*.c),$(basename $(token)))
TESTS := $(wildcard $(TESTS_SRC_DIR)/test_*.c)
TESTS_BIN := $(patsubst $(TESTS_SRC_DIR)/%.c,$(TESTS_BIN_DIR)/%,$(TESTS))
TESTS_OBJ := $(patsubst $(TESTS_SRC_DIR)/%.c,$(TESTS_OBJ_DIR)/%.o,$(TESTS))
TESTS_OUT := $(patsubst $(TESTS_SRC_DIR)/%.c,$(TESTS_OUT_DIR)/%.out,$(TESTS))

# List of .sage scripts for verification
VERIFIERS := $(wildcard $(VERIFIERS_SRC_DIR)/verify_*.sage)
VECTORS := $(patsubst $(VERIFIERS_SRC_DIR)/%.sage,$(VECTORS_DIR)/%.out,$(VERIFIERS))

# Optional CPP flags: -MMD -MP 
CPPFLAGS := -Iinclude 
CFLAGS   := -Wall -Wextra
LDFLAGS  := -Llib
LDLIBS   := -lgmp

.PHONY: all clean check vectors

all: tests

tests: $(TESTS_BIN)

# Create object files based on the 'src/*.c' files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	@echo "Compiling: '$@'"
	@$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# Create the directories if they do not exist
$(TESTS_BIN_DIR) $(OBJ_DIR) $(TESTS_OBJ_DIR) $(TESTS_OUT_DIR) $(VECTORS_DIR):
	@echo "Creating Directory: '$@'"
	@mkdir -p $@

# Create obj/tests/test_name_something.o <-- from tests/test_name_something.c
# $(TESTS_OBJ_DIR) is order-only prerequisite, meaning it must exist before compiling, but change does not impose rebuilding
$(TESTS_OBJ_DIR)/%.o: $(TESTS_SRC_DIR)/%.c | $(TESTS_OBJ_DIR)
	@echo "Compiling test: '$@'"
	@$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

# Create the executable 'test_something' based on the compiled 'src/*.c' -> 'obj/*.c' objects and 'obj/tests/test_something.o'
$(TESTS_BIN_DIR)/%: $(TESTS_OBJ_DIR)/%.o $(OBJ) | $(TESTS_BIN_DIR)
	@echo "Linking: '$@'"
	@$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

# benches:

# Run unit tests to see if return codes are 0 
check: $(TESTS_OUT)

# Run each of the tests, send output to .out file
$(TESTS_OUT_DIR)/%.out: $(TESTS_BIN_DIR)/% FORCE | $(TESTS_OUT_DIR)
	@echo "------------------------------------"
	@echo "> Run: $(notdir $<)"
	@echo "------------------------------------"
	@$< 1> $@ 

vectors: $(VECTORS) | $(VECTORS_DIR)

# Run sage scripts to generate the test vectors
$(VECTORS_DIR)/%.out: $(VERIFIERS_SRC_DIR)/%.sage FORCE | $(VECTORS_DIR)
	@echo "------------------------------------"
	@echo "> Run: $(notdir $<)"
	@echo "------------------------------------"
	@(sage $< > $@ && diff -q $@ $(subst $(VECTORS_DIR),$(VECTORS_SRC_DIR),$@) --color && echo "OK: Output equal to Test Vectors") || echo "ERROR: Output differ from Test Vectors"

clean:
	@echo "Removing: $(BUILD_DIR)"
	@$(RM) -r $(BUILD_DIR)

# https://www.gnu.org/software/make/manual/make.html#Force-Targets
FORCE:

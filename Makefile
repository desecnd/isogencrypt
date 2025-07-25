# Helpful references:
# https://makefiletutorial.com
# https://stackoverflow.com/a/30602701
# https://stackoverflow.com/a/27838026 

# Run make DEBUG=0 to turn off debugging build
DEBUG ?= 1

SRC_DIR := src

BUILD_DIR := build
OBJ_DIR := $(BUILD_DIR)/obj
OUT_DIR := $(BUILD_DIR)/out
TMP_DIR := $(BUILD_DIR)/tmp

TESTS_SRC_DIR := tests
TESTS_OBJ_DIR := $(OBJ_DIR)/tests
TESTS_BIN_DIR := $(BUILD_DIR)/tests
TESTS_OUT_DIR := $(OUT_DIR)/tests

BENCHES_SRC_DIR := benches
BENCHES_OBJ_DIR := $(OBJ_DIR)/benches
BENCHES_BIN_DIR := $(BUILD_DIR)/benches
BENCHES_OUT_DIR := $(OUT_DIR)/benches

VECTORS_SRC_DIR := assets/test_vectors
DIFFS_OUT_DIR := $(OUT_DIR)/diffs

EXAMPLE_SRC_DIR := example
EXAMPLE_OBJ_DIR := $(OBJ_DIR)/example
EXAMPLE_BIN_DIR := $(BUILD_DIR)/example

# Variables for general file 
SRC := $(wildcard $(SRC_DIR)/*.c)
OBJ := $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC))

# TESTS_NAMES := $(foreach token,$(wildcard $(TESTS_DIR)/test_*.c),$(basename $(token)))
TESTS := $(wildcard $(TESTS_SRC_DIR)/test_*.c)
TESTS_BIN := $(patsubst $(TESTS_SRC_DIR)/%.c,$(TESTS_BIN_DIR)/%,$(TESTS))
TESTS_OBJ := $(patsubst $(TESTS_SRC_DIR)/%.c,$(TESTS_OBJ_DIR)/%.o,$(TESTS))
TESTS_OUT := $(patsubst $(TESTS_SRC_DIR)/%.c,$(TESTS_OUT_DIR)/%.out,$(TESTS))
TESTS_TMP := $(patsubst $(TESTS_SRC_DIR)/%.c,$(TESTS_OUT_DIR)/%.tmp,$(TESTS))

# Diffs must match up with each of the test_vectors
VECTORS := $(wildcard $(VECTORS_SRC_DIR)/*.out)
DIFFS_OUT := $(patsubst $(VECTORS_SRC_DIR)/%.out,$(DIFFS_OUT_DIR)/%.diff,$(VECTORS))

# Same for benchmarks
BENCHES := $(wildcard $(BENCHES_SRC_DIR)/bench_*.c)
BENCHES_BIN := $(patsubst $(BENCHES_SRC_DIR)/%.c,$(BENCHES_BIN_DIR)/%,$(BENCHES))
BENCHES_OBJ := $(patsubst $(BENCHES_SRC_DIR)/%.c,$(BENCHES_OBJ_DIR)/%.o,$(BENCHES))
BENCHES_OUT := $(patsubst $(BENCHES_SRC_DIR)/%.c,$(BENCHES_OUT_DIR)/%.out,$(BENCHES))

# Build example
EXAMPLE := $(wildcard $(EXAMPLE_SRC_DIR)/*.c)
EXAMPLE_BIN := $(patsubst $(EXAMPLE_SRC_DIR)/%.c,$(EXAMPLE_BIN_DIR)/%,$(EXAMPLE))
EXAMPLE_OBJ := $(patsubst $(EXAMPLE_SRC_DIR)/%.c,$(EXAMPLE_OBJ_DIR)/%.o,$(EXAMPLE))

# Optional CPP flags: -MMD -MP 
CPPFLAGS := -Iinclude 
LDFLAGS  := -Llib

# Compilation flags for debug build and release build
# -pg: add profiler data (gprof)
# -g: generate debugging symbols
# -O0: do not optimize 
# -fno-inline: do not inline functions (symbols should be present)
# -fsanitize=address: enable asan
ifeq ($(DEBUG),1)
	LDLIBS := -lgmp -pg -fsanitize=address
	CFLAGS := -Wall -Wextra -O0 -pg -g -fno-inline
else
	LDLIBS := -lgmp
	CFLAGS := -Wall -Wextra -O2
endif


.PHONY: tests benches example all clean run-tests run-diffs
# This allows for calling run-diffs without running run-tests
.NOTINTERMEDIATE: $(TESTS_OUT)

all: tests benches example

# Create object files based on the 'src/*.c' files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	@echo "Compiling: '$@'"
	@$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# Create the directories if they do not exist
$(TESTS_BIN_DIR) $(OBJ_DIR) $(TESTS_OBJ_DIR) $(TESTS_OUT_DIR) $(VECTORS_DIR) $(BENCHES_OBJ_DIR) $(BENCHES_BIN_DIR) $(DIFFS_OUT_DIR) $(EXAMPLE_BIN_DIR) $(EXAMPLE_OBJ_DIR):
	@echo "Creating Directory: '$@'"
	@mkdir -p $@

example: $(EXAMPLE_BIN)

$(EXAMPLE_OBJ_DIR)/%.o: $(EXAMPLE_SRC_DIR)/%.c | $(EXAMPLE_OBJ_DIR)
	@echo "Compiling example: '$@'"
	@$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

$(EXAMPLE_BIN_DIR)/%: $(EXAMPLE_OBJ_DIR)/%.o $(OBJ) | $(EXAMPLE_BIN_DIR)
	@echo "Linking example: '$@'"
	@$(CC) $(LDFLAGS) $^ $(LDLIBS) -lcrypto -lssl -o $@

tests: $(TESTS_BIN)

# Create obj/tests/test_name_something.o <-- from tests/test_name_something.c
# $(TESTS_OBJ_DIR) is order-only prerequisite, meaning it must exist before compiling, but change does not impose rebuilding
$(TESTS_OBJ_DIR)/%.o: $(TESTS_SRC_DIR)/%.c | $(TESTS_OBJ_DIR)
	@echo "Compiling test: '$@'"
	@$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

# Create the executable 'test_something' based on the compiled 'src/*.c' -> 'obj/*.c' objects and 'obj/tests/test_something.o'
$(TESTS_BIN_DIR)/%: $(TESTS_OBJ_DIR)/%.o $(OBJ) | $(TESTS_BIN_DIR)
	@echo "Linking test: '$@'"
	@$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

# benches:
benches: $(BENCHES_BIN)

# Create obj/benches/bench_name_something.o <-- from benches/bench_name_something.c
# $(TESTS_OBJ_DIR) is order-only prerequisite, meaning it must exist before compiling, but change does not impose rebuilding
$(BENCHES_OBJ_DIR)/%.o: $(BENCHES_SRC_DIR)/%.c | $(BENCHES_OBJ_DIR)
	@echo "Compiling benchmark: '$@'"
	@$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

# Create the executable 'bench_something' based on the compiled 'src/*.c' -> 'obj/*.c' objects and 'obj/tests/test_something.o'
$(BENCHES_BIN_DIR)/%: $(BENCHES_OBJ_DIR)/%.o $(OBJ) | $(BENCHES_BIN_DIR)
	@echo "Linking benchmark: '$@'"
	@$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

# Ephemeral target to "always" run
run-tests: $(TESTS_TMP)

# Target for "temporary" files, but actually crate .out files
# It must be this way, ohterwise: run-diff will always rerun the tests (run-tests) (not wanted)
$(TESTS_OUT_DIR)/%.tmp: $(TESTS_BIN_DIR)/% FORCE | $(TESTS_OUT_DIR)
	@echo "------------------------------------"
	@echo "> Run: $(notdir $<)"
	@echo "------------------------------------"
	@$< > $(patsubst $(TESTS_OUT_DIR)/%.tmp,$(TESTS_OUT_DIR)/%.out,$@)

run-diffs: $(DIFFS_OUT) 

# Target for real "out" files, in order to recreate them as neccessary 
$(TESTS_OUT_DIR)/%.out: $(TESTS_BIN_DIR)/%  | $(TESTS_OUT_DIR)
	@echo "------------------------------------"
	@echo "> Run for diffs: $(notdir $<)"
	@echo "------------------------------------"
	$< > $@

# Run each of the tests, send output to .out file
# @echo "> $(notdir $@)"
# && echo "Check diff: $(notdir $@) (\033[0;32mPASSED\033[0m)" || (echo "[\033[0;31mFAILED\033[0m]: $(notdir $@) (see: $@)")
$(DIFFS_OUT_DIR)/%.diff: $(TESTS_OUT_DIR)/%.out $(VECTORS_SRC_DIR)/%.out $(TESTS_BIN_DIR)/% FORCE | $(DIFFS_OUT_DIR)
	@diff -Z $< $(patsubst $(DIFFS_OUT_DIR)/%.diff,$(VECTORS_SRC_DIR)/%.out,$@) --color > $@ \
	&& echo "Check diff: $(notdir $@) (\033[0;32mPASSED\033[0m)" || (echo "Check diff: $(notdir $@) (\033[0;31mFAILED\033[0m) - see: $@")

clean:
	@echo "Removing: $(BUILD_DIR)"
	@$(RM) -r $(BUILD_DIR)

# https://www.gnu.org/software/make/manual/make.html#Force-Targets
FORCE:

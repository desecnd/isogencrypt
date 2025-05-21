# Isogencrypt

C library for isogeny-based cryptographic research. 

Requires `GMP` (Gnu Multiprecision Library) https://gmplib.org/.

## Compilation

```bash
# Compile all unit test executables
make tests

# Compile all benchmark executables
make benches
```

## Testing and Verifying

```bash
# Run all unit tests under `tests` directory. Compile if neccessary
$ make run-tests

# Compare output of .c unit tests with prepared test_vectors in assets dir. 
# Execute 'run-tests' target if .out files are missing in build/tests.
# Compile if necessary. Results of diffs are stored in build directory for inspection
$ make run-diffs
```

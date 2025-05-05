# Isogencrypt

C library for isogeny-based cryptographic research. 

Requires `GMP` (Gnu Multiprecision Library) https://gmplib.org/.

## Near Roadmap

### Theoretical
- [x] Implement `fp2_mul_int` for easier mulitplying complex a+bi by scalar c 
- [x] Implement `A_from_A24p()` and `A24p_from_A()` and write tests for it
- [x] Implement `xLADDER()` algorithm 
    * [x] Test `xLADDER_int()` algorithm
    * [x] Implement `xLADDER()` for mpz_t
- [x] Implement `iso_curve()` function for the codomain in the Montgomery Model
    * [x] Implemented `aISOG_from_KPS` method
    * [x] Write `verify.sage` tests and `ec_mont` tests
- [x] Implement `pow2-isogeny` and test it
- [x] Implement chain isogeny
    * [x] Implement 2-isogeny codomain curve
    * [x] Implement 2-isogeny point image
    * [x] Add **even**-degree isogeny chaining
    * [x] Add list of "push-through" points
- [x] Add j-invariant calculation in montgomery model
- [ ] Add MSIDH with initialized torsion basis
    * [x] Add msidh_gen_pubkey
    * [ ] Add msidh_key_exchange
- [ ] Implement mechanism behind the order2 of the point (x_push isog)
    * implement a_from_alpha function
    * implement root alpha finding (?)


### Technical

- [ ] Implement benchmark mechanism
- [ ] Implement `fp2_equal` for comparsion and replace all calls to `mpz_cmp_ui`
- [ ] Implement **debug** mechanism (when to print specific message of values)
- [ ] Find out the best way to create `const` typedef pointers (make point_t not a const pointer but const point)
- [ ] Replace `assert` in some functions to "error" checking
- [ ] Add `two_exp` ? or something that means number of 2 powers into `pprod_t` (simple counter of the power of 2), this will make some of the algorithms easier and more optimized
- [ ] Make prints during tests and verifies as "Debug" and print them on the "error" stream, this way we can compare the results between the Python and C using simple "diff"
- [ ] Add variadic function for fp2 init? also for point_init 
- [ ] Add clangd formating with `.clang-format` 
- [ ] Setup Tests for Sage and [Makefile](https://stackoverflow.com/questions/4927676/implementing-make-check-or-make-test) for comparison with the expected output

## Future

Run benchmarks:
- Test with and without `-O2` flag 
- Compare performance with using `fp` and `fp2` (and `ec_mont`?) registers instead of allocating  
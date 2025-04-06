# Isogencrypt

C library for isogeny-based cryptographic research. 

Requires `GMP` (Gnu Multiprecision Library) https://gmplib.org/.

## Near Roadmap

- [x] Implement `fp2_mul_int` for easier mulitplying complex a+bi by scalar c 
- [x] Implement `A_from_A24p()` and `A24p_from_A()` and write tests for it
- [x] Implement `xLADDER()` algorithm 
    * [x] Test xLADDER algoirhtm
- [x] Implement `iso_curve()` function for the codomain in the Montgomery Model
    * [x] Implemented `aISOG_from_KPS` method
    * [x] Write `verify.sage` tests and `ec_mont` tests
- [x] Implement chain isogeny
    * Add **even**-degree isogeny chaining
    * Add list of "push-through" points
- [ ] Implement mechanism behind the order2 of the point (x_push isog)
    * implement a_from_alpha function
    * implement root alpha finding (?)
- [ ] Implement benchmark mechanism
- [ ] Implement `fp2_equal` for comparsion and replace all calls to `mpz_cmp_ui`
- [ ] Implement **debug** mechanism (when to print specific message of values)

## Future

Run benchmarks:
- Test with and without `-O2` flag 
- Compare performance with using `fp` and `fp2` (and `ec_mont`?) registers instead of allocating  
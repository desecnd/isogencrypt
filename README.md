# Isogencrypt

C library for isogeny-based cryptographic research. 

Requires `GMP` (Gnu Multiprecision Library) https://gmplib.org/.

## Near Roadmap

- [x] Implement `fp2_mul_int` for easier mulitplying complex a+bi by scalar c 
- [x] Implement `A_from_A24p()` and `A24p_from_A()` and write tests for it
- [ ] Implement `fp2_equal` for comparsion and replace all calls to `mpz_cmp_ui`
- [] Implement `iso_curve()` function for the codomain in the Montgomery Model
    * [x] Implemented `aISOG_from_KPS` method
    * [] Write `verify.sage` tests and `ec_mont` tests
- [] Implement mechanism behind the order2 of the point (x_push isog)
    * implement a_from_alpha function
    * implement root alpha finding (?)

## Future

Run benchmarks:
- Test with and without `-O2` flag 
- Compare performance with using `fp` and `fp2` (and `ec_mont`?) registers instead of allocating  
all: test_fp2 test_fp test_ec_mont test_msidh

test_fp2: test_fp2.c fp2.c fp2.h
	gcc test_fp2.c fp2.c -g -o test_fp2 -Wextra -Wall -lgmp

test_fp: test_fp.c fp2.c fp2.h
	gcc test_fp.c fp2.c -g -o test_fp -Wextra -Wall -lgmp

test_ec_mont: test_ec_mont.c fp2.c fp2.h ec_mont.c ec_mont.h pprod.h pprod.c
	gcc fp2.c ec_mont.c test_ec_mont.c pprod.c -g -o test_ec_mont -Wextra -Wall -lgmp

test_msidh: test_msidh.c proto_msidh.c proto_msidh.h pprod.h pprod.c
	gcc fp2.c test_msidh.c proto_msidh.c pprod.c -g -o test_msidh -Wextra -Wall -lgmp

clean: 
	rm -f test_fp test_fp2 test_ec_mont test_msidh

all: test_fp2 test_fp test_mont

test_fp2: test_fp2.c fp2.c fp2.h
	gcc test_fp2.c fp2.c -g -o test_fp2 -Wextra -Wall -lgmp

test_fp: test_fp.c fp2.c fp2.h
	gcc test_fp.c fp2.c -g -o test_fp -Wextra -Wall -lgmp

test_mont: test_mont.c fp2.c fp2.h mont.c mont.h
	gcc fp2.c mont.c test_mont.c -g -o test_mont -Wextra -Wall -lgmp

clean: 
	rm -f test_fp test_fp2 test_mont
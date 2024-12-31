all: test_fp2 test_fp

test_fp2: test_fp2.c fp2.c fp2.h
	gcc test_fp2.c fp2.c -o test_fp2 -Wextra -Wall -lgmp

test_fp: test_fp.c fp2.c fp2.h
	gcc test_fp.c fp2.c -o test_fp -Wextra -Wall -lgmp

clean: 
	rm -f test_fp test_fp2
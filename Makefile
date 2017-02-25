
ldacenc:ldacenc.c
	gcc -O2 -Werror -Ilibldac/inc -o ldacenc ldacenc.c libldac/src/ldaclib.c libldac/src/ldacBT.c
clean:
	rm -f ldacenc.o ldaclib.o ldacBT.o ldacenc


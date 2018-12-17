
all:ldacenc ldacinfo

ldacenc:ldacenc.c
	gcc -O2 -Werror -Ilibldac/inc -o ldacenc ldacenc.c libldac/src/ldaclib.c libldac/src/ldacBT.c

ldacinfo:ldacinfo.c
	gcc -O2 -Werror -o ldacinfo ldacinfo.c

clean:
	rm -f ldacenc.o ldaclib.o ldacBT.o ldacenc ldacinfo


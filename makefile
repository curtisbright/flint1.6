ifdef PLEASE_DO_NOT_ERASE_DAVIDS_64_BIT_FLAG_AGAIN_THANKS_VERY_MUCH
	BITFLAG = -m64
else
	BITFLAG = 
endif

CC = gcc -std=c99 $(BITFLAG)

ifndef FLINT_GMP_INCLUDE_DIR
	FLINT_GMP_INCLUDE_DIR = "/home/dmharvey/gmp/install/include"
endif

ifndef FLINT_GMP_LIB_DIR
	FLINT_GMP_LIB_DIR = "/home/dmharvey/gmp/install/lib"
endif

LIBS = -L$(FLINT_GMP_LIB_DIR) -lgmp -lpthread -lm
INCS =  -I"/usr/include" -I$(FLINT_GMP_INCLUDE_DIR)

CFLAGS = $(INCS) -funroll-loops -fexpensive-optimizations -O3

RM = rm -f

.PHONY: all all-before all-after clean clean-custom driver test

all: all-before ssmul all-after

clean: clean-custom
	${RM} $(FLINTOBJ) $(PROFOBJ) $(TESTOBJ) $(DRIVEROBJ) $(BIN) $(FBIN)

HEADERS = Z.h Z_mpn.h ZmodF.h ZmodFpoly.h Zpoly.h Zpoly_mpn.h \
			 extras.h flint-manager.h flint.h longlong.h longlong_wrapper.h \
			 mpn_extras.h profiler-main.h profiler.h


mpn_extras.o: mpn_extras.c $(HEADERS)
	$(CC) -c mpn_extras.c -o mpn_extras.o $(CFLAGS)

Z.o: Z.c $(HEADERS)
	$(CC) -c Z.c -o Z.o $(CFLAGS)

flint-manager.o: flint-manager.c $(HEADERS)
	$(CC) -c flint-manager.c -o flint-manager.o $(CFLAGS)



profiler.o: profiler.c $(HEADERS)
	$(CC) -c profiler.c -o profiler.o $(CFLAGS)
	
profiler-main.o: profiler-main.c $(HEADERS)
	$(CC) -c profiler-main.c -o profiler-main.o $(CFLAGS)



Z_mpn.o: Z_mpn.c $(HEADERS)
	$(CC) -c Z_mpn.c -o Z_mpn.o $(CFLAGS)



Zpoly.o: Zpoly.c $(HEADERS)
	$(CC) -c Zpoly.c -o Zpoly.o $(CFLAGS)
	
Zpoly-test.o: Zpoly-test.c $(HEADERS)
	$(CC) -c Zpoly-test.c -o Zpoly-test.o $(CFLAGS)
	
Zpoly-test: Zpoly-test.o Zpoly.o flint-manager.o
	$(CC) Zpoly.o Zpoly-test.o flint-manager.o -o Zpoly-test $(CFLAGS) $(LIBS)
	


Zpoly_mpn.o: Zpoly_mpn.c $(HEADERS)
	$(CC) -c Zpoly_mpn.c -o Zpoly_mpn.o $(CFLAGS)
	
Zpoly_mpn-test.o: Zpoly_mpn-test.c $(HEADERS)
	$(CC) -c Zpoly_mpn-test.c -o Zpoly_mpn-test.o $(CFLAGS)
	
Zpoly_mpn-test: Zpoly_mpn-test.o ZmodFpoly.o ZmodF.o Z_mpn.o mpn_extras.o Zpoly.o Zpoly_mpn.o Zpoly.h flint-manager.o 
	$(CC) Zpoly_mpn.o ZmodFpoly.o ZmodF.o Z_mpn.o mpn_extras.o Zpoly.o Zpoly_mpn-test.o flint-manager.o -o Zpoly_mpn-test $(CFLAGS) $(LIBS)

Zpoly_mpn-profile-tables.o: Zpoly_mpn-profile.c $(HEADERS)
	python make-profile-tables.py Zpoly_mpn
	$(CC) -c Zpoly_mpn-profile-tables.c -o Zpoly_mpn-profile-tables.o $(CFLAGS)
	rm Zpoly_mpn-profile-tables.c
	
Zpoly_mpn-profile.o: Zpoly_mpn-profile.c $(HEADERS)
	$(CC) -c Zpoly_mpn-profile.c -o Zpoly_mpn-profile.o $(CFLAGS)

Zpoly_mpn-profile: Z_mpn.o ZmodFpoly.o flint-manager.o Zpoly_mpn-profile.o Zpoly_mpn-profile-tables.o Zpoly_mpn.o profiler-main.o profiler.o ZmodF.o Zpoly_mpn.o Zpoly.o mpn_extras.o
	$(CC) -o Zpoly_mpn-profile Z_mpn.o Zpoly_mpn-profile.o Zpoly_mpn-profile-tables.o profiler.o profiler-main.o Zpoly_mpn.o Zpoly.o flint-manager.o mpn_extras.o ZmodFpoly.o ZmodF.o $(CFLAGS) $(LIBS)



ZmodF.o: ZmodF.c $(HEADERS)
	$(CC) -c ZmodF.c -o ZmodF.o $(CFLAGS)
	
ZmodF-test.o: ZmodF-test.c $(HEADERS)
	$(CC) -c ZmodF-test.c -o ZmodF-test.o $(CFLAGS)

ZmodF-test: ZmodF-test.o ZmodF.o
	$(CC) ZmodF-test.o ZmodF.o -o ZmodF-test $(CFLAGS) $(LIBS)



ZmodFpoly.o: ZmodFpoly.c $(HEADERS)
	$(CC) -c ZmodFpoly.c -o ZmodFpoly.o $(CFLAGS)
	
ZmodFpoly-test.o: ZmodFpoly-test.c $(HEADERS)
	$(CC) -c ZmodFpoly-test.c -o ZmodFpoly-test.o $(CFLAGS)

ZmodFpoly-test: Z_mpn.o Zpoly_mpn.o Zpoly.o flint-manager.o mpn_extras.o ZmodFpoly-test.o ZmodFpoly.o ZmodF.o
	$(CC) Z_mpn.o Zpoly_mpn.o Zpoly.o flint-manager.o mpn_extras.o ZmodFpoly-test.o ZmodFpoly.o ZmodF.o -o ZmodFpoly-test $(CFLAGS) $(LIBS)

ZmodFpoly-profile-tables.o: ZmodFpoly-profile.c $(HEADERS)
	python make-profile-tables.py ZmodFpoly
	$(CC) -c ZmodFpoly-profile-tables.c -o ZmodFpoly-profile-tables.o $(CFLAGS)
	rm ZmodFpoly-profile-tables.c

ZmodFpoly-profile.o: ZmodFpoly-profile.c $(HEADERS)
	$(CC) -c ZmodFpoly-profile.c -o ZmodFpoly-profile.o $(CFLAGS)

ZmodFpoly-profile: Z_mpn.o flint-manager.o ZmodFpoly-profile.o ZmodFpoly-profile-tables.o ZmodFpoly.o profiler-main.o profiler.o ZmodF.o Zpoly_mpn.o Zpoly.o mpn_extras.o
	$(CC) -o ZmodFpoly-profile Z_mpn.o ZmodFpoly-profile.o ZmodFpoly-profile-tables.o profiler.o profiler-main.o Zpoly_mpn.o Zpoly.o flint-manager.o mpn_extras.o ZmodFpoly.o ZmodF.o $(CFLAGS) $(LIBS)
	

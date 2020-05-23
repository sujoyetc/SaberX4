#CC 		  = /usr/bin/gcc-4.8
CC 		  = /usr/bin/gcc
CFLAGS  = -Wall -Wextra -O3 -fomit-frame-pointer -msse2avx -mavx2 -march=native -lcrypto
CFLAGS2  = -Wall -Wextra -O3 -fomit-frame-pointer -msse2avx -mavx2 -march=native -std=c11 
CFLAGS3  = -Wall -Wextra -O3 -fomit-frame-pointer -msse2avx -mavx2 -march=native
CFLAGS4  = -Wall -Wextra -O3 -fomit-frame-pointer -msse2avx -mavx2 -march=native -std=c99 
CLANG   = clang -march=native -O3 -fomit-frame-pointer -lcrypto -fwrapv -Qunused-arguments -w
RM 		  = /bin/rm


all: saberX4Test/PQCgenKAT_kem \
     saberX4Test/test_kex \
     saberX4Test/test_sample_matrix 

SOURCES = pack_unpack.c poly.c verify.c fips202.o fips202x4.c SABER_indcpa.c kem.c \
					cbd.c crypto_stream_aes256ctr.c printBits.c andodoAUX.c\
					keccak4x/KeccakP-1600-times4-SIMD256.o rng.o
#rng.o
#randombytes.c

SOURCES2 = pack_unpack.c rng.o poly.c verify.c fips202.o fips202x4.c SABER_indcpa.c kem.c \
					cbd.c crypto_stream_aes256ctr.c printBits.c andodoAUX.c\
					keccak4x/KeccakP-1600-times4-SIMD256.o

HEADERS = SABER_params.h pack_unpack.h poly.h verify.h cbd.h SABER_indcpa.h kem.h rng.h fips202.h fips202x4.h api.h randombytes.h \
					printBits.h andodoAUX.h


saberX4Test/test_kex: $(SOURCES) $(HEADERS) saberX4Test/test_kex.c 
	$(CC) $(CFLAGS3) -o $@ $(SOURCES) saberX4Test/test_kex.c -lcrypto

saberX4Test/test_sample_matrix: $(SOURCES) $(HEADERS) saberX4Test/sample_matrix_test.c 
	$(CC) $(CFLAGS3) -o $@ $(SOURCES) saberX4Test/sample_matrix_test.c -lcrypto

saberX4Test/PQCgenKAT_kem: $(SOURCES2) $(HEADERS) saberX4Test/PQCgenKAT_kem.c 
	$(CC) $(CFLAGS3) -o $@ $(SOURCES2) saberX4Test/PQCgenKAT_kem.c -lcrypto  
	
# saberX4Test/test_BS2POLVEC: $(SOURCES) $(HEADERS) saberX4Test/test_BS2POLVEC.c 
	# $(CC) $(CFLAGS3) -o $@ $(SOURCES) saberX4Test/test_BS2POLVEC.c -lcrypto


keccak4x/KeccakP-1600-times4-SIMD256.o: keccak4x/KeccakP-1600-times4-SIMD256.c \
					keccak4x/align.h \
					keccak4x/brg_endian.h \
					keccak4x/KeccakP-1600-times4-SIMD256.c \
					keccak4x/KeccakP-1600-times4-SnP.h \
					keccak4x/KeccakP-1600-unrolling.macros \
					keccak4x/SIMD256-config.h

	$(CC) $(CFLAGS) -c keccak4x/KeccakP-1600-times4-SIMD256.c -o $@

rng.o: rng.c
	$(CC) $(CFLAGS2) -c rng.c -lcrypto -o $@ 

fips202.o: fips202.c
	$(CLANG) -c $^ -o $@

.PHONY: clean saberX4Test

saberX4Test:
	./saberX4Test/test_kex
	./saberX4Test/PQCgenKAT_kem
	# ./saberX4Test/test_BS2POLVEC


clean:
	-$(RM) *.o
	-$(RM) -r saberX4Test/test_sample_matrix
	-$(RM) -r saberX4Test/test_kex
	-$(RM) -r saberX4Test/PQCgenKAT_kem
	# -$(RM) -r saberX4Test/test_BS2POLVEC
	-$(RM) -r keccak4x/KeccakP-1600-times4-SIMD256.o



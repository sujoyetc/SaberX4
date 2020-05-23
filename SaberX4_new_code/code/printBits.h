#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <arpa/inet.h>
#include <errno.h>
#include <immintrin.h>
char  buffer[20];

FILE* openfile(char* filename);

void print_byte_as_bits(char val, FILE* fp);

void print_bits(char * ty, char * val, unsigned char * bytes, size_t num_bytes, char* filename);



#define SHOW(T,V, filename) do { T x = V; print_bits(#T, #V, (unsigned char*) &x, sizeof(x), filename); } while(0)



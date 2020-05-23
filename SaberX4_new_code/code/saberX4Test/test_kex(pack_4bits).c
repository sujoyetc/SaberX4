#include "../SABER_indcpa.h"
#include "../kem.h"
#include "../api.h"
#include "../poly.h"
//#include "../randombytes.h"
#include "../rng.h"

#include "../cpucycles.c"
#include "../verify.h"

#include<stdio.h>
#include<stdint.h>
#include<stdlib.h>
#include<time.h>
#include<immintrin.h>
#include<string.h>


uint64_t clock_S, clock_E;
uint64_t clock_S0, clock_E0,clock_S1, clock_E1,clock_S2, clock_E2;
long long int clock_T0,clock_T1,clock_T2, clock_T;
unsigned char *bs0, *bs1, *bs2, *bs3;
uint8_t bytes[128];
uint8_t bytes0[128], bytes1[128], bytes2[128], bytes3[128];


__m256i data0[16], data1[16], data2[16], data3[16];
__m256i avx_epi64_01, avx_epi64_03, avx_epi64_07, avx_epi64_0f, avx_epi64_1f, avx_epi64_3f, avx_epi64_7f, avx_epi64_ff;
uint16_t dataN0[256], dataN1[256], dataN2[256], dataN3[256];
unsigned char* stream;

extern int crypto_kem_keypair(
							unsigned char *pk0, unsigned char *sk0,
							unsigned char *pk1, unsigned char *sk1,
							unsigned char *pk2, unsigned char *sk2,
							unsigned char *pk3, unsigned char *sk3);

//extern int crypto_kem_enc(unsigned char *ct, unsigned char *ss, const unsigned char *pk);
extern int crypto_kem_enc(
									unsigned char *c0, unsigned char *k0, const unsigned char *pk0,
									unsigned char *c1, unsigned char *k1, const unsigned char *pk1,
									unsigned char *c2, unsigned char *k2, const unsigned char *pk2,
									unsigned char *c3, unsigned char *k3, const unsigned char *pk3);


//extern int crypto_kem_dec(unsigned char *ss, const unsigned char *ct, const unsigned char *sk);
extern int crypto_kem_dec(
							unsigned char *ss0, const unsigned char *ct0, const unsigned char *sk0,
							unsigned char *ss1, const unsigned char *ct1, const unsigned char *sk1,
							unsigned char *ss2, const unsigned char *ct2, const unsigned char *sk2,
							unsigned char *ss3, const unsigned char *ct3, const unsigned char *sk3);
void init_avx(){
	avx_epi64_01 = _mm256_set1_epi64x(0x01);
	avx_epi64_03 = _mm256_set1_epi64x(0x03);
	avx_epi64_07 = _mm256_set1_epi64x(0x07);
	avx_epi64_0f = _mm256_set1_epi64x(0x0f);
	avx_epi64_1f = _mm256_set1_epi64x(0x1f);
	avx_epi64_3f = _mm256_set1_epi64x(0x3f);
	avx_epi64_7f = _mm256_set1_epi64x(0x7f);
	avx_epi64_ff = _mm256_set1_epi64x(0xff);
}

unsigned char* gene_bytestream(size_t num_bytes){
  stream = malloc(num_bytes);
  int i;
  for (i = 0; i < num_bytes; i++)
  {
    //stream[i]=255;
    stream[i] = rand ();
  }
  return stream;
}

void gene_4bs(){
  bs0 = gene_bytestream(256*16/8);
  bs1 = gene_bytestream(256*16/8);
  bs2 = gene_bytestream(256*16/8);
  bs3 = gene_bytestream(256*16/8);
}

void test_pack_4bit(uint8_t *bytes, uint16_t *data){
	clock_S = cpucycles();
	uint32_t i,j;
	uint32_t offset_data=0;
	
	for(j=0;j<SABER_N/2;j++)
	{
		offset_data=2*j;
		bytes[j]= (data[offset_data] & 0x0f) | ( (data[offset_data + 1] & 0x0f)<<4 );
	}
	clock_E = cpucycles();
	clock_T += clock_E - clock_S;
}


void test_pack_4bitX4(uint8_t *bytes0, uint8_t *bytes1, uint8_t *bytes2, uint8_t *bytes3, uint16_t *data0, uint16_t *data1, uint16_t *data2, uint16_t *data3){
	//printf("#011 Original bytes0 address is: OX%p\n",bytes0);
	
uint32_t i,j,k;
	uint32_t offset_data=0;
	clock_S0 = cpucycles();
	//uint64_t tmp_data[256][4];
	__m256i avx_data[256];
	for(i=0;i<256;++i){
		avx_data[i] = _mm256_set_epi64x(data3[i], data2[i], data1[i], data0[i]);
	}
	clock_E0 = cpucycles();
	clock_T0 += clock_E0 - clock_S0;

	clock_S1 = cpucycles();
	__m256i avx_bytes[128];
	__m256i avx_epi64_0f = _mm256_set1_epi64x(0x0f);
	for(j=0;j<SABER_N/2;++j)
	{
		offset_data=2*j;
		avx_bytes[j]= _mm256_or_si256((_mm256_and_si256(avx_data[offset_data],avx_epi64_0f)),(_mm256_slli_epi64(_mm256_and_si256(avx_data[offset_data + 1],avx_epi64_0f),4)));
	}
	clock_E1 = cpucycles();
	clock_T1 += clock_E1 - clock_S1;

	clock_S2 = cpucycles();
	__m256i *avx_byte0_ptr = (__m256i*)bytes0;
	__m256i *avx_byte1_ptr = (__m256i*)bytes1;
	__m256i *avx_byte2_ptr = (__m256i*)bytes2;
	__m256i *avx_byte3_ptr = (__m256i*)bytes3;
	
	//k<4
	uint32_t offset_bytes = 0;
	//256/8=everytime store 32 value
	//totally have SABER_N/2 = 128 value need to be stored
	//(SABER_N/2)/(256/8)
	
	__m256i temp_byte0, temp_byte1, temp_byte2, temp_byte3;
	//printf("#AVX pointer size is: %d\n",sizeof(avx_byte0_ptr));
	for(k=0;k<4;++k){

		offset_bytes = (256/8)*k;

		temp_byte0 = _mm256_set_epi8(avx_bytes[offset_bytes + 31][0], avx_bytes[offset_bytes + 30][0], avx_bytes[offset_bytes + 29][0], avx_bytes[offset_bytes + 28][0],
						avx_bytes[offset_bytes + 27][0], avx_bytes[offset_bytes + 26][0], avx_bytes[offset_bytes + 25][0], avx_bytes[offset_bytes + 24][0],
						avx_bytes[offset_bytes + 23][0], avx_bytes[offset_bytes + 22][0], avx_bytes[offset_bytes + 21][0], avx_bytes[offset_bytes + 20][0],
						avx_bytes[offset_bytes + 19][0], avx_bytes[offset_bytes + 18][0], avx_bytes[offset_bytes + 17][0], avx_bytes[offset_bytes + 16][0],
						avx_bytes[offset_bytes + 15][0], avx_bytes[offset_bytes + 14][0], avx_bytes[offset_bytes + 13][0], avx_bytes[offset_bytes + 12][0],
						avx_bytes[offset_bytes + 11][0], avx_bytes[offset_bytes + 10][0], avx_bytes[offset_bytes + 9][0], avx_bytes[offset_bytes + 8][0],
						avx_bytes[offset_bytes + 7][0], avx_bytes[offset_bytes + 6][0], avx_bytes[offset_bytes + 5][0], avx_bytes[offset_bytes + 4][0],
						avx_bytes[offset_bytes + 3][0], avx_bytes[offset_bytes + 2][0], avx_bytes[offset_bytes + 1][0], avx_bytes[offset_bytes + 0][0]);
				_mm256_storeu_si256(avx_byte0_ptr+k,temp_byte0);
		
		temp_byte1 = _mm256_set_epi8(avx_bytes[offset_bytes + 31][1], avx_bytes[offset_bytes + 30][1], avx_bytes[offset_bytes + 29][1], avx_bytes[offset_bytes + 28][1],
						avx_bytes[offset_bytes + 27][1], avx_bytes[offset_bytes + 26][1], avx_bytes[offset_bytes + 25][1], avx_bytes[offset_bytes + 24][1],
						avx_bytes[offset_bytes + 23][1], avx_bytes[offset_bytes + 22][1], avx_bytes[offset_bytes + 21][1], avx_bytes[offset_bytes + 20][1],
						avx_bytes[offset_bytes + 19][1], avx_bytes[offset_bytes + 18][1], avx_bytes[offset_bytes + 17][1], avx_bytes[offset_bytes + 16][1],
						avx_bytes[offset_bytes + 15][1], avx_bytes[offset_bytes + 14][1], avx_bytes[offset_bytes + 13][1], avx_bytes[offset_bytes + 12][1],
						avx_bytes[offset_bytes + 11][1], avx_bytes[offset_bytes + 10][1], avx_bytes[offset_bytes + 9][1], avx_bytes[offset_bytes + 8][1],
						avx_bytes[offset_bytes + 7][1], avx_bytes[offset_bytes + 6][1], avx_bytes[offset_bytes + 5][1], avx_bytes[offset_bytes + 4][1],
						avx_bytes[offset_bytes + 3][1], avx_bytes[offset_bytes + 2][1], avx_bytes[offset_bytes + 1][1], avx_bytes[offset_bytes + 0][1]);
				_mm256_storeu_si256(avx_byte1_ptr+k,temp_byte1);
						
		temp_byte2 = _mm256_set_epi8(avx_bytes[offset_bytes + 31][2], avx_bytes[offset_bytes + 30][2], avx_bytes[offset_bytes + 29][2], avx_bytes[offset_bytes + 28][2],
						avx_bytes[offset_bytes + 27][2], avx_bytes[offset_bytes + 26][2], avx_bytes[offset_bytes + 25][2], avx_bytes[offset_bytes + 24][2],
						avx_bytes[offset_bytes + 23][2], avx_bytes[offset_bytes + 22][2], avx_bytes[offset_bytes + 21][2], avx_bytes[offset_bytes + 20][2],
						avx_bytes[offset_bytes + 19][2], avx_bytes[offset_bytes + 18][2], avx_bytes[offset_bytes + 17][2], avx_bytes[offset_bytes + 16][2],
						avx_bytes[offset_bytes + 15][2], avx_bytes[offset_bytes + 14][2], avx_bytes[offset_bytes + 13][2], avx_bytes[offset_bytes + 12][2],
						avx_bytes[offset_bytes + 11][2], avx_bytes[offset_bytes + 10][2], avx_bytes[offset_bytes + 9][2], avx_bytes[offset_bytes + 8][2],
						avx_bytes[offset_bytes + 7][2], avx_bytes[offset_bytes + 6][2], avx_bytes[offset_bytes + 5][2], avx_bytes[offset_bytes + 4][2],
						avx_bytes[offset_bytes + 3][2], avx_bytes[offset_bytes + 2][2], avx_bytes[offset_bytes + 1][2], avx_bytes[offset_bytes + 0][2]);
				_mm256_storeu_si256(avx_byte2_ptr+k,temp_byte2);
						
		temp_byte3 = _mm256_set_epi8(avx_bytes[offset_bytes + 31][3], avx_bytes[offset_bytes + 30][3], avx_bytes[offset_bytes + 29][3], avx_bytes[offset_bytes + 28][3],
						avx_bytes[offset_bytes + 27][3], avx_bytes[offset_bytes + 26][3], avx_bytes[offset_bytes + 25][3], avx_bytes[offset_bytes + 24][3],
						avx_bytes[offset_bytes + 23][3], avx_bytes[offset_bytes + 22][3], avx_bytes[offset_bytes + 21][3], avx_bytes[offset_bytes + 20][3],
						avx_bytes[offset_bytes + 19][3], avx_bytes[offset_bytes + 18][3], avx_bytes[offset_bytes + 17][3], avx_bytes[offset_bytes + 16][3],
						avx_bytes[offset_bytes + 15][3], avx_bytes[offset_bytes + 14][3], avx_bytes[offset_bytes + 13][3], avx_bytes[offset_bytes + 12][3],
						avx_bytes[offset_bytes + 11][3], avx_bytes[offset_bytes + 10][3], avx_bytes[offset_bytes + 9][3], avx_bytes[offset_bytes + 8][3],
						avx_bytes[offset_bytes + 7][3], avx_bytes[offset_bytes + 6][3], avx_bytes[offset_bytes + 5][3], avx_bytes[offset_bytes + 4][3],
						avx_bytes[offset_bytes + 3][3], avx_bytes[offset_bytes + 2][3], avx_bytes[offset_bytes + 1][3], avx_bytes[offset_bytes + 0][3]);
				_mm256_storeu_si256(avx_byte3_ptr+k,temp_byte3);
	}
	clock_E2 = cpucycles();
	clock_T2 += clock_E2 - clock_S2;
}



int test_kem_cca()
{


  uint8_t pk0[SABER_PUBLICKEYBYTES], sk0[SABER_SECRETKEYBYTES];
  uint8_t pk1[SABER_PUBLICKEYBYTES], sk1[SABER_SECRETKEYBYTES];
  uint8_t pk2[SABER_PUBLICKEYBYTES], sk2[SABER_SECRETKEYBYTES];
  uint8_t pk3[SABER_PUBLICKEYBYTES], sk3[SABER_SECRETKEYBYTES];


	uint8_t c0[SABER_BYTES_CCA_DEC], c1[SABER_BYTES_CCA_DEC], c2[SABER_BYTES_CCA_DEC], c3[SABER_BYTES_CCA_DEC];	
	uint8_t k_a0[SABER_KEYBYTES], k_b0[SABER_KEYBYTES], k_a1[SABER_KEYBYTES], k_b1[SABER_KEYBYTES];
	uint8_t k_a2[SABER_KEYBYTES], k_b2[SABER_KEYBYTES], k_a3[SABER_KEYBYTES], k_b3[SABER_KEYBYTES];

  unsigned char entropy_input[48];
	
  uint64_t i, j, repeat;
  repeat=10000;
  //repeat = 1;

  uint64_t CLOCK1,CLOCK2;
  uint64_t CLOCK_kp,CLOCK_enc,CLOCK_dec;

  	CLOCK1 = 0;
        CLOCK2 = 0;
	CLOCK_kp = CLOCK_enc = CLOCK_dec = 0;
        clock_arith = clock_samp = clock_load = 0;

	time_t t;
   	// Intializes random number generator
   	srand((unsigned) time(&t));

    	for (i=0; i<48; i++){
        	entropy_input[i] = i;
        	//entropy_input[i] = rand()%256;
	}
    	randombytes_init(entropy_input, NULL, 256);


  	for(i=0; i<repeat; i++)
  	{
	    //printf("i : %lu\n",i);

	    //Generation of secret key sk and public key pk pair
	    CLOCK1=cpucycles();	
	    crypto_kem_keypair(pk0,sk0, pk1,sk1, pk2,sk2, pk3,sk3);
	    CLOCK2=cpucycles();	
	    CLOCK_kp=CLOCK_kp+(CLOCK2-CLOCK1);	
	
	    //Key-Encapsulation call; input: pk; output: ciphertext c, shared-secret k_a;	
	    CLOCK1=cpucycles();
			crypto_kem_enc(
										c0, k_a0, pk0,
										c1, k_a1, pk1,
										c2, k_a2, pk2,
										c3, k_a3, pk3);


	    CLOCK2=cpucycles();	
	    CLOCK_enc=CLOCK_enc+(CLOCK2-CLOCK1);	

	
	    //Key-Decapsulation call; input: sk, c; output: shared-secret k_b;	
	    CLOCK1=cpucycles();
			crypto_kem_dec(
										k_b0, c0, sk0,
										k_b1, c1, sk1,
										k_b2, c2, sk2,
										k_b3, c3, sk3);

	    CLOCK2=cpucycles();	
	    CLOCK_dec=CLOCK_dec+(CLOCK2-CLOCK1);	
  

	    		
	    // Functional verification: check if k_a == k_b?
	    for(j=0; j<SABER_KEYBYTES; j++)
	    {

				if(k_a0[j] != k_b0[j])
				{
					printf("----- ERR CCA KEM0 ------\n");
					return 0;		
					break;
				}
				if(k_a1[j] != k_b1[j])
				{
					printf("----- ERR CCA KEM1 ------\n");
					return 0;		
					break;
				}
				if(k_a2[j] != k_b2[j])
				{
					printf("----- ERR CCA KEM2 ------\n");
					return 0;		
					break;
				}
				if(k_a3[j] != k_b3[j])
				{
					printf("----- ERR CCA KEM3 ------\n");
					return 0;		
					break;
				}


	    }
   		
  	}
	
	
	printf("Repeat is : %ld\n",repeat);
	printf("Average times key_pair: \t %lu \n",CLOCK_kp/repeat);

	printf("Average times enc: \t %lu \n",CLOCK_enc/repeat);
	printf("Average times dec: \t %lu \n",CLOCK_dec/repeat);
	printf("Average time sample_matrix: \t %lu \n",clock_matrix/repeat);
	printf("Average times sample_secret: \t %lu \n",clock_secret/repeat);
	printf("Average times polynomial mul: \t %lu \n",clock_mul/(3*repeat));
	printf("Average times polynomial mul: \t %lu \n",clock_mul/(3*count_mul));
	printf("Number of times polynomial mul: \t %lu \n",count_mul);


  	return 0;
}



int main()
{
	
	//test_kem_cca();
	init_avx();
	gene_4bs();
	int i;
	for(i=0;i<10000;i++){
		test_pack_4bit(bs0,dataN0);
		test_pack_4bit(bs1,dataN1);
		test_pack_4bit(bs2,dataN2);
		test_pack_4bit(bs3,dataN3);
		test_pack_4bitX4(bs0,bs1,bs2,bs3,dataN0,dataN1,dataN2,dataN3);
	}
	printf("Average load: \t %lu \n",clock_T0/10000);
	printf("Average operate: \t %lu \n",clock_T1/10000);
	printf("Average store: \t %lu \n",clock_T2/10000);
	printf("Average X1: \t %lu \n",clock_T/10000);
	
	return 0;
}

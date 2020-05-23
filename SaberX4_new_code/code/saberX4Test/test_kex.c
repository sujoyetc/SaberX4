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
  repeat=1000000;

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

  	return 0;
}



int main()
{
	test_kem_cca();
	return 0;
}

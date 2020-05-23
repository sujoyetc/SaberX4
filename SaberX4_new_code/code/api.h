//
//  api.h
//
//  Created by Bassham, Lawrence E (Fed) on 9/6/17.
//  Copyright © 2017 Bassham, Lawrence E (Fed). All rights reserved.
//


//   This is a sample 'api.h' for use 'sign.c'

#ifndef api_h
#define api_h


// Available algorithms for different security levels
#define Saber 2

// Change the algorithm name 
#define SABER_TYPE Saber

//  Set these three values apropriately for your algorithm
#if SABER_TYPE == Saber
	#define CRYPTO_ALGNAME "Saber"
	#define CRYPTO_SECRETKEYBYTES 2304
	#define CRYPTO_PUBLICKEYBYTES (3*320+32)
	#define CRYPTO_BYTES 32
	#define CRYPTO_CIPHERTEXTBYTES 1088
	#define Saber_type 2
#endif
int crypto_kem_keypair(
						unsigned char *pk0, unsigned char *sk0,
						unsigned char *pk1, unsigned char *sk1,
						unsigned char *pk2, unsigned char *sk2,
						unsigned char *pk3, unsigned char *sk3);

extern int crypto_kem_enc(
									unsigned char *c0, unsigned char *k0, const unsigned char *pk0,
									unsigned char *c1, unsigned char *k1, const unsigned char *pk1,
									unsigned char *c2, unsigned char *k2, const unsigned char *pk2,
									unsigned char *c3, unsigned char *k3, const unsigned char *pk3);


int crypto_kem_dec(
				   	unsigned char *ss0, const unsigned char *ct0, const unsigned char *sk0,
					unsigned char *ss1, const unsigned char *ct1, const unsigned char *sk1,
					unsigned char *ss2, const unsigned char *ct2, const unsigned char *sk2,
					unsigned char *ss3, const unsigned char *ct3, const unsigned char *sk3);

#endif /* api_h */

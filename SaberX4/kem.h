#ifndef INDCPA_H
#define INDCPA_H

//void indcpa_keypair(uint8_t *pk, uint8_t *sk);
void indcpa_keypair(
					uint8_t *pk0, uint8_t *sk0,
					uint8_t *pk1, uint8_t *sk1,
					uint8_t *pk2, uint8_t *sk2,
					uint8_t *pk3, uint8_t *sk3);

void indcpa_client(uint8_t *pk, uint8_t *b_prime, uint8_t *c, uint8_t *key);

void indcpa_server(uint8_t *pk, uint8_t *b_prime, uint8_t *c, uint8_t *key);

void indcpa_kem_keypair( 
						uint8_t *pk0, uint8_t *sk0,
						uint8_t *pk1, uint8_t *sk1,
						uint8_t *pk2, uint8_t *sk2,
						uint8_t *pk3, uint8_t *sk3);

//void indcpa_kem_enc(uint8_t *message, uint8_t *noiseseed, uint8_t *pk,  uint8_t *ciphertext);
void indcpa_kem_enc(unsigned char *m0, unsigned char *m1, unsigned char *m2, unsigned char *m3,
				    const unsigned char *pk0, unsigned char *c0,
					const unsigned char *pk1, unsigned char *c1,
					const unsigned char *pk2, unsigned char *c2,
					const unsigned char *pk3, unsigned char *c3);


void indcpa_kem_dec(uint8_t *sk, uint8_t *ciphertext, uint8_t message_dec[]);

int crypto_kem_keypair(
						unsigned char *pk0, unsigned char *sk0,
						unsigned char *pk1, unsigned char *sk1,
						unsigned char *pk2, unsigned char *sk2,
						unsigned char *pk3, unsigned char *sk3);

//int crypto_kem_enc(unsigned char *c, unsigned char *k, const unsigned char *pk);
int crypto_kem_enc(
									unsigned char *c0, unsigned char *k0, const unsigned char *pk0,
									unsigned char *c1, unsigned char *k1, const unsigned char *pk1,
									unsigned char *c2, unsigned char *k2, const unsigned char *pk2,
									unsigned char *c3, unsigned char *k3, const unsigned char *pk3);

//int crypto_kem_dec(unsigned char *k, const unsigned char *c, const unsigned char *sk);
int crypto_kem_dec(
				   	unsigned char *k0, const unsigned char *c0, const unsigned char *sk0,
					unsigned char *k1, const unsigned char *c1, const unsigned char *sk1,
					unsigned char *k2, const unsigned char *c2, const unsigned char *sk2,
					unsigned char *k3, const unsigned char *c3, const unsigned char *sk3);

uint64_t clock1,clock2;
uint64_t clock_kp_mv,clock_cl_mv, clock_kp_sm, clock_cl_sm;


#endif


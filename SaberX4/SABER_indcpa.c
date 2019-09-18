#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "api.h"
#include "SABER_indcpa.h"
#include "pack_unpack.h"
//#include "randombytes.h"
#include "rng.h"
#include "cbd.h"
#include "SABER_params.h"
#include "./polymul/toom_cook_4/toom-cook_4way.c"
#include "fips202.h"
#include "fips202x4.h"

#define h1 4 //2^(EQ-EP-1)

#define h2 ( (1<<(SABER_EP-2)) - (1<<(SABER_EP-SABER_ET-1)) + (1<<(SABER_EQ-SABER_EP-1)) )



uint64_t mask_ar[4]={0xFFFFFFFFUL,0xFFFFFFFFUL,0xFFFFFFFFUL,0xFFFFFFFFUL};
__m256i mask_load;
__m256i floor_round;
__m256i H1_avx;	
__m256i H2_avx;

void POL2MSG(uint16_t *message_dec_unpacked, unsigned char *message_dec);

/*--------------------------------------------------------------------------------------
	This routine loads the constant values for Toom-Cook multiplication 
----------------------------------------------------------------------------------------*/
void load_values(){

	int64_t i;

	int64_t inv3=43691;
	int64_t inv9=36409;
	int64_t inv15=61167;

	int64_t int45=45;
	int64_t int30=30;
	int64_t int0=0;


	int16_t inv3_avx_load[16],inv9_avx_load[16],inv15_avx_load[16],int45_avx_load[16],int30_avx_load[16],int0_avx_load[16];

	for(i=0;i<16;i++){
		inv3_avx_load[i]=inv3;
		inv9_avx_load[i]=inv9;
		inv15_avx_load[i]=inv15;
		int45_avx_load[i]=int45;
		int30_avx_load[i]=int30;
		int0_avx_load[i]=int0;
	}

	inv3_avx = _mm256_loadu_si256 ((__m256i const *) (&inv3_avx_load));
	inv9_avx = _mm256_loadu_si256 ((__m256i const *) (&inv9_avx_load));
	inv15_avx = _mm256_loadu_si256 ((__m256i const *) (&inv15_avx_load));
	int45_avx = _mm256_loadu_si256 ((__m256i const *) (&int45_avx_load));
	int30_avx = _mm256_loadu_si256 ((__m256i const *) (&int30_avx_load));
	int0_avx = _mm256_loadu_si256 ((__m256i const *) (&int0_avx_load));
	mask = _mm256_loadu_si256 ((__m256i const *)mask_ar);	
}



/*-----------------------------------------------------------------------------------
	This routine generates a=[Matrix K x K] of 256-coefficient polynomials 
-------------------------------------------------------------------------------------*/



void BS2POLq(const unsigned char *bytes, uint16_t data[SABER_N]){
	
	uint32_t j;
	uint32_t offset_data=0,offset_byte=0;	
	
	offset_byte=0;

		for(j=0;j<SABER_N/8;j++){
			offset_byte=13*j;
			offset_data=8*j;
			data[offset_data + 0]= ( bytes[ offset_byte + 0 ] & (0xff)) | ((bytes[offset_byte + 1] & 0x1f)<<8);
			data[offset_data + 1]= ( bytes[ offset_byte + 1 ]>>5 & (0x07)) | ((bytes[offset_byte + 2] & 0xff)<<3) | ((bytes[offset_byte + 3] & 0x03)<<11);
			data[offset_data + 2]= ( bytes[ offset_byte + 3 ]>>2 & (0x3f)) | ((bytes[offset_byte + 4] & 0x7f)<<6);
			data[offset_data + 3]= ( bytes[ offset_byte + 4 ]>>7 & (0x01)) | ((bytes[offset_byte + 5] & 0xff)<<1) | ((bytes[offset_byte + 6] & 0x0f)<<9);
			data[offset_data + 4]= ( bytes[ offset_byte + 6 ]>>4 & (0x0f)) | ((bytes[offset_byte + 7] & 0xff)<<4) | ((bytes[offset_byte + 8] & 0x01)<<12);
			data[offset_data + 5]= ( bytes[ offset_byte + 8]>>1 & (0x7f)) | ((bytes[offset_byte + 9] & 0x3f)<<7);
			data[offset_data + 6]= ( bytes[ offset_byte + 9]>>6 & (0x03)) | ((bytes[offset_byte + 10] & 0xff)<<2) | ((bytes[offset_byte + 11] & 0x07)<<10);
			data[offset_data + 7]= ( bytes[ offset_byte + 11]>>3 & (0x1f)) | ((bytes[offset_byte + 12] & 0xff)<<5);
		}


}


void GenMatrix(polyvec *a, const unsigned char *seed) 
{
  unsigned int one_vector=13*SABER_N/8;
  unsigned int byte_bank_length=SABER_K*SABER_K*one_vector;
  unsigned char buf[byte_bank_length];

  uint16_t temp_ar[SABER_N];

  int i,j,k;
  uint16_t mod = (SABER_Q-1);


  shake128(buf,byte_bank_length,seed,SABER_SEEDBYTES);


  for(i=0;i<SABER_K;i++)
  {
    for(j=0;j<SABER_K;j++)
    {
			BS2POLq(buf+(i*SABER_K+j)*one_vector,temp_ar);
			for(k=0;k<SABER_N;k++){
				a[i].vec[j].coeffs[k] = (temp_ar[k])& mod ;
			}
    }
  }
}

void GenMatrix4x(polyvec *a0, polyvec *a1, polyvec *a2, polyvec *a3, 
								const unsigned char *seed0, const unsigned char *seed1, 
								const unsigned char *seed2, const unsigned char *seed3) 
{
  unsigned int one_vector=13*SABER_N/8;
  unsigned int byte_bank_length=SABER_K*SABER_K*one_vector;

  unsigned char buf0[byte_bank_length], buf1[byte_bank_length], buf2[byte_bank_length], buf3[byte_bank_length];

  int i,j;

  shake128x4(buf0,buf1,buf2,buf3, byte_bank_length, seed0, seed1, seed2, seed3, SABER_SEEDBYTES);

  for(i=0;i<SABER_K;i++)
  {
    for(j=0;j<SABER_K;j++)
    {
			//BS2POLq(buf0+(i*SABER_K+j)*one_vector,temp_ar);
			BS2POLq(buf0+(i*SABER_K+j)*one_vector, &a0[i].vec[j].coeffs[0]);
			BS2POLq(buf1+(i*SABER_K+j)*one_vector, &a1[i].vec[j].coeffs[0]);
			BS2POLq(buf2+(i*SABER_K+j)*one_vector, &a2[i].vec[j].coeffs[0]);
			BS2POLq(buf3+(i*SABER_K+j)*one_vector, &a3[i].vec[j].coeffs[0]);

    }
  }
}

void GenSecret(uint16_t r[SABER_K][SABER_N],const unsigned char *seed){


		uint32_t i;

		int32_t buf_size= SABER_MU*SABER_N*SABER_K/8;

		uint8_t buf[buf_size];

		shake128(buf, buf_size, seed,SABER_NOISESEEDBYTES);

		for(i=0;i<SABER_K;i++)
		{
			cbd(r[i],buf+i*SABER_MU*SABER_N/8);
		}
}

void GenSecret4x(uint16_t r0[SABER_K][SABER_N], uint16_t r1[SABER_K][SABER_N], 
								uint16_t r2[SABER_K][SABER_N], uint16_t r3[SABER_K][SABER_N],
								const unsigned char *seed0, const unsigned char *seed1,
								const unsigned char *seed2, const unsigned char *seed3)
{

		uint32_t i;

		int32_t buf_size= SABER_MU*SABER_N*SABER_K/8;

		uint8_t buf0[buf_size], buf1[buf_size], buf2[buf_size], buf3[buf_size];

  	shake128x4(buf0, buf1, buf2, buf3, buf_size, seed0, seed1, seed2, seed3, SABER_SEEDBYTES);

		for(i=0;i<SABER_K;i++)
		{
			cbd(r0[i],buf0+i*SABER_MU*SABER_N/8);
			cbd(r1[i],buf1+i*SABER_MU*SABER_N/8);
			cbd(r2[i],buf2+i*SABER_MU*SABER_N/8);
			cbd(r3[i],buf3+i*SABER_MU*SABER_N/8);
		}
}

void indcpa_kem_keypair(unsigned char *pk0, unsigned char *sk0,
												unsigned char *pk1, unsigned char *sk1,
												unsigned char *pk2, unsigned char *sk2,
												unsigned char *pk3, unsigned char *sk3 )
{
 
	polyvec a0[SABER_K], a1[SABER_K], a2[SABER_K], a3[SABER_K];

  uint16_t skpv1_0[SABER_K][SABER_N], skpv1_1[SABER_K][SABER_N], skpv1_2[SABER_K][SABER_N], skpv1_3[SABER_K][SABER_N];
	uint16_t (*skpv1_ptr)[SABER_N];
	unsigned char *sk_ptr, *pk_ptr, *seed_ptr;

  unsigned char seed0[SABER_SEEDBYTES], seed1[SABER_SEEDBYTES], seed2[SABER_SEEDBYTES], seed3[SABER_SEEDBYTES];
  unsigned char noiseseed0[SABER_COINBYTES], noiseseed1[SABER_COINBYTES], noiseseed2[SABER_COINBYTES], noiseseed3[SABER_COINBYTES];

  int32_t i,j,k;


	//--------------AVX declaration------------------
	
  __m256i res_avx[SABER_K][SABER_N/16];
  __m256i acc[2*SABER_N/16];

  __m256i sk0_avx[SABER_K][SABER_N/16], sk1_avx[SABER_K][SABER_N/16];
  __m256i sk2_avx[SABER_K][SABER_N/16], sk3_avx[SABER_K][SABER_N/16];
	__m256i a0_avx[SABER_K][SABER_K][SABER_N/16], a1_avx[SABER_K][SABER_K][SABER_N/16];
	__m256i a2_avx[SABER_K][SABER_K][SABER_N/16], a3_avx[SABER_K][SABER_K][SABER_N/16];
	__m256i (*sk_avx_ptr)[SABER_N/16], (*a_avx_ptr)[SABER_K][SABER_N/16];


	mask_ar[0]=~(0UL);mask_ar[1]=~(0UL);mask_ar[2]=~(0UL);mask_ar[3]=~(0UL);
	mask_load = _mm256_loadu_si256 ((__m256i const *)mask_ar);

	floor_round=_mm256_set1_epi16(4);

	H1_avx=_mm256_set1_epi16(h1);


	//--------------AVX declaration ends------------------

	load_values();


	// randombytes0-3 are used to pass the KAT test. randombytes() is a deterministic function in rng.c
	// In a real application randombytes0-3 should be replaced by randombytes()
	randombytes0(seed0, SABER_SEEDBYTES);
	randombytes0(noiseseed0, SABER_COINBYTES);
	randombytes1(seed1, SABER_SEEDBYTES);
	randombytes1(noiseseed1, SABER_COINBYTES);
	randombytes2(seed2, SABER_SEEDBYTES);
	randombytes2(noiseseed2, SABER_COINBYTES);
	randombytes3(seed3, SABER_SEEDBYTES);
	randombytes3(noiseseed3, SABER_COINBYTES);
	

	shake128x4(seed0, seed1, seed2, seed3, SABER_SEEDBYTES, seed0, seed1, seed2, seed3, SABER_SEEDBYTES);

	GenMatrix4x(a0, a1, a2, a3, seed0, seed1, seed2, seed3);

	GenSecret4x(skpv1_0, skpv1_1, skpv1_2, skpv1_3, noiseseed0, noiseseed1, noiseseed2, noiseseed3);


 	// Load sk into avx vectors		
 	for(i=0;i<SABER_K;i++)
 	{
		for(j=0; j<SABER_N/16; j++){
  		//sk_avx[i][j] = _mm256_loadu_si256 ((__m256i const *) (&skpv1[i][j*16]));
  		sk0_avx[i][j] = _mm256_loadu_si256 ((__m256i const *) (&skpv1_0[i][j*16]));
  		sk1_avx[i][j] = _mm256_loadu_si256 ((__m256i const *) (&skpv1_1[i][j*16]));
  		sk2_avx[i][j] = _mm256_loadu_si256 ((__m256i const *) (&skpv1_2[i][j*16]));
  		sk3_avx[i][j] = _mm256_loadu_si256 ((__m256i const *) (&skpv1_3[i][j*16]));
		}
	}

  // Load a into avx vectors	
  for(i=0;i<SABER_K;i++){ 
	  for(j=0;j<SABER_K;j++){
		  for(k=0;k<SABER_N/16;k++){
				//a_avx[i][j][k]=_mm256_loadu_si256 ((__m256i const *) (&a[i].vec[j].coeffs[k*16]));
				a0_avx[i][j][k]=_mm256_loadu_si256 ((__m256i const *) (&a0[i].vec[j].coeffs[k*16]));
				a1_avx[i][j][k]=_mm256_loadu_si256 ((__m256i const *) (&a1[i].vec[j].coeffs[k*16]));
				a2_avx[i][j][k]=_mm256_loadu_si256 ((__m256i const *) (&a2[i].vec[j].coeffs[k*16]));
				a3_avx[i][j][k]=_mm256_loadu_si256 ((__m256i const *) (&a3[i].vec[j].coeffs[k*16]));
		  }
	  }
  }	

	int SABER_SERIAL=0;
	
	//////////////////////////////////////////////////////////////////////////////
	// Serially Generate {sk0,pk0}, {sk1,pk1}, {sk2,pk2}, {sk3,pk3}. 

	for(SABER_SERIAL=0; SABER_SERIAL<4; SABER_SERIAL++)
	{
		if(SABER_SERIAL==0) {sk_avx_ptr = sk0_avx; a_avx_ptr=a0_avx; skpv1_ptr=skpv1_0; sk_ptr=sk0; pk_ptr=pk0; seed_ptr=seed0;}
		if(SABER_SERIAL==1) {sk_avx_ptr = sk1_avx; a_avx_ptr=a1_avx; skpv1_ptr=skpv1_1; sk_ptr=sk1; pk_ptr=pk1; seed_ptr=seed1;}
		if(SABER_SERIAL==2) {sk_avx_ptr = sk2_avx; a_avx_ptr=a2_avx; skpv1_ptr=skpv1_2; sk_ptr=sk2; pk_ptr=pk2; seed_ptr=seed2;}
		if(SABER_SERIAL==3) {sk_avx_ptr = sk3_avx; a_avx_ptr=a3_avx; skpv1_ptr=skpv1_3; sk_ptr=sk3; pk_ptr=pk3; seed_ptr=seed3;}


		for(i=0;i<SABER_K;i++){
			for(j=0;j<SABER_N/16;j++){
				res_avx[i][j]=_mm256_xor_si256(res_avx[i][j],res_avx[i][j]);
			}
		}

		// Matrix-vector multiplication; Matrix in transposed order
		for(i=0;i<SABER_K;i++){
			for(j=0;j<SABER_K;j++){
				toom_cook_4way_avx(a_avx_ptr[j][i], sk_avx_ptr[j], SABER_Q, acc);

				for(k=0;k<SABER_N/16;k++){
					res_avx[i][k]=_mm256_add_epi16(res_avx[i][k],acc[k]);
				}
			}
		}

		// Now truncation
		for(i=0;i<SABER_K;i++){ //shift right EQ-EP bits
			for(j=0;j<SABER_N/16;j++){
				res_avx[i][j]=_mm256_add_epi16 (res_avx[i][j], H1_avx);
				res_avx[i][j]=_mm256_srli_epi16 (res_avx[i][j], (SABER_EQ-SABER_EP) );
			}
		}

		//------------------Pack sk into byte string-------
		
		POLVEC2BS(sk_ptr,skpv1_ptr,SABER_Q);

		//------------------Pack pk into byte string-------
	
		for(i=0;i<SABER_K;i++){ // reuses skpv1[] for unpacking avx of public-key
				for(j=0;j<SABER_N/16;j++){
					_mm256_maskstore_epi32 ((int *) (skpv1_ptr[i]+j*16), mask_load, res_avx[i][j]);
				}
			}
		POLVEC2BS(pk_ptr,skpv1_ptr,SABER_P); // load the public-key into pk byte string 	


		for(i=0;i<SABER_SEEDBYTES;i++){ // now load the seedbytes in PK. Easy since seed bytes are kept in byte format.
			pk_ptr[SABER_POLYVECCOMPRESSEDBYTES + i]=seed_ptr[i]; 
		}
	}
	//////////////////////////////////////////////////////////////////////////////
	// End: Serially Generate {sk0,pk0}, {sk1,pk1}, {sk2,pk2}, {sk3,pk3}. 
	

}


void indcpa_kem_enc(
					unsigned char *message_received0, unsigned char *message_received1,
					unsigned char *message_received2, unsigned char *message_received3,
					unsigned char *noiseseed0, unsigned char *noiseseed1, 
					unsigned char *noiseseed2, unsigned char *noiseseed3,
					const unsigned char *pk0, unsigned char *ciphertext0,
					const unsigned char *pk1, unsigned char *ciphertext1,
					const unsigned char *pk2, unsigned char *ciphertext2,
					const unsigned char *pk3, unsigned char *ciphertext3)
{ 


	uint32_t i,j,k;
	polyvec a0[SABER_K], a1[SABER_K], a2[SABER_K], a3[SABER_K];
	unsigned char seed0[SABER_SEEDBYTES], seed1[SABER_SEEDBYTES], seed2[SABER_SEEDBYTES], seed3[SABER_SEEDBYTES];

	uint16_t pkcl[SABER_K][SABER_N]; 	//public key of received by the client

	uint16_t skpv1_0[SABER_K][SABER_N], skpv1_1[SABER_K][SABER_N], skpv1_2[SABER_K][SABER_N], skpv1_3[SABER_K][SABER_N];

	uint16_t temp[SABER_K][SABER_N];
	uint16_t message[SABER_KEYBYTES*8];

	unsigned char *message_received_ptr, *ciphertext_ptr;
	const unsigned char *pk_ptr;

	//--------------AVX declaration------------------
	
	  __m256i sk0_avx[SABER_K][SABER_N/16], sk1_avx[SABER_K][SABER_N/16], sk2_avx[SABER_K][SABER_N/16], sk3_avx[SABER_K][SABER_N/16];
	  __m256i mod_p;
	  __m256i res_avx[SABER_K][SABER_N/16];
	  __m256i vprime_avx[SABER_N/16];
	  __m256i a0_avx[SABER_K][SABER_K][SABER_N/16], a1_avx[SABER_K][SABER_K][SABER_N/16];
		__m256i a2_avx[SABER_K][SABER_K][SABER_N/16], a3_avx[SABER_K][SABER_K][SABER_N/16];
	  __m256i acc[2*SABER_N/16];

	  __m256i pkcl_avx[SABER_K][SABER_N/16];

		__m256i message_avx[SABER_N/16];

		__m256i (*sk_avx_ptr)[SABER_N/16], (*a_avx_ptr)[SABER_K][SABER_N/16];

		
	  mask_ar[0]=~(0UL);mask_ar[1]=~(0UL);mask_ar[2]=~(0UL);mask_ar[3]=~(0UL);
	  mask_load = _mm256_loadu_si256 ((__m256i const *)mask_ar);

	  mod_p=_mm256_set1_epi16(SABER_P-1);

	  

	floor_round=_mm256_set1_epi16(4);

	H1_avx=_mm256_set1_epi16(h1);
 
	//--------------AVX declaration ends------------------
	load_values();
      
	for(i=0;i<SABER_SEEDBYTES;i++){ // Load the seedbytes in the client seed from PK.
		//seed[i]=pk[ SABER_POLYVECCOMPRESSEDBYTES + i];
		seed0[i]=pk0[ SABER_POLYVECCOMPRESSEDBYTES + i];
		seed1[i]=pk1[ SABER_POLYVECCOMPRESSEDBYTES + i];
		seed2[i]=pk2[ SABER_POLYVECCOMPRESSEDBYTES + i];
		seed3[i]=pk3[ SABER_POLYVECCOMPRESSEDBYTES + i]; 
	}

	GenMatrix4x(a0,a1,a2,a3, seed0, seed1, seed2, seed3);				

	GenSecret4x(skpv1_0,skpv1_1,skpv1_2,skpv1_3, noiseseed0, noiseseed1, noiseseed2, noiseseed3);

	// ----------- Load skpv1 into avx vectors ---------- 
	for(i=0;i<SABER_K;i++){ 
		for(j=0; j<SABER_N/16; j++){
		    //sk_avx[i][j] = _mm256_loadu_si256 ((__m256i const *) (&skpv1[i][j*16]));
		    sk0_avx[i][j] = _mm256_loadu_si256 ((__m256i const *) (&skpv1_0[i][j*16]));
		    sk1_avx[i][j] = _mm256_loadu_si256 ((__m256i const *) (&skpv1_1[i][j*16]));
		    sk2_avx[i][j] = _mm256_loadu_si256 ((__m256i const *) (&skpv1_2[i][j*16]));
		    sk3_avx[i][j] = _mm256_loadu_si256 ((__m256i const *) (&skpv1_3[i][j*16]));
		}
  	}

	// ----------- Load skpv1 into avx vectors ---------- 
	  for(i=0;i<SABER_K;i++){ 
		  for(j=0;j<SABER_K;j++){
			  for(k=0;k<SABER_N/16;k++){
				//a_avx[i][j][k]=_mm256_loadu_si256 ((__m256i const *) (&a[i].vec[j].coeffs[k*16]));
				a0_avx[i][j][k]=_mm256_loadu_si256 ((__m256i const *) (&a0[i].vec[j].coeffs[k*16]));
				a1_avx[i][j][k]=_mm256_loadu_si256 ((__m256i const *) (&a1[i].vec[j].coeffs[k*16]));
				a2_avx[i][j][k]=_mm256_loadu_si256 ((__m256i const *) (&a2[i].vec[j].coeffs[k*16]));
				a3_avx[i][j][k]=_mm256_loadu_si256 ((__m256i const *) (&a3[i].vec[j].coeffs[k*16]));
			  }
		  }
 	 }


	int SABER_SERIAL=0;
	
	//////////////////////////////////////////////////////////////////////////////
	// Serially compute (ciphertext0,ss0), (ciphertext1,ss1), (ciphertext2,ss2), (ciphertext3,ss3)

	for(SABER_SERIAL=0; SABER_SERIAL<4; SABER_SERIAL++)
	{
		if(SABER_SERIAL==0) {
			sk_avx_ptr = sk0_avx; a_avx_ptr=a0_avx; 
			pk_ptr=pk0;  
			ciphertext_ptr=ciphertext0; message_received_ptr=message_received0;
		}
		if(SABER_SERIAL==1) {
			sk_avx_ptr = sk1_avx; a_avx_ptr=a1_avx;  
			pk_ptr=pk1;  
			ciphertext_ptr=ciphertext1; message_received_ptr=message_received1;
		}
		if(SABER_SERIAL==2) {
			sk_avx_ptr = sk2_avx; a_avx_ptr=a2_avx;  
			pk_ptr=pk2; 
			ciphertext_ptr=ciphertext2; message_received_ptr=message_received2;
		}
		if(SABER_SERIAL==3) {
			sk_avx_ptr = sk3_avx; a_avx_ptr=a3_avx;  
			pk_ptr=pk3;  
			ciphertext_ptr=ciphertext3; message_received_ptr=message_received3;
		}

		for(i=0;i<SABER_K;i++){
			for(j=0;j<SABER_N/16;j++){
				res_avx[i][j]=_mm256_xor_si256(res_avx[i][j],res_avx[i][j]);
			}
		}	

		// Matrix-vector multiplication; 
		for(i=0;i<SABER_K;i++){
			for(j=0;j<SABER_K;j++){

				toom_cook_4way_avx(a_avx_ptr[i][j], sk_avx_ptr[j], SABER_Q, acc);

				for(k=0;k<SABER_N/16;k++){
					res_avx[i][k]=_mm256_add_epi16(res_avx[i][k],acc[k]);
				}
			
			}
		}
	
		// Now truncation

		for(i=0;i<SABER_K;i++){ //shift right EQ-EP bits
			for(j=0;j<SABER_N/16;j++){
				res_avx[i][j]=_mm256_add_epi16 (res_avx[i][j], H1_avx);
				res_avx[i][j]=_mm256_srli_epi16 (res_avx[i][j], (SABER_EQ-SABER_EP) );
			}
		}

		//-----this result should be put in b_prime for later use in server.
		for(i=0;i<SABER_K;i++){ // first store in 16 bit arrays
				for(j=0;j<SABER_N/16;j++){
				_mm256_maskstore_epi32 ((int *)(temp[i]+j*16), mask_load, res_avx[i][j]);
				}
			}

		POLVEC2BS(ciphertext_ptr, temp, SABER_P); // Pack b_prime into ciphertext byte string

		// client matrix-vector multiplication ends

		//------now calculate the v'

		//-------unpack the public_key
		BS2POLVEC(pk_ptr, pkcl, SABER_P);

		for(i=0;i<SABER_K;i++){
			for(j=0; j<SABER_N/16; j++){
				  pkcl_avx[i][j] = _mm256_loadu_si256 ((__m256i const *) (&pkcl[i][j*16]));
			}
		}

		// InnerProduct
		for(k=0;k<SABER_N/16;k++){
			vprime_avx[k]=_mm256_xor_si256(vprime_avx[k],vprime_avx[k]);
		}

		// vector-vector scalar multiplication with mod p
		clock1=cpucycles();
		count_mul++;
		for(j=0;j<SABER_K;j++){
			toom_cook_4way_avx(pkcl_avx[j], sk_avx_ptr[j], SABER_P, acc);

				for(k=0;k<SABER_N/16;k++){
					vprime_avx[k]=_mm256_add_epi16(vprime_avx[k],acc[k]);
				}
		}
		clock2=cpucycles();
		clock_mul=clock_mul+(clock2-clock1);

		// Computation of v'+h1 
		for(i=0;i<SABER_N/16;i++){//adding h1
	 		vprime_avx[i]=_mm256_add_epi16(vprime_avx[i], H1_avx);
		}

		// unpack message_received;
		for(j=0; j<SABER_KEYBYTES; j++)
		{
			for(i=0; i<8; i++)
			{
				message[8*j+i] = ((message_received_ptr[j]>>i) & 0x01);
			}
		}
		// message encoding
		for(i=0; i<SABER_N/16; i++)
		{
			message_avx[i] = _mm256_loadu_si256 ((__m256i const *) (&message[i*16]));
			message_avx[i] = _mm256_slli_epi16 (message_avx[i], (SABER_EP-1) );
		}	

		// SHIFTRIGHT(v'+h1-m mod p, EP-ET)
		for(k=0;k<SABER_N/16;k++)
		{
			vprime_avx[k]=_mm256_sub_epi16(vprime_avx[k], message_avx[k]);
			vprime_avx[k]=_mm256_and_si256(vprime_avx[k], mod_p);
			vprime_avx[k]=_mm256_srli_epi16 (vprime_avx[k], (SABER_EP-SABER_ET) );
		}

		// Unpack avx
		for(j=0;j<SABER_N/16;j++)
		{
				_mm256_maskstore_epi32 ((int *) (temp[0]+j*16), mask_load, vprime_avx[j]);
		}
		#if Saber_type == 1
			SABER_pack_3bit(&ciphertext_ptr[SABER_CIPHERTEXTBYTES], temp[0]);
		#elif Saber_type == 2
			SABER_pack_4bit(&ciphertext_ptr[SABER_CIPHERTEXTBYTES], temp[0]);
		#elif Saber_type == 3
			SABER_pack_6bit(&ciphertext_ptr[SABER_CIPHERTEXTBYTES], temp[0]);
		#endif


	}
	// End: (ciphertext0,ss0), (ciphertext1,ss1), (ciphertext2,ss2), (ciphertext3,ss3)
	//////////////////////////////////////////////////////////////////////////////////////

}


void indcpa_kem_dec(
					const unsigned char *sk0, const unsigned char *ciphertext0, unsigned char message0_dec[],
					const unsigned char *sk1, const unsigned char *ciphertext1, unsigned char message1_dec[],
					const unsigned char *sk2, const unsigned char *ciphertext2, unsigned char message2_dec[],
					const unsigned char *sk3, const unsigned char *ciphertext3, unsigned char message3_dec[])
{

	uint32_t i,j,k;
	
	
	uint16_t sksv[SABER_K][SABER_N]; //secret key of the server

	uint16_t pksv[SABER_K][SABER_N];
	
	uint16_t message_dec_unpacked[SABER_KEYBYTES*8];	// one element containes on decrypted bit;

	uint16_t op[SABER_N];

	const unsigned char *sk_ptr, *ciphertext_ptr;
	unsigned char *message_dec_ptr;

	//--------------AVX declaration------------------
	
	__m256i v_avx[SABER_N/16];
	__m256i acc[2*SABER_N/16];
	__m256i sksv_avx[SABER_K][SABER_N/16];
	__m256i pksv_avx[SABER_K][SABER_N/16];
	  
	mask_ar[0]=~(0UL);mask_ar[1]=~(0UL);mask_ar[2]=~(0UL);mask_ar[3]=~(0UL);
	mask_load = _mm256_loadu_si256 ((__m256i const *)mask_ar);


	H2_avx=_mm256_set1_epi16(h2);

	//--------------AVX declaration ends------------------
	
 	load_values();

	// Compute serially;


	// Start: Decrypt ciphertext
	////////////////////////////////////////////////////////////////////////////

	int SABER_SERIAL=0;
	
	//////////////////////////////////////////////////////////////////////////////
	// Serially decrypt (ciphertext0, ciphertext1 ...)

	for(SABER_SERIAL=0; SABER_SERIAL<4; SABER_SERIAL++)
	{
		if(SABER_SERIAL==0) { sk_ptr=sk0; ciphertext_ptr=ciphertext0; message_dec_ptr=message0_dec; }
		if(SABER_SERIAL==1) { sk_ptr=sk1; ciphertext_ptr=ciphertext1; message_dec_ptr=message1_dec; }
		if(SABER_SERIAL==2) { sk_ptr=sk2; ciphertext_ptr=ciphertext2; message_dec_ptr=message2_dec; }
		if(SABER_SERIAL==3) { sk_ptr=sk3; ciphertext_ptr=ciphertext3; message_dec_ptr=message3_dec; }
		//if(SABER_SERIAL==4) { sk_ptr=sk; ciphertext_ptr=ciphertext; message_dec_ptr=message_dec; }
	
		//-------unpack the public_key
		BS2POLVEC(sk_ptr, sksv, SABER_Q); //sksv is the secret-key
		BS2POLVEC(ciphertext_ptr, pksv, SABER_P); //pksv is the ciphertext

		for(i=0;i<SABER_K;i++){
			for(j=0; j<SABER_N/16; j++){
				  sksv_avx[i][j] = _mm256_loadu_si256 ((__m256i const *) (&sksv[i][j*16]));
				  pksv_avx[i][j] = _mm256_loadu_si256 ((__m256i const *) (&pksv[i][j*16]));
			}
		}

		for(i=0;i<SABER_N/16;i++){
			v_avx[i]=_mm256_xor_si256(v_avx[i],v_avx[i]);
		}

		// InnerProduct(b', s, mod p)
		for(j=0;j<SABER_K;j++){

			toom_cook_4way_avx(pksv_avx[j], sksv_avx[j], SABER_P, acc);

				for(k=0;k<SABER_N/16;k++){
					v_avx[k]=_mm256_add_epi16(v_avx[k],acc[k]);
				}
		}

		for(i=0; i<SABER_N/16; i++){
			_mm256_maskstore_epi32 ((int *)(message_dec_unpacked+i*16), mask_load, v_avx[i]);
		}


		#if Saber_type == 1
			SABER_un_pack3bit(&ciphertext_ptr[SABER_CIPHERTEXTBYTES], op);
		#elif Saber_type == 2
			SABER_un_pack4bit(&ciphertext_ptr[SABER_CIPHERTEXTBYTES], op);
		#elif Saber_type == 3
			SABER_un_pack6bit(&ciphertext_ptr[SABER_CIPHERTEXTBYTES], op);
		#endif

		//addition of h2
		for(i=0;i<SABER_N;i++){
			message_dec_unpacked[i]= ( ( message_dec_unpacked[i] + h2 - (op[i]<<(SABER_EP-SABER_ET)) ) & (SABER_P-1) ) >> (SABER_EP-1);
		}

		POL2MSG(message_dec_unpacked, message_dec_ptr);
	}
	// End: Decrypt ciphertext
	////////////////////////////////////////////////////////////////////////////

}

void POL2MSG(uint16_t *message_dec_unpacked, unsigned char *message_dec){

	int32_t i,j;

	for(j=0; j<SABER_KEYBYTES; j++)
	{
		message_dec[j] = 0;
		for(i=0; i<8; i++)
		message_dec[j] = message_dec[j] | (message_dec_unpacked[j*8 + i] <<i);
	} 

}


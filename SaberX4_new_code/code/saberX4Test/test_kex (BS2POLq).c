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
uint64_t clock_S0,clock_S1,clock_S2;
uint64_t clock_E0,clock_E1,clock_E2;
uint64_t clock_S10,clock_S11,clock_S12;
uint64_t clock_E10,clock_E11,clock_E12;

long long int clock_T;
long long int clock_T0,clock_T1,clock_T2;
long long int clock_T10,clock_T11,clock_T12;

unsigned char *bs0, *bs1, *bs2, *bs3;
__m256i data0[16], data1[16], data2[16], data3[16];
__m256i avx_epi64_01, avx_epi64_03, avx_epi64_07, avx_epi64_0f, avx_epi64_1f, avx_epi64_3f, avx_epi64_7f, avx_epi64_ff;
__m256i avx_epi32_01, avx_epi32_03, avx_epi32_07, avx_epi32_0f, avx_epi32_1f, avx_epi32_3f, avx_epi32_7f, avx_epi32_ff;
uint16_t dataN0[256], dataN1[256], dataN2[256], dataN3[256];
unsigned char* stream;


void init_avx(){
	avx_epi64_01 = _mm256_set1_epi64x(0x01);
	avx_epi64_03 = _mm256_set1_epi64x(0x03);
	avx_epi64_07 = _mm256_set1_epi64x(0x07);
	avx_epi64_0f = _mm256_set1_epi64x(0x0f);
	avx_epi64_1f = _mm256_set1_epi64x(0x1f);
	avx_epi64_3f = _mm256_set1_epi64x(0x3f);
	avx_epi64_7f = _mm256_set1_epi64x(0x7f);
	avx_epi64_ff = _mm256_set1_epi64x(0xff);

	avx_epi32_01 = _mm256_set1_epi32(0x01);
	avx_epi32_03 = _mm256_set1_epi32(0x03);
	avx_epi32_07 = _mm256_set1_epi32(0x07);
	avx_epi32_0f = _mm256_set1_epi32(0x0f);
	avx_epi32_1f = _mm256_set1_epi32(0x1f);
	avx_epi32_3f = _mm256_set1_epi32(0x3f);
	avx_epi32_7f = _mm256_set1_epi32(0x7f);
	avx_epi32_ff = _mm256_set1_epi32(0xff);
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

void test_BS2POLq(const unsigned char *bytes, uint16_t data[SABER_N]){
	
	clock_S = cpucycles();

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

	clock_E = cpucycles();
	clock_T = clock_T+clock_E-clock_S;

}


void test_bs2polqX4(unsigned char *bytes0, unsigned char *bytes1, unsigned char *bytes2, unsigned char *bytes3, 
				__m256i data0[SABER_N/16], __m256i data1[SABER_N/16], __m256i data2[SABER_N/16], __m256i data3[SABER_N/16]){
	//clock_inBS2POLq = cpucycles();

	clock_S0 = cpucycles();
	__m256i avx_bytes[32*13];
	__m256i avx_data[256];
	//uint64_t tmp_a[32*13][4];
	int i = 0;
	
	//clock_before_temp_loop = cpucycles();
	
	for(i=0;i<32*13;++i){
		avx_bytes[i] = _mm256_set_epi64x(bytes3[i], bytes2[i], bytes1[i], bytes0[i]);
	}
	
	clock_E0 = cpucycles();
	clock_T0 = clock_T0+clock_E0-clock_S0;

	//clock_after_temp_loop = cpucycles();
	//printf("BS2POLq loop CPUcycle: %d\n",clock_after_temp_loop-clock_before_temp_loop);
	
	//clock_before_shift = cpucycles();
	clock_S1 = cpucycles();
	uint32_t j;
	uint32_t offset_data=0,offset_byte=0;
	offset_byte=0;
	for(j=0;j<SABER_N/8;j++){
		offset_byte=13*j;
		offset_data=8*j;
		avx_data[offset_data + 0] = _mm256_or_si256(_mm256_and_si256(avx_bytes[offset_byte+0],avx_epi64_ff),
									_mm256_slli_epi64(_mm256_and_si256(avx_bytes[offset_byte+1],avx_epi64_1f),8));
		avx_data[offset_data + 1] = _mm256_or_si256(_mm256_srli_epi64(avx_bytes[offset_byte+1],5),
									_mm256_or_si256(_mm256_slli_epi64(avx_bytes[offset_byte+2],3),
									_mm256_slli_epi64(_mm256_and_si256(avx_bytes[offset_byte+3],avx_epi64_03),11)));
		avx_data[offset_data + 2] = _mm256_or_si256(_mm256_srli_epi64(avx_bytes[offset_byte+3],2),
									_mm256_slli_epi64(_mm256_and_si256(avx_bytes[offset_byte+4],avx_epi64_7f),6));
		avx_data[offset_data + 3] = _mm256_or_si256(_mm256_and_si256(_mm256_srli_epi64(avx_bytes[offset_byte+4],7),avx_epi64_01),
									_mm256_or_si256(_mm256_slli_epi64(avx_bytes[offset_byte+5],1),
									_mm256_slli_epi64(_mm256_and_si256(avx_bytes[offset_byte+6],avx_epi64_0f),9)));
		avx_data[offset_data + 4] = _mm256_or_si256(_mm256_and_si256(_mm256_srli_epi64(avx_bytes[offset_byte+6],4),avx_epi64_0f),
									_mm256_or_si256(_mm256_slli_epi64(_mm256_and_si256(avx_bytes[offset_byte+7],avx_epi64_ff),4),
									_mm256_slli_epi64(_mm256_and_si256(avx_bytes[offset_byte+8],avx_epi64_01),12)));
		avx_data[offset_data + 5] = _mm256_or_si256(_mm256_and_si256(_mm256_srli_epi64(avx_bytes[offset_byte+8],1),avx_epi64_7f),
									_mm256_slli_epi64(_mm256_and_si256(avx_bytes[offset_byte+9],avx_epi64_3f),7));
		avx_data[offset_data + 6] = _mm256_or_si256(_mm256_and_si256(_mm256_srli_epi64(avx_bytes[offset_byte+9],6),avx_epi64_03),
									_mm256_or_si256(_mm256_slli_epi64(_mm256_and_si256(avx_bytes[offset_byte+10],avx_epi64_ff),2),
									_mm256_slli_epi64(_mm256_and_si256(avx_bytes[offset_byte+11],avx_epi64_07),10)));
		avx_data[offset_data + 7] = _mm256_or_si256(_mm256_and_si256(_mm256_srli_epi64(avx_bytes[offset_byte+11],3),avx_epi64_1f),
									_mm256_slli_epi64(_mm256_and_si256(avx_bytes[offset_byte+12],avx_epi64_ff),5));
	}

	clock_E1 = cpucycles();
	clock_T1 = clock_T1+clock_E1-clock_S1;
	

	//clock_after_shift = cpucycles();
	//printf("BS2POLq shift CPUcycle: %d\n",clock_after_shift-clock_before_shift);
	clock_S2 = cpucycles();
	//clock_before_store = cpucycles();
	uint32_t k;
	uint32_t offset_extract;
	//int flag = 1;
	// SHOW(uint64_t, avx_data[0][0],"bs2polq");
	// SHOW(uint64_t, avx_data[0][1],"bs2polq");
	// SHOW(uint64_t, avx_data[0][2],"bs2polq");
	// SHOW(uint64_t, avx_data[0][3],"bs2polq");
	

	
	for(k=0;k<SABER_N/16;++k){
		offset_extract = 16*k;
		
		
		data0[k] = _mm256_set_epi16(avx_data[offset_extract + 15][0], avx_data[offset_extract + 14][0], avx_data[offset_extract + 13][0], avx_data[offset_extract + 12][0],
						avx_data[offset_extract + 11][0], avx_data[offset_extract + 10][0], avx_data[offset_extract + 9][0], avx_data[offset_extract + 8][0],
						avx_data[offset_extract + 7][0], avx_data[offset_extract + 6][0], avx_data[offset_extract + 5][0], avx_data[offset_extract + 4][0],
						avx_data[offset_extract + 3][0], avx_data[offset_extract + 2][0], avx_data[offset_extract + 1][0], avx_data[offset_extract + 0][0]);				
		
		data1[k] = _mm256_set_epi16(avx_data[offset_extract + 15][1], avx_data[offset_extract + 14][1], avx_data[offset_extract + 13][1], avx_data[offset_extract + 12][1],
						avx_data[offset_extract + 11][1], avx_data[offset_extract + 10][1], avx_data[offset_extract + 9][1], avx_data[offset_extract + 8][1],
						avx_data[offset_extract + 7][1], avx_data[offset_extract + 6][1], avx_data[offset_extract + 5][1], avx_data[offset_extract + 4][1],
						avx_data[offset_extract + 3][1], avx_data[offset_extract + 2][1], avx_data[offset_extract + 1][1], avx_data[offset_extract + 0][1]);				
		
		data2[k] =_mm256_set_epi16(avx_data[offset_extract + 15][2], avx_data[offset_extract + 14][2], avx_data[offset_extract + 13][2], avx_data[offset_extract + 12][2],
						avx_data[offset_extract + 11][2], avx_data[offset_extract + 10][2], avx_data[offset_extract + 9][2], avx_data[offset_extract + 8][2],
						avx_data[offset_extract + 7][2], avx_data[offset_extract + 6][2], avx_data[offset_extract + 5][2], avx_data[offset_extract + 4][2],
						avx_data[offset_extract + 3][2], avx_data[offset_extract + 2][2], avx_data[offset_extract + 1][2], avx_data[offset_extract + 0][2]);				
		
		data3[k] = _mm256_set_epi16(avx_data[offset_extract + 15][3], avx_data[offset_extract + 14][3], avx_data[offset_extract + 13][3], avx_data[offset_extract + 12][3],
						avx_data[offset_extract + 11][3], avx_data[offset_extract + 10][3], avx_data[offset_extract + 9][3], avx_data[offset_extract + 8][3],
						avx_data[offset_extract + 7][3], avx_data[offset_extract + 6][3], avx_data[offset_extract + 5][3], avx_data[offset_extract + 4][3],
						avx_data[offset_extract + 3][3], avx_data[offset_extract + 2][3], avx_data[offset_extract + 1][3], avx_data[offset_extract + 0][3]);				
		
	}
		
		clock_E2 = cpucycles();
		clock_T2 = clock_T2+clock_E2-clock_S2;

}


void test_improved_BS2POLq(unsigned char * bs0, unsigned char * bs1, unsigned char * bs2, unsigned char * bs3,
				__m256i data0[SABER_N/16], __m256i data1[SABER_N/16], __m256i data2[SABER_N/16], __m256i data3[SABER_N/16]){
	
	/*--------------------------------------------------------*/
	/*------------------Load to avx_bytes---------------------*/
	/*--------------------------------------------------------*/
	// int flag=1;
	// int flag0=1;
	// int flag1=1;
	clock_S10 = cpucycles();

	const float * cfbs0 = (const float *)bs0;
	const float * cfbs1 = (const float *)bs1;
	const float * cfbs2 = (const float *)bs2;
	const float * cfbs3 = (const float *)bs3;

	const float * cfbs4, * cfbs5, * cfbs6, * cfbs7;
	
	__m128 l0,l1,l2,l3,l4,l5,l6,l7;
	__m128 l8,l9,l10,l11,l12,l13,l14,l15;
	__m256i l0_256,l1_256,l2_256,l3_256;
	__m256i avx_bytes[32*13/2];
	__m256i avx_data[256/2];
	
	
	// if(flag0){
		// SHOW(float,bs0[0]);
		// flag0=0;
	// }
	
	const float *bs4, *bs5, *bs6, *bs7;

	cfbs4 = cfbs0+13*8/2;
	cfbs5 = cfbs1+13*8/2;
	cfbs6 = cfbs2+13*8/2;
	cfbs7 = cfbs3+13*8/2;
	int p=0;
	int i;
	for(i = 0;i<13;i++){
		l0 = _mm_load_ps(cfbs0+p*4); //0
		l1 = _mm_load_ps(cfbs1+p*4); //160
		l2 = _mm_load_ps(cfbs2+p*4); //64
		l3 = _mm_load_ps(cfbs3+p*4); //224

		l4 = _mm_load_ps(cfbs4+p*4); //208
		l5 = _mm_load_ps(cfbs5+p*4); //112
		l6 = _mm_load_ps(cfbs6+p*4); //16
		l7 = _mm_load_ps(cfbs7+p*4); //176

		_MM_TRANSPOSE4_PS(l0,l1,l2,l3);
		_MM_TRANSPOSE4_PS(l4,l5,l6,l7);

		l0_256 = _mm256_insertf128_si256(_mm256_castsi128_si256((__m128i)l0),(__m128i)l4,1);
		l1_256 = _mm256_insertf128_si256(_mm256_castsi128_si256((__m128i)l1),(__m128i)l5,1);
		l2_256 = _mm256_insertf128_si256(_mm256_castsi128_si256((__m128i)l2),(__m128i)l6,1);
		l3_256 = _mm256_insertf128_si256(_mm256_castsi128_si256((__m128i)l3),(__m128i)l7,1);

		avx_bytes[i*16+3] = _mm256_srli_epi32(l0_256,24);
		avx_bytes[i*16+2] = _mm256_srli_epi32(_mm256_slli_epi32(l0_256,8),24);
		avx_bytes[i*16+1] = _mm256_srli_epi32(_mm256_slli_epi32(l0_256,16),24);
		avx_bytes[i*16+0] = _mm256_srli_epi32(_mm256_slli_epi32(l0_256,24),24);

		avx_bytes[i*16+7] = _mm256_srli_epi32(l1_256,24);
		avx_bytes[i*16+6] = _mm256_srli_epi32(_mm256_slli_epi32(l1_256,8),24);
		avx_bytes[i*16+5] = _mm256_srli_epi32(_mm256_slli_epi32(l1_256,16),24);
		avx_bytes[i*16+4] = _mm256_srli_epi32(_mm256_slli_epi32(l1_256,24),24);

		avx_bytes[i*16+11] = _mm256_srli_epi32(l2_256,24);
		avx_bytes[i*16+10] = _mm256_srli_epi32(_mm256_slli_epi32(l2_256,8),24);
		avx_bytes[i*16+9] = _mm256_srli_epi32(_mm256_slli_epi32(l2_256,16),24);
		avx_bytes[i*16+8] = _mm256_srli_epi32(_mm256_slli_epi32(l2_256,24),24);

		avx_bytes[i*16+15] = _mm256_srli_epi32(l3_256,24);
		avx_bytes[i*16+14] = _mm256_srli_epi32(_mm256_slli_epi32(l3_256,8),24);
		avx_bytes[i*16+13] = _mm256_srli_epi32(_mm256_slli_epi32(l3_256,16),24);
		avx_bytes[i*16+12] = _mm256_srli_epi32(_mm256_slli_epi32(l3_256,24),24);
		p++;
	}
	//printf("Print?");
	//SHOW(__m256i,avx_bytes[0]);
	//SHOW(__m256i,avx_bytes[(32*13/2)-1]);
	//printf("%lld\n",avx_bytes[0][0]);
	clock_E10 = cpucycles();
	clock_T10 += clock_E10 - clock_S10;
	
	// SHOW(__m256i,avx_bytes[0]);
	// SHOW(__m256i,avx_bytes[1]);

	/*--------------------------------------------------------*/
	/*---------------------Do the shifts----------------------*/
	/*--------------------------------------------------------*/
	clock_S11 = cpucycles();
	
	// SHOW(__m256i,avx_bytes[0]);
	// SHOW(__m256i,_mm256_and_si256(avx_bytes[0],avx_epi32_ff));
	// SHOW(__m256i,_mm256_and_si256(avx_bytes[1],avx_epi32_1f));
	// SHOW(__m256i,_mm256_slli_epi32(_mm256_and_si256(avx_bytes[1],avx_epi32_1f),8));
	
	uint32_t j;
	uint32_t offset_data=0,offset_byte=0;
	offset_byte=0;
	for(j=0;j<(SABER_N/8)/2;j++){
		offset_byte=13*j;
		offset_data=8*j;
		avx_data[offset_data + 0] = _mm256_or_si256(_mm256_and_si256(avx_bytes[offset_byte+0],avx_epi32_ff),
									_mm256_slli_epi32(_mm256_and_si256(avx_bytes[offset_byte+1],avx_epi32_1f),8));
		
		avx_data[offset_data + 1] = _mm256_or_si256(_mm256_srli_epi32(avx_bytes[offset_byte+1],5),
									_mm256_or_si256(_mm256_slli_epi32(avx_bytes[offset_byte+2],3),
									_mm256_slli_epi32(_mm256_and_si256(avx_bytes[offset_byte+3],avx_epi32_03),11)));
		avx_data[offset_data + 2] = _mm256_or_si256(_mm256_srli_epi32(avx_bytes[offset_byte+3],2),
									_mm256_slli_epi32(_mm256_and_si256(avx_bytes[offset_byte+4],avx_epi32_7f),6));
		avx_data[offset_data + 3] = _mm256_or_si256(_mm256_and_si256(_mm256_srli_epi32(avx_bytes[offset_byte+4],7),avx_epi32_01),
									_mm256_or_si256(_mm256_slli_epi32(avx_bytes[offset_byte+5],1),
									_mm256_slli_epi32(_mm256_and_si256(avx_bytes[offset_byte+6],avx_epi32_0f),9)));
		avx_data[offset_data + 4] = _mm256_or_si256(_mm256_and_si256(_mm256_srli_epi32(avx_bytes[offset_byte+6],4),avx_epi32_0f),
									_mm256_or_si256(_mm256_slli_epi32(_mm256_and_si256(avx_bytes[offset_byte+7],avx_epi32_ff),4),
									_mm256_slli_epi32(_mm256_and_si256(avx_bytes[offset_byte+8],avx_epi32_01),12)));
		avx_data[offset_data + 5] = _mm256_or_si256(_mm256_and_si256(_mm256_srli_epi32(avx_bytes[offset_byte+8],1),avx_epi32_7f),
									_mm256_slli_epi32(_mm256_and_si256(avx_bytes[offset_byte+9],avx_epi32_3f),7));
		avx_data[offset_data + 6] = _mm256_or_si256(_mm256_and_si256(_mm256_srli_epi32(avx_bytes[offset_byte+9],6),avx_epi32_03),
									_mm256_or_si256(_mm256_slli_epi32(_mm256_and_si256(avx_bytes[offset_byte+10],avx_epi32_ff),2),
									_mm256_slli_epi32(_mm256_and_si256(avx_bytes[offset_byte+11],avx_epi32_07),10)));
		avx_data[offset_data + 7] = _mm256_or_si256(_mm256_and_si256(_mm256_srli_epi32(avx_bytes[offset_byte+11],3),avx_epi32_1f),
									_mm256_slli_epi32(_mm256_and_si256(avx_bytes[offset_byte+12],avx_epi32_ff),5));
	}
	clock_E11 = cpucycles();
	clock_T11 += clock_E11 - clock_S11;
	
	// SHOW(__m256i,avx_data[0]);
	// SHOW(__m256i,avx_data[1]);
	// SHOW(__m256i,avx_data[2]);
	// SHOW(__m256i,avx_data[3]);
	
	/*--------------------------------------------------------*/
	/*---------------------Store back-------------------------*/
	/*--------------------------------------------------------*/
	// _mm256_cvtps_ph
	// _mm_cvtps_ph(l0,0x08)
	// _mm_cvtps_pi16
	__m128i ph0,ph1,ph2,ph3;
	__m128i ph4,ph5,ph6,ph7;
	//__m128i l0i,l1/i,l2i,l3i;
	//__m64 ph0,ph1,ph2,ph3;
	
	clock_S12 = cpucycles();
	
	int offset0 = 0;
	int k;
	for(k=0;k<8;k++){
		l0 = _mm256_extractf128_ps((__m256)avx_data[offset0+0],0);
		l1 = _mm256_extractf128_ps((__m256)avx_data[offset0+1],0);
		l2 = _mm256_extractf128_ps((__m256)avx_data[offset0+2],0);
		l3 = _mm256_extractf128_ps((__m256)avx_data[offset0+3],0);
		
		l4 = _mm256_extractf128_ps((__m256)avx_data[offset0+4],0);
		l5 = _mm256_extractf128_ps((__m256)avx_data[offset0+5],0);
		l6 = _mm256_extractf128_ps((__m256)avx_data[offset0+6],0);
		l7 = _mm256_extractf128_ps((__m256)avx_data[offset0+7],0);
		
		l8 = _mm256_extractf128_ps((__m256)avx_data[offset0+8],0);
		l9 = _mm256_extractf128_ps((__m256)avx_data[offset0+9],0);
		l10 = _mm256_extractf128_ps((__m256)avx_data[offset0+10],0);
		l11 = _mm256_extractf128_ps((__m256)avx_data[offset0+11],0);
		
		l12 = _mm256_extractf128_ps((__m256)avx_data[offset0+12],0);
		l13 = _mm256_extractf128_ps((__m256)avx_data[offset0+13],0);
		l14 = _mm256_extractf128_ps((__m256)avx_data[offset0+14],0);
		l15 = _mm256_extractf128_ps((__m256)avx_data[offset0+15],0);

		_MM_TRANSPOSE4_PS(l0,l1,l2,l3);
		_MM_TRANSPOSE4_PS(l4,l5,l6,l7);
		_MM_TRANSPOSE4_PS(l8,l9,l10,l11);
		_MM_TRANSPOSE4_PS(l12,l13,l14,l15);
		
		ph0 = _mm_packs_epi32((__m128i)l0,(__m128i)l4);
		ph1 = _mm_packs_epi32((__m128i)l1,(__m128i)l5);
		ph2 = _mm_packs_epi32((__m128i)l2,(__m128i)l6);
		ph3 = _mm_packs_epi32((__m128i)l3,(__m128i)l7);
		ph4 = _mm_packs_epi32((__m128i)l8,(__m128i)l12);
		ph5 = _mm_packs_epi32((__m128i)l9,(__m128i)l13);
		ph6 = _mm_packs_epi32((__m128i)l10,(__m128i)l14);
		ph7 = _mm_packs_epi32((__m128i)l11,(__m128i)l15);
		
		//SHOW(__m256i,_mm256_insertf128_si256(_mm256_castsi128_si256(ph0),ph4,1));
		data0[k]   = _mm256_insertf128_si256(_mm256_castsi128_si256(ph0),ph4,1);
		data1[k] = _mm256_insertf128_si256(_mm256_castsi128_si256(ph1),ph5,1);
		data2[k] = _mm256_insertf128_si256(_mm256_castsi128_si256(ph2),ph6,1);
		data3[k] = _mm256_insertf128_si256(_mm256_castsi128_si256(ph3),ph7,1);
		
		l0 = _mm256_extractf128_ps((__m256)avx_data[offset0+0],1);
		l1 = _mm256_extractf128_ps((__m256)avx_data[offset0+1],1);
		l2 = _mm256_extractf128_ps((__m256)avx_data[offset0+2],1);
		l3 = _mm256_extractf128_ps((__m256)avx_data[offset0+3],1);
		
		l4 = _mm256_extractf128_ps((__m256)avx_data[offset0+4],1);
		l5 = _mm256_extractf128_ps((__m256)avx_data[offset0+5],1);
		l6 = _mm256_extractf128_ps((__m256)avx_data[offset0+6],1);
		l7 = _mm256_extractf128_ps((__m256)avx_data[offset0+7],1);
		
		l8 = _mm256_extractf128_ps((__m256)avx_data[offset0+8],1);
		l9 = _mm256_extractf128_ps((__m256)avx_data[offset0+9],1);
		l10 = _mm256_extractf128_ps((__m256)avx_data[offset0+10],1);
		l11 = _mm256_extractf128_ps((__m256)avx_data[offset0+11],1);
		
		l12 = _mm256_extractf128_ps((__m256)avx_data[offset0+12],1);
		l13 = _mm256_extractf128_ps((__m256)avx_data[offset0+13],1);
		l14 = _mm256_extractf128_ps((__m256)avx_data[offset0+14],1);
		l15 = _mm256_extractf128_ps((__m256)avx_data[offset0+15],1);

		_MM_TRANSPOSE4_PS(l0,l1,l2,l3);
		_MM_TRANSPOSE4_PS(l4,l5,l6,l7);
		_MM_TRANSPOSE4_PS(l8,l9,l10,l11);
		_MM_TRANSPOSE4_PS(l12,l13,l14,l15);
		
		ph0 = _mm_packs_epi32((__m128i)l0,(__m128i)l4);
		ph1 = _mm_packs_epi32((__m128i)l1,(__m128i)l5);
		ph2 = _mm_packs_epi32((__m128i)l2,(__m128i)l6);
		ph3 = _mm_packs_epi32((__m128i)l3,(__m128i)l7);
		ph4 = _mm_packs_epi32((__m128i)l8,(__m128i)l12);
		ph5 = _mm_packs_epi32((__m128i)l9,(__m128i)l13);
		ph6 = _mm_packs_epi32((__m128i)l10,(__m128i)l14);
		ph7 = _mm_packs_epi32((__m128i)l11,(__m128i)l15);

		
		data0[k+8] = _mm256_insertf128_si256(_mm256_castsi128_si256(ph0),ph4,1);
		data1[k+8] = _mm256_insertf128_si256(_mm256_castsi128_si256(ph1),ph5,1);
		data2[k+8] = _mm256_insertf128_si256(_mm256_castsi128_si256(ph2),ph6,1);
		data3[k+8] = _mm256_insertf128_si256(_mm256_castsi128_si256(ph3),ph7,1);
		
		offset0 += 16;
	}
	
	
	//SHOW(__m64,ph0);

	clock_E12 = cpucycles();
	clock_T12 += clock_E12 - clock_S12;
}


int main()
{
	
	//test_kem_cca();
	init_avx();
	gene_4bs();
	int i;
	//const float * cfbs0 = (const float *)bs0;
	//const float * cfbs1 = (const float *)bs1;
	//const float * cfbs2 = (const float *)bs2;
	//const float * cfbs3 = (const float *)bs3;
	for(i=0;i<100000;i++){
		test_BS2POLq(bs0,dataN0);
		test_BS2POLq(bs1,dataN1);
		test_BS2POLq(bs2,dataN2);
		test_BS2POLq(bs3,dataN3);
		test_bs2polqX4(bs0,bs1,bs2,bs3,data0,data1,data2,data3);
		test_improved_BS2POLq(bs0,bs1,bs2,bs3,data0,data1,data2,data3);
		
	}
	printf("Average times store(X1): \t %lu \n",clock_T/100000);

	printf("Average times load(X4): \t %lu \n",clock_T0/100000);
	printf("Average times shift(X4): \t %lu \n",clock_T1/100000);
	printf("Average times store(X4): \t %lu \n",clock_T2/100000);

	printf("Average times load(mo): \t %lu \n",clock_T10/100000);
	printf("Average times shift(mo): \t %lu \n",clock_T11/100000);
	printf("Average times store(mo): \t %lu \n",clock_T12/100000);
	
	return 0;
}

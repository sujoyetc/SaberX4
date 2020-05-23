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

uint32_t SABER_N = 256;
uint32_t SABER_K = 3;

uint64_t clock_S, clock_E;
uint64_t clock_S0, clock_E0,clock_S1, clock_E1,clock_S2, clock_E2;
long long int clock_T0,clock_T1,clock_T2, clock_T;
unsigned char *bs0, *bs1, *bs2, *bs3;
uint8_t bytes[128];
uint8_t bytes0[10000], bytes1[10000], bytes2[10000], bytes3[10000];

uint16_t data0[SABER_K][SABER_N], data1[SABER_K][SABER_N], data2[SABER_K][SABER_N], data3[SABER_K][SABER_N];

__m256i data0[16], data1[16], data2[16], data3[16];
__m256i avx_epi64_01, avx_epi64_03, avx_epi64_07, avx_epi64_0f, avx_epi64_1f, avx_epi64_3f, avx_epi64_7f, avx_epi64_ff;
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
  int i,j;
  for(i=0;i<SABER_K;i++){
    for(j=0;j<SABER_N;j++){
        data0[i][j] = rand ();
        data1[i][j] = rand ();
        data2[i][j] = rand ();
        data3[i][j] = rand ();
    }
  }
  //bs0 = gene_bytestream(256*16/8);
  //bs1 = gene_bytestream(256*16/8);
  //bs2 = gene_bytestream(256*16/8);
  //bs3 = gene_bytestream(256*16/8);
}


void test_pack_13bit(uint8_t *bytes, uint16_t data[SABER_K][SABER_N]){
	clock_S = cpucycles();
	uint32_t i,j;
	uint32_t offset_data=0,offset_byte=0,offset_byte1=0;	
	
	offset_byte=0;
	for(i=0;i<SABER_K;i++){
		offset_byte1=i*(SABER_N*13)/8;
		for(j=0;j<SABER_N/8;j++){
			offset_byte=offset_byte1+13*j;
			offset_data=8*j;
			bytes[offset_byte + 0]= ( data[i][ offset_data + 0 ] & (0xff));

			bytes[offset_byte + 1]= ( (data[i][ offset_data + 0 ] >>8) & 0x1f ) | ((data[i][ offset_data + 1 ] & 0x07) << 5);

			bytes[offset_byte + 2]= ( (data[i][ offset_data + 1 ] >>3) & 0xff );

			bytes[offset_byte + 3]= ( (data[i][ offset_data + 1 ] >>11) & 0x03 ) | ((data[i][ offset_data + 2 ] & 0x3f) << 2);

			bytes[offset_byte + 4]= ( (data[i][ offset_data + 2 ] >>6) & 0x7f ) | ( (data[i][ offset_data + 3 ] & 0x01) << 7 );

			bytes[offset_byte + 5]= ( (data[i][ offset_data + 3 ] >>1) & 0xff );

			bytes[offset_byte + 6]= ( (data[i][ offset_data + 3 ] >>9) & 0x0f ) | ( (data[i][ offset_data + 4 ] & 0x0f) << 4 );

			bytes[offset_byte + 7]= ( (data[i][ offset_data + 4] >>4) & 0xff );

			bytes[offset_byte + 8]= ( (data[i][ offset_data + 4 ] >>12) & 0x01 ) | ( (data[i][ offset_data + 5 ] & 0x7f) << 1 );

			bytes[offset_byte + 9]= ( (data[i][ offset_data + 5 ] >>7) & 0x3f ) | ( (data[i][ offset_data + 6 ] & 0x03) << 6 );

			bytes[offset_byte + 10]= ( (data[i][ offset_data + 6 ] >>2) & 0xff );

			bytes[offset_byte + 11]= ( (data[i][ offset_data + 6 ] >>10) & 0x07 ) | ( (data[i][ offset_data + 7 ] & 0x1f) << 3 );

			bytes[offset_byte + 12]= ( (data[i][ offset_data + 7 ] >>5) & 0xff );

		}
	}
	clock_E = cpucycles();
	clock_T += clock_E - clock_S;

}

void test_pack_13bitX4(uint8_t *bytes0, uint8_t *bytes1, uint8_t *bytes2, uint8_t *bytes3, uint16_t data0[SABER_K][SABER_N], uint16_t data1[SABER_K][SABER_N], uint16_t data2[SABER_K][SABER_N], uint16_t data3[SABER_K][SABER_N]){

	clock_S0 = cpucycles();
	uint32_t i,j,m;
	uint32_t offset_data=0,offset_byte=0,offset_byte1=0;	
	
	__m256i avx_data[SABER_K*SABER_N];
	for(i=0;i<SABER_K;++i){
		for(j=0;j<SABER_N;++j){
			avx_data[i*256+j] = _mm256_set_epi64x(data3[i][j], data2[i][j], data1[i][j], data0[i][j]);
		}
	}
	clock_E0 = cpucycles();
	clock_T0 += clock_E0 - clock_S0;

	clock_S1 = cpucycles();
	__m256i *avx_bytes0_ptr = (__m256i*)bytes0;
	__m256i *avx_bytes1_ptr = (__m256i*)bytes1;
	__m256i *avx_bytes2_ptr = (__m256i*)bytes2;
	__m256i *avx_bytes3_ptr = (__m256i*)bytes3;
	
	offset_byte=0;
	//para=1248
	uint16_t para = (SABER_K-1)*((SABER_N*13)/8)+13*(SABER_N/8);
	
	__m256i avx_bytes[para];
	
	for(i=0;i<SABER_K;i++){
		offset_byte1=i*(SABER_N*13)/8;
		for(j=0;j<SABER_N/8;j++){
			offset_byte=offset_byte1+13*j;
			offset_data=8*j;
			avx_bytes[offset_byte + 0 ]= _mm256_and_si256( avx_data[i*SABER_N + offset_data + 0 ],(avx_epi64_ff));

			avx_bytes[offset_byte + 1 ]= _mm256_or_si256(_mm256_and_si256( _mm256_srli_epi64(avx_data[i*SABER_N + offset_data + 0 ] ,8 ),avx_epi64_1f ) , _mm256_slli_epi64(_mm256_and_si256(avx_data[i*SABER_N + offset_data + 1 ],avx_epi64_07) , 5));

			avx_bytes[offset_byte + 2 ]= _mm256_and_si256( _mm256_srli_epi64(avx_data[i*SABER_N + offset_data + 1 ] ,3 ),avx_epi64_ff );

			avx_bytes[offset_byte + 3 ]= _mm256_or_si256(_mm256_and_si256( _mm256_srli_epi64(avx_data[i*SABER_N + offset_data + 1 ] ,11),avx_epi64_03 ) , _mm256_slli_epi64(_mm256_and_si256(avx_data[i*SABER_N + offset_data + 2 ],avx_epi64_3f) , 2));

			avx_bytes[offset_byte + 4 ]= _mm256_or_si256(_mm256_and_si256( _mm256_srli_epi64(avx_data[i*SABER_N + offset_data + 2 ] ,6 ),avx_epi64_7f ) , _mm256_slli_epi64(_mm256_and_si256(avx_data[i*SABER_N + offset_data + 3 ],avx_epi64_01) , 7 ));

			avx_bytes[offset_byte + 5 ]= _mm256_and_si256( _mm256_srli_epi64(avx_data[i*SABER_N + offset_data + 3 ] ,1 ),avx_epi64_ff );

			avx_bytes[offset_byte + 6 ]= _mm256_or_si256(_mm256_and_si256( _mm256_srli_epi64(avx_data[i*SABER_N + offset_data + 3 ] ,9 ),avx_epi64_0f ) , _mm256_slli_epi64(_mm256_and_si256(avx_data[i*SABER_N + offset_data + 4 ],avx_epi64_0f) , 4 ));

			avx_bytes[offset_byte + 7 ]= _mm256_and_si256( _mm256_srli_epi64(avx_data[i*SABER_N + offset_data + 4 ] ,4 ),avx_epi64_ff );

			avx_bytes[offset_byte + 8 ]= _mm256_or_si256(_mm256_and_si256( _mm256_srli_epi64(avx_data[i*SABER_N + offset_data + 4 ] ,12),avx_epi64_01 ) , _mm256_slli_epi64(_mm256_and_si256(avx_data[i*SABER_N + offset_data + 5 ],avx_epi64_7f) , 1 ));

			avx_bytes[offset_byte + 9 ]= _mm256_or_si256(_mm256_and_si256( _mm256_srli_epi64(avx_data[i*SABER_N + offset_data + 5 ] ,7 ),avx_epi64_3f ) , _mm256_slli_epi64(_mm256_and_si256(avx_data[i*SABER_N + offset_data + 6 ],avx_epi64_03) , 6 ));

			avx_bytes[offset_byte + 10]= _mm256_and_si256( _mm256_srli_epi64(avx_data[i*SABER_N + offset_data + 6 ] ,2 ),avx_epi64_ff );

			avx_bytes[offset_byte + 11]= _mm256_or_si256(_mm256_and_si256( _mm256_srli_epi64(avx_data[i*SABER_N + offset_data + 6 ] ,10),avx_epi64_07 ) , _mm256_slli_epi64(_mm256_and_si256(avx_data[i*SABER_N + offset_data + 7 ],avx_epi64_1f) , 3 ));

			avx_bytes[offset_byte + 12]= _mm256_and_si256( _mm256_srli_epi64(avx_data[i*SABER_N + offset_data + 7 ] ,5 ),avx_epi64_ff );

		}
	}
	clock_E1 = cpucycles();
	clock_T1 += clock_E1 - clock_S1;

	clock_S2 = cpucycles();
	__m256i temp_byte0, temp_byte1, temp_byte2, temp_byte3;
	uint16_t counter = para/32;
	uint16_t offset_bytes;
	for(m=0;m<counter;++m){
		offset_bytes = (256/8)*m;

		temp_byte0 = _mm256_set_epi8(avx_bytes[offset_bytes + 31][0], avx_bytes[offset_bytes + 30][0], avx_bytes[offset_bytes + 29][0], avx_bytes[offset_bytes + 28][0],
						avx_bytes[offset_bytes + 27][0], avx_bytes[offset_bytes + 26][0], avx_bytes[offset_bytes + 25][0], avx_bytes[offset_bytes + 24][0],
						avx_bytes[offset_bytes + 23][0], avx_bytes[offset_bytes + 22][0], avx_bytes[offset_bytes + 21][0], avx_bytes[offset_bytes + 20][0],
						avx_bytes[offset_bytes + 19][0], avx_bytes[offset_bytes + 18][0], avx_bytes[offset_bytes + 17][0], avx_bytes[offset_bytes + 16][0],
						avx_bytes[offset_bytes + 15][0], avx_bytes[offset_bytes + 14][0], avx_bytes[offset_bytes + 13][0], avx_bytes[offset_bytes + 12][0],
						avx_bytes[offset_bytes + 11][0], avx_bytes[offset_bytes + 10][0], avx_bytes[offset_bytes + 9][0], avx_bytes[offset_bytes + 8][0],
						avx_bytes[offset_bytes + 7][0], avx_bytes[offset_bytes + 6][0], avx_bytes[offset_bytes + 5][0], avx_bytes[offset_bytes + 4][0],
						avx_bytes[offset_bytes + 3][0], avx_bytes[offset_bytes + 2][0], avx_bytes[offset_bytes + 1][0], avx_bytes[offset_bytes + 0][0]);
				_mm256_storeu_si256(avx_bytes0_ptr+m,temp_byte0);
				
				

		temp_byte1 = _mm256_set_epi8(avx_bytes[offset_bytes + 31][1], avx_bytes[offset_bytes + 30][1], avx_bytes[offset_bytes + 29][1], avx_bytes[offset_bytes + 28][1],
						avx_bytes[offset_bytes + 27][1], avx_bytes[offset_bytes + 26][1], avx_bytes[offset_bytes + 25][1], avx_bytes[offset_bytes + 24][1],
						avx_bytes[offset_bytes + 23][1], avx_bytes[offset_bytes + 22][1], avx_bytes[offset_bytes + 21][1], avx_bytes[offset_bytes + 20][1],
						avx_bytes[offset_bytes + 19][1], avx_bytes[offset_bytes + 18][1], avx_bytes[offset_bytes + 17][1], avx_bytes[offset_bytes + 16][1],
						avx_bytes[offset_bytes + 15][1], avx_bytes[offset_bytes + 14][1], avx_bytes[offset_bytes + 13][1], avx_bytes[offset_bytes + 12][1],
						avx_bytes[offset_bytes + 11][1], avx_bytes[offset_bytes + 10][1], avx_bytes[offset_bytes + 9][1], avx_bytes[offset_bytes + 8][1],
						avx_bytes[offset_bytes + 7][1], avx_bytes[offset_bytes + 6][1], avx_bytes[offset_bytes + 5][1], avx_bytes[offset_bytes + 4][1],
						avx_bytes[offset_bytes + 3][1], avx_bytes[offset_bytes + 2][1], avx_bytes[offset_bytes + 1][1], avx_bytes[offset_bytes + 0][1]);
				_mm256_storeu_si256(avx_bytes1_ptr+m,temp_byte1);
						

		temp_byte2 = _mm256_set_epi8(avx_bytes[offset_bytes + 31][2], avx_bytes[offset_bytes + 30][2], avx_bytes[offset_bytes + 29][2], avx_bytes[offset_bytes + 28][2],
						avx_bytes[offset_bytes + 27][2], avx_bytes[offset_bytes + 26][2], avx_bytes[offset_bytes + 25][2], avx_bytes[offset_bytes + 24][2],
						avx_bytes[offset_bytes + 23][2], avx_bytes[offset_bytes + 22][2], avx_bytes[offset_bytes + 21][2], avx_bytes[offset_bytes + 20][2],
						avx_bytes[offset_bytes + 19][2], avx_bytes[offset_bytes + 18][2], avx_bytes[offset_bytes + 17][2], avx_bytes[offset_bytes + 16][2],
						avx_bytes[offset_bytes + 15][2], avx_bytes[offset_bytes + 14][2], avx_bytes[offset_bytes + 13][2], avx_bytes[offset_bytes + 12][2],
						avx_bytes[offset_bytes + 11][2], avx_bytes[offset_bytes + 10][2], avx_bytes[offset_bytes + 9][2], avx_bytes[offset_bytes + 8][2],
						avx_bytes[offset_bytes + 7][2], avx_bytes[offset_bytes + 6][2], avx_bytes[offset_bytes + 5][2], avx_bytes[offset_bytes + 4][2],
						avx_bytes[offset_bytes + 3][2], avx_bytes[offset_bytes + 2][2], avx_bytes[offset_bytes + 1][2], avx_bytes[offset_bytes + 0][2]);
				_mm256_storeu_si256(avx_bytes2_ptr+m,temp_byte2);
						

		temp_byte3 = _mm256_set_epi8(avx_bytes[offset_bytes + 31][3], avx_bytes[offset_bytes + 30][3], avx_bytes[offset_bytes + 29][3], avx_bytes[offset_bytes + 28][3],
						avx_bytes[offset_bytes + 27][3], avx_bytes[offset_bytes + 26][3], avx_bytes[offset_bytes + 25][3], avx_bytes[offset_bytes + 24][3],
						avx_bytes[offset_bytes + 23][3], avx_bytes[offset_bytes + 22][3], avx_bytes[offset_bytes + 21][3], avx_bytes[offset_bytes + 20][3],
						avx_bytes[offset_bytes + 19][3], avx_bytes[offset_bytes + 18][3], avx_bytes[offset_bytes + 17][3], avx_bytes[offset_bytes + 16][3],
						avx_bytes[offset_bytes + 15][3], avx_bytes[offset_bytes + 14][3], avx_bytes[offset_bytes + 13][3], avx_bytes[offset_bytes + 12][3],
						avx_bytes[offset_bytes + 11][3], avx_bytes[offset_bytes + 10][3], avx_bytes[offset_bytes + 9][3], avx_bytes[offset_bytes + 8][3],
						avx_bytes[offset_bytes + 7][3], avx_bytes[offset_bytes + 6][3], avx_bytes[offset_bytes + 5][3], avx_bytes[offset_bytes + 4][3],
						avx_bytes[offset_bytes + 3][3], avx_bytes[offset_bytes + 2][3], avx_bytes[offset_bytes + 1][3], avx_bytes[offset_bytes + 0][3]);
				_mm256_storeu_si256(avx_bytes3_ptr+m,temp_byte3);
	}
	clock_E2 = cpucycles();
	clock_T2 += clock_E2 - clock_S2;
	
}





int main()
{
	
	//test_kem_cca();
	init_avx();
	gene_4bs();
	int i;
	for(i=0;i<10000;i++){
		test_pack_13bit(bs0,data0);
		test_pack_13bit(bs1,data1);
		test_pack_13bit(bs2,data2);
		test_pack_13bit(bs3,data3);
		test_pack_4bitX4(bs0,bs1,bs2,bs3,data0,data1,data2,data3);
	}
	printf("Average load: \t %lu \n",clock_T0/10000);
	printf("Average operate: \t %lu \n",clock_T1/10000);
	printf("Average store: \t %lu \n",clock_T2/10000);
	printf("Average X1: \t %lu \n",clock_T/10000);
	

}

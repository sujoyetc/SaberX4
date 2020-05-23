#include <immintrin.h>
#include <tmmintrin.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#include "../SABER_indcpa.h"
#include "../kem.h"
#include "../api.h"
#include "../poly.h"
//#include "../randombytes.h"
#include "../rng.h"

#include "../cpucycles.c"
#include "../verify.h"


#define SABER_K 3
#define SABER_N 256

unsigned char* stream;
//__m256i pkcl[3][256];
const unsigned char *bs0,*bs1,*bs2,*bs3;
const unsigned char bstring0[960],bstring1[960],bstring2[960],bstring3[960];
uint16_t m=0;
uint16_t data[SABER_K][SABER_N];
__m256i data0[SABER_K][SABER_N/16], data1[SABER_K][SABER_N/16], data2[SABER_K][SABER_N/16], data3[SABER_K][SABER_N/16];
uint16_t dataN0[SABER_K][SABER_N], dataN1[SABER_K][SABER_N], dataN2[SABER_K][SABER_N], dataN3[SABER_K][SABER_N];

uint64_t clock_S, clock_E;
uint64_t clock_S0, clock_E0,clock_S1, clock_E1,clock_S2, clock_E2;
long long int clock_T0,clock_T1,clock_T2, clock_T;


uint64_t clockSum,round;

const unsigned char* gene_bytestream(size_t num_bytes){
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
  bs0 = gene_bytestream(960);
  bs1 = gene_bytestream(960);
  bs2 = gene_bytestream(960);
  bs3 = gene_bytestream(960);
}


void BPVpX4(const unsigned char *bytes0, const unsigned char *bytes1, const unsigned char *bytes2, const unsigned char *bytes3, 
__m256i data0[SABER_K][SABER_N/16], __m256i data1[SABER_K][SABER_N/16], __m256i data2[SABER_K][SABER_N/16], __m256i data3[SABER_K][SABER_N/16]){
	// clock_in = cpucycles();
	clock_S0 = cpucycles();
	__m256i avx_epi64_03, avx_epi64_0f, avx_epi64_3f, avx_epi64_ff;
	avx_epi64_03 = _mm256_set1_epi64x(0x03);
	avx_epi64_0f = _mm256_set1_epi64x(0x0f);
	avx_epi64_3f = _mm256_set1_epi64x(0x3f);
	avx_epi64_ff = _mm256_set1_epi64x(0xff);
	
	
		//printf("BS2POLq store CPUcycle: %d\n",clock_after_store-clock_before_store);
	// clock_out = cpucycles();
	// printf("Load Const CPUcycle: %d\n",clock_out-clock_in);

	// clock_in = cpucycles();
	uint32_t para = (SABER_K-1)*(SABER_N*10)/8+5*(SABER_N/4);
	//uint64_t temp_bytes[para][4];
	__m256i avx_bytes[para];
	__m256i avx_data[SABER_K*SABER_N];
	uint32_t k;
	for(k=0;k<para;++k){
		avx_bytes[k] = _mm256_set_epi64x(bytes3[k], bytes2[k], bytes1[k], bytes0[k]);
	}

	clock_E0 = cpucycles();
	clock_T0 += clock_E0 - clock_S0;
	// clock_out = cpucycles();
	// printf("Load avx CPUcycle: %d\n",clock_out-clock_in);

	// clock_in = cpucycles();

	clock_S1 = cpucycles();
	uint32_t i,j;
	uint32_t offset_data=0,offset_byte=0,offset_byte1=0;
	
	offset_byte=0;
	for(i=0;i<SABER_K;++i){
		offset_byte1=i*(SABER_N*10)/8;
		for(j=0;j<SABER_N/4;++j){
			offset_byte=offset_byte1+5*j;
			offset_data=4*j;
			avx_data[i*SABER_N+offset_data + 0] = _mm256_or_si256(_mm256_and_si256(avx_bytes[ offset_byte + 0 ],avx_epi64_ff), _mm256_slli_epi64(_mm256_and_si256(avx_bytes[ offset_byte + 1 ],avx_epi64_03),8));
			avx_data[i*SABER_N+offset_data + 1] = _mm256_or_si256(_mm256_and_si256(_mm256_srli_epi64(avx_bytes[ offset_byte + 1 ],2),(avx_epi64_3f)),_mm256_slli_epi64(_mm256_and_si256(avx_bytes[ offset_byte + 2 ],avx_epi64_0f),6));
			avx_data[i*SABER_N+offset_data + 2] = _mm256_or_si256(_mm256_and_si256(_mm256_srli_epi64(avx_bytes[ offset_byte + 2 ],4),(avx_epi64_0f)),_mm256_slli_epi64(_mm256_and_si256(avx_bytes[ offset_byte + 3 ],avx_epi64_3f),4));
			avx_data[i*SABER_N+offset_data + 3] = _mm256_or_si256(_mm256_and_si256(_mm256_srli_epi64(avx_bytes[ offset_byte + 3 ],6),(avx_epi64_03)),_mm256_slli_epi64(_mm256_and_si256(avx_bytes[ offset_byte + 4 ],avx_epi64_ff),2));

		}
	}
	clock_E1 = cpucycles();
	clock_T1 += clock_E1 - clock_S1;

	// clock_out = cpucycles();
	// printf("Shift CPUcycle: %d\n",clock_out-clock_in);

	// clock_in = cpucycles();
	clock_S2 = cpucycles();
	uint32_t m,n;
	uint32_t offset_bytes;
	for(m=0;m<SABER_K;++m){
		for(n=0;n<(SABER_N/16);++n){
		offset_bytes = m*SABER_N+16*n;
		data0[m][n] = _mm256_set_epi16(avx_data[offset_bytes + 15][0], avx_data[offset_bytes + 14][0], avx_data[offset_bytes + 13][0], avx_data[offset_bytes + 12][0],
						avx_data[offset_bytes + 11][0], avx_data[offset_bytes + 10][0], avx_data[offset_bytes + 9][0], avx_data[offset_bytes + 8][0],
						avx_data[offset_bytes + 7][0], avx_data[offset_bytes + 6][0], avx_data[offset_bytes + 5][0], avx_data[offset_bytes + 4][0],
						avx_data[offset_bytes + 3][0], avx_data[offset_bytes + 2][0], avx_data[offset_bytes + 1][0], avx_data[offset_bytes + 0][0]);
				//_mm256_storeu_si256(avx_data0_ptr+n+m*16+1,temp_data0);
				
				
		data1[m][n] = _mm256_set_epi16(avx_data[offset_bytes + 15][1], avx_data[offset_bytes + 14][1], avx_data[offset_bytes + 13][1], avx_data[offset_bytes + 12][1],
						avx_data[offset_bytes + 11][1], avx_data[offset_bytes + 10][1], avx_data[offset_bytes + 9][1], avx_data[offset_bytes + 8][1],
						avx_data[offset_bytes + 7][1], avx_data[offset_bytes + 6][1], avx_data[offset_bytes + 5][1], avx_data[offset_bytes + 4][1],
						avx_data[offset_bytes + 3][1], avx_data[offset_bytes + 2][1], avx_data[offset_bytes + 1][1], avx_data[offset_bytes + 0][1]);
						

		data2[m][n] = _mm256_set_epi16(avx_data[offset_bytes + 15][2], avx_data[offset_bytes + 14][2], avx_data[offset_bytes + 13][2], avx_data[offset_bytes + 12][2],
						avx_data[offset_bytes + 11][2], avx_data[offset_bytes + 10][2], avx_data[offset_bytes + 9][2], avx_data[offset_bytes + 8][2],
						avx_data[offset_bytes + 7][2], avx_data[offset_bytes + 6][2], avx_data[offset_bytes + 5][2], avx_data[offset_bytes + 4][2],
						avx_data[offset_bytes + 3][2], avx_data[offset_bytes + 2][2], avx_data[offset_bytes + 1][2], avx_data[offset_bytes + 0][2]);
						

		data3[m][n] = _mm256_set_epi16(avx_data[offset_bytes + 15][3], avx_data[offset_bytes + 14][3], avx_data[offset_bytes + 13][3], avx_data[offset_bytes + 12][3],
						avx_data[offset_bytes + 11][3], avx_data[offset_bytes + 10][3], avx_data[offset_bytes + 9][3], avx_data[offset_bytes + 8][3],
						avx_data[offset_bytes + 7][3], avx_data[offset_bytes + 6][3], avx_data[offset_bytes + 5][3], avx_data[offset_bytes + 4][3],
						avx_data[offset_bytes + 3][3], avx_data[offset_bytes + 2][3], avx_data[offset_bytes + 1][3], avx_data[offset_bytes + 0][3]);
		}
	}
	
	clock_E2 = cpucycles();
	clock_T2 += clock_E2 - clock_S2;

}

void BPVp(const unsigned char *bytes, uint16_t data[SABER_K][SABER_N]){
	
	clock_S = cpucycles();
	uint32_t i,j;
	uint32_t offset_data=0,offset_byte=0,offset_byte1=0;	
	
	offset_byte=0;
	//printf("SABER_K:%d\n",SABER_K);
	for(i=0;i<SABER_K;i++){
		//printf("offset_byte1=i*(SABER_N*10)/8 = %d*(%d*10)/8 = %d\n",i,SABER_N,i*(SABER_N*10)/8);
		offset_byte1=i*(SABER_N*10)/8;
		for(j=0;j<SABER_N/4;j++){
			offset_byte=offset_byte1+5*j;
			offset_data=4*j;
			data[i][offset_data + 0]= ( bytes[ offset_byte + 0 ] & (0xff)) |  ((bytes[ offset_byte + 1 ] & 0x03)<<8);
			data[i][offset_data + 1]= ( (bytes[ offset_byte + 1 ]>>2) & (0x3f)) |  ((bytes[ offset_byte + 2 ] & 0x0f)<<6);		
			data[i][offset_data + 2]= ( (bytes[ offset_byte + 2 ]>>4) & (0x0f)) |  ((bytes[ offset_byte + 3 ] & 0x3f)<<4);
			data[i][offset_data + 3]= ( (bytes[ offset_byte + 3 ]>>6) & (0x03)) |  ((bytes[ offset_byte + 4 ] & 0xff)<<2);		

		}
	}
	clock_E = cpucycles();
	clock_T += clock_E - clock_S;
}


int test_BS2POLVEC(){
	uint64_t i, j, repeat;
	repeat=1000;

	//uint64_t CLOCK1,CLOCK2;
	//uint64_t CLOCK_kp=0;

	for(i=0; i<repeat; i++)
  	{
	    //CLOCK1=cpucycles();	
	    BPVpX4(bs0, bs1, bs2, bs3, data0 ,data1, data2, data3);
	    BPVp(bs0, dataN0);
	    BPVp(bs1, dataN1);
	    BPVp(bs2, dataN2);
	    BPVp(bs3, dataN3);
	    //CLOCK2=cpucycles();	
	    //CLOCK_kp=CLOCK_kp+(CLOCK2-CLOCK1);	
	}

	//clockSum +=CLOCK_kp;
	printf("Average load times: \t %lu \n",clock_T0/1000);
	printf("Average shift times: \t %lu \n",clock_T1/1000);
	printf("Average store times: \t %lu \n",clock_T2/1000);
	printf("Average 1X4 times: \t %lu \n",clock_T/1000);

	
}

int main()
{
	gene_4bs();
	test_BS2POLVEC();
	return 0;
}

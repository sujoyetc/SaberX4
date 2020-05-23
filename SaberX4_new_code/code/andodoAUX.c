#include <stdio.h>
#include <string.h>
#include <immintrin.h>
#include <stdint.h>

void fprint_m256i_16(__m256i avx){
	int i = 0;
	uint16_t* sho = (uint16_t*) &avx;
	FILE *fp;
	fp=fopen("data0.txt","a+");
	if (fp!=NULL){
		for(i=0;i<16;i++){
			fprintf(fp,"[%d]:%d  ",i,sho[i]);
		}	
		fprintf(fp,"\n");
	}
}

/*
FILE write2file(char* filename){
	FILE* fp;
	fp=fopen(strcat(filename,".txt"),"a+");
	if (fp!=NULL) return *fp;
	return NULL;
}
*/
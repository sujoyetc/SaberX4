#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <arpa/inet.h>
#include <errno.h>
#include <immintrin.h>
#define TXT ".txt"
char  buffer[20];

FILE* openfile(char* filename){
	FILE* fp;
	sprintf(buffer,"%s.txt",filename);
	//printf("%s",buffer);
	fp=fopen(buffer,"a+");
	//fp=fopen(strcat(filename,".txt"),"a+");
	if (fp!=NULL) return fp;
	return NULL;
}

void print_byte_as_bits(char val, FILE* fp) {
  int i = 0;
  //FILE* fp = openfile(filename);
  for (i = 7; 0 <= i; i--) {
    fprintf(fp, "%c", (val & (1 << i)) ? '1' : '0');
  }
  //fclose(fp);
}

void print_bits(char * ty, char * val, unsigned char * bytes, size_t num_bytes, char* filename) {
  FILE* fp = openfile(filename);
  if(fp==NULL) printf("End"); 
  fprintf(fp,"(%*s) %*s = [ ", 4, ty, 8, val);
  fprintf(fp,"[ ");
  size_t i = 0;
  for (i = 0; i < num_bytes; i++) {
    print_byte_as_bits(bytes[i],fp);
    fprintf(fp," ");
  }
  fprintf(fp,"]\n");
  fclose(fp);

}



#define SHOW(T,V, filename) do { T x = V; print_bits(#T, #V, (unsigned char*) &x, sizeof(x), filename); } while(0)

/*
int main(){
	SHOW(int,5,"anan");
}
*/

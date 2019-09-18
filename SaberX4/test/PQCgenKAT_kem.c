
//
//  PQCgenKAT_kem.c
//
//  Created by Bassham, Lawrence E (Fed) on 8/29/17.
//  Copyright © 2017 Bassham, Lawrence E (Fed). All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "../rng.h"
#include "../api.h"
#include "../SABER_params.h"
#include "../cpucycles.c"

#define	MAX_MARKER_LEN		50
#define KAT_SUCCESS          0
#define KAT_FILE_OPEN_ERROR -1
#define KAT_DATA_ERROR      -3
#define KAT_CRYPTO_FAILURE  -4


int		FindMarker(FILE *infile, const char *marker);
int		ReadHex(FILE *infile, unsigned char *A, int Length, char *str);
void	fprintBstr(FILE *fp, char *S, unsigned char *A, unsigned long long L);

int
main()
{
    char                fn_req[32], fn_rsp[32];
    FILE                *fp_req, *fp_rsp;
    unsigned char       seed[48];
    unsigned char       entropy_input[48];
    unsigned char       ct_0[CRYPTO_CIPHERTEXTBYTES], ss_0[CRYPTO_BYTES], ss1_0[CRYPTO_BYTES];
    unsigned char       ct_1[CRYPTO_CIPHERTEXTBYTES], ss_1[CRYPTO_BYTES], ss1_1[CRYPTO_BYTES];
    unsigned char       ct_2[CRYPTO_CIPHERTEXTBYTES], ss_2[CRYPTO_BYTES], ss1_2[CRYPTO_BYTES];
    unsigned char       ct_3[CRYPTO_CIPHERTEXTBYTES], ss_3[CRYPTO_BYTES], ss1_3[CRYPTO_BYTES];
    int                 count;
    int                 done;
    unsigned char       pk0[CRYPTO_PUBLICKEYBYTES], sk0[CRYPTO_SECRETKEYBYTES];
    unsigned char       pk1[CRYPTO_PUBLICKEYBYTES], sk1[CRYPTO_SECRETKEYBYTES];
    unsigned char       pk2[CRYPTO_PUBLICKEYBYTES], sk2[CRYPTO_SECRETKEYBYTES];
    unsigned char       pk3[CRYPTO_PUBLICKEYBYTES], sk3[CRYPTO_SECRETKEYBYTES];

    int                 ret_val;
    int i;

    // Create the REQUEST file
    sprintf(fn_req, "PQCkemKAT_%d.req", CRYPTO_SECRETKEYBYTES);
    if ( (fp_req = fopen(fn_req, "w")) == NULL ) {
        printf("Couldn't open <%s> for write\n", fn_req);
        return KAT_FILE_OPEN_ERROR;
    }
    sprintf(fn_rsp, "PQCkemKAT_%d.rsp", CRYPTO_SECRETKEYBYTES);
    if ( (fp_rsp = fopen(fn_rsp, "w")) == NULL ) {
        printf("Couldn't open <%s> for write\n", fn_rsp);
        return KAT_FILE_OPEN_ERROR;
    }
    
    for (i=0; i<48; i++)
        entropy_input[i] = i;

    randombytes_init(entropy_input, NULL, 256);
    for (i=0; i<100; i++) {
        fprintf(fp_req, "count = %d\n", i);
        randombytes(seed, 48);
        fprintBstr(fp_req, "seed = ", seed, 48);
        fprintf(fp_req, "pk =\n");
        fprintf(fp_req, "sk =\n");
        fprintf(fp_req, "ct =\n");
        fprintf(fp_req, "ss =\n\n");
    }
    fclose(fp_req);
    
    //Create the RESPONSE file based on what's in the REQUEST file
    if ( (fp_req = fopen(fn_req, "r")) == NULL ) {
        printf("Couldn't open <%s> for read\n", fn_req);
        return KAT_FILE_OPEN_ERROR;
    }
    
    fprintf(fp_rsp, "# %s\n\n", CRYPTO_ALGNAME);
    done = 0;
    do {
        if ( FindMarker(fp_req, "count = ") )
            fscanf(fp_req, "%d", &count);
        else {
            done = 1;
            break;
        }
        fprintf(fp_rsp, "count = %d\n", count);
        
        if ( !ReadHex(fp_req, seed, 48, "seed = ") ) {
            printf("ERROR: unable to read 'seed' from <%s>\n", fn_req);
            return KAT_DATA_ERROR;
        }
        fprintBstr(fp_rsp, "seed = ", seed, 48);
        
        randombytes_init(seed, NULL, 256);
        
        // Generate the public/private keypair
        if ( (ret_val = crypto_kem_keypair(pk0,sk0, pk1,sk1, pk2,sk2, pk3,sk3)) != 0) {
            printf("crypto_kem_keypair returned <%d>\n", ret_val);
            return KAT_CRYPTO_FAILURE;
        }
        fprintBstr(fp_rsp, "pk0 = ", pk0, CRYPTO_PUBLICKEYBYTES);
        fprintBstr(fp_rsp, "sk0 = ", sk0, CRYPTO_SECRETKEYBYTES);
        fprintBstr(fp_rsp, "pk1 = ", pk1, CRYPTO_PUBLICKEYBYTES);
        fprintBstr(fp_rsp, "sk1 = ", sk1, CRYPTO_SECRETKEYBYTES);
        fprintBstr(fp_rsp, "pk2 = ", pk1, CRYPTO_PUBLICKEYBYTES);
        fprintBstr(fp_rsp, "sk2 = ", sk1, CRYPTO_SECRETKEYBYTES);
        fprintBstr(fp_rsp, "pk3 = ", pk1, CRYPTO_PUBLICKEYBYTES);
        fprintBstr(fp_rsp, "sk3 = ", sk1, CRYPTO_SECRETKEYBYTES);

        
        if ( (ret_val = crypto_kem_enc(ct_0,ss_0,pk0, ct_1,ss_1,pk1, ct_2,ss_2,pk2, ct_3,ss_3,pk3)) != 0) {
            printf("crypto_kem_enc returned <%d>\n", ret_val);
            return KAT_CRYPTO_FAILURE;
        }
        fprintBstr(fp_rsp, "ct0 = ", ct_0, CRYPTO_CIPHERTEXTBYTES);
        fprintBstr(fp_rsp, "ss0 = ", ss_0, CRYPTO_BYTES);
        fprintBstr(fp_rsp, "ct1 = ", ct_0, CRYPTO_CIPHERTEXTBYTES);
        fprintBstr(fp_rsp, "ss1 = ", ss_0, CRYPTO_BYTES);
        fprintBstr(fp_rsp, "ct2 = ", ct_0, CRYPTO_CIPHERTEXTBYTES);
        fprintBstr(fp_rsp, "ss2 = ", ss_0, CRYPTO_BYTES);
        fprintBstr(fp_rsp, "ct3 = ", ct_0, CRYPTO_CIPHERTEXTBYTES);
        fprintBstr(fp_rsp, "ss3 = ", ss_0, CRYPTO_BYTES);
        
        fprintf(fp_rsp, "\n");
        
        if ( (ret_val = crypto_kem_dec(ss1_0,ct_0,sk0, ss1_1,ct_1,sk1, ss1_2,ct_2,sk2, ss1_3,ct_3,sk3)) != 0) {
            printf("crypto_kem_dec returned <%d>\n", ret_val);
            return KAT_CRYPTO_FAILURE;
        }
        
        if ( memcmp(ss_0, ss1_0, CRYPTO_BYTES) || memcmp(ss_1, ss1_1, CRYPTO_BYTES) || memcmp(ss_2, ss1_2, CRYPTO_BYTES) || memcmp(ss_3, ss1_3, CRYPTO_BYTES)) {
            printf("crypto_kem_dec returned bad 'ss' value\n");
            return KAT_CRYPTO_FAILURE;
        }

    } while ( !done );
    
    fclose(fp_req);
    fclose(fp_rsp);

    return KAT_SUCCESS;
}



//
// ALLOW TO READ HEXADECIMAL ENTRY (KEYS, DATA, TEXT, etc.)
//
//
// ALLOW TO READ HEXADECIMAL ENTRY (KEYS, DATA, TEXT, etc.)
//
int
FindMarker(FILE *infile, const char *marker)
{
	char	line[MAX_MARKER_LEN];
	int		i, len;
	int curr_line;

	len = (int)strlen(marker);
	if ( len > MAX_MARKER_LEN-1 )
		len = MAX_MARKER_LEN-1;

	for ( i=0; i<len; i++ )
	  {
	    curr_line = fgetc(infile);
	    line[i] = curr_line;
	    if (curr_line == EOF )
	      return 0;
	  }
	line[len] = '\0';

	while ( 1 ) {
		if ( !strncmp(line, marker, len) )
			return 1;

		for ( i=0; i<len-1; i++ )
			line[i] = line[i+1];
		curr_line = fgetc(infile);
		line[len-1] = curr_line;
		if (curr_line == EOF )
		    return 0;
		line[len] = '\0';
	}

	// shouldn't get here
	return 0;
}

//
// ALLOW TO READ HEXADECIMAL ENTRY (KEYS, DATA, TEXT, etc.)
//
int
ReadHex(FILE *infile, unsigned char *A, int Length, char *str)
{
	int			i, ch, started;
	unsigned char	ich;

	if ( Length == 0 ) {
		A[0] = 0x00;
		return 1;
	}
	memset(A, 0x00, Length);
	started = 0;
	if ( FindMarker(infile, str) )
		while ( (ch = fgetc(infile)) != EOF ) {
			if ( !isxdigit(ch) ) {
				if ( !started ) {
					if ( ch == '\n' )
						break;
					else
						continue;
				}
				else
					break;
			}
			started = 1;
			if ( (ch >= '0') && (ch <= '9') )
				ich = ch - '0';
			else if ( (ch >= 'A') && (ch <= 'F') )
				ich = ch - 'A' + 10;
			else if ( (ch >= 'a') && (ch <= 'f') )
				ich = ch - 'a' + 10;
            else // shouldn't ever get here
                ich = 0;
			
			for ( i=0; i<Length-1; i++ )
				A[i] = (A[i] << 4) | (A[i+1] >> 4);
			A[Length-1] = (A[Length-1] << 4) | ich;
		}
	else
		return 0;

	return 1;
}

void
fprintBstr(FILE *fp, char *S, unsigned char *A, unsigned long long L)
{
	unsigned long long  i;

	fprintf(fp, "%s", S);

	for ( i=0; i<L; i++ )
		fprintf(fp, "%02X", A[i]);

	if ( L == 0 )
		fprintf(fp, "00");

	fprintf(fp, "\n");
}


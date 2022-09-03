#include <stdio.h>
#include <string.h>
#include <inttypes.h>

// ENDIANNESS SWAP
void uint64swe(uint64_t  * inttoswap){

	char *p= (char*) inttoswap; char t;
	t = p[7]; p[7] = p[0] ; p[0] = t;
	t = p[6]; p[6] = p[1] ; p[1] = t;
	t = p[5]; p[5] = p[2] ; p[2] = t;
	t = p[4]; p[4] = p[3] ; p[3] = t;

}
void uint32swe(uint32_t  * inttoswap){

	char *p= (char*) inttoswap; char t;
	t = p[3]; p[3] = p[0] ; p[0] = t;
	t = p[2]; p[2] = p[1] ; p[1] = t;


}
void uint16swe(uint16_t  * inttoswap){

	char *p= (char*) inttoswap; char t;
	t = p[1]; p[1] = p[0] ; p[0] = t;

}
typedef struct ILDG_header_t
{
	uint32_t magic_no;
	uint16_t version;
	uint16_t mbme_flag;
	uint64_t data_length;
	char type[128];
}ILDG_header;


ILDG_header headers[15];

int main(){
	FILE *fg;

	const int len=10000;
	char header[len];
	fg = fopen("bioparco","r+");
	int allfound = 0;

	while(fgets(header,len,fg)!= NULL && !allfound){
		char * strfnx = strstr(header,"<lx>");
		char * strfny = strstr(header,"<ly>");
		char * strfnz = strstr(header,"<lz>");
		char * strfnt = strstr(header,"<lt>");
		int nx,ny,nz,nt;
		allfound = 1;
		if(strfnx!= NULL){sscanf(strfnx," <lx> %d </lx> ",&nx);printf("nx: %d\n",nx);}else{allfound = 0;}
		if(strfny!= NULL){sscanf(strfny," <ly> %d </ly> ",&ny);printf("ny: %d\n",ny);}else{allfound = 0;}
		if(strfnz!= NULL){sscanf(strfnz," <lz> %d </lz> ",&nz);printf("nz: %d\n",nz);}else{allfound = 0;}
		if(strfnt!= NULL){sscanf(strfnt," <lt> %d </lt> ",&nt);printf("nt: %d\n",nt);}else{allfound = 0;}

	}

	fseek(fg,0,SEEK_SET);


	int i;
	for(i=0; i<5; i++){
    fread(&headers[i],sizeof(ILDG_header),1,fg);

    uint64swe(&headers[i].data_length);
    uint32swe(&headers[i].magic_no);
    printf("header type: %s, data length:%" PRIu64 ", magic number: %" PRIu32 "\n", headers[i].type, headers[i].data_length,headers[i].magic_no);
    fseek(fg,headers[i].data_length,SEEK_CUR);
	}

	fclose(fg);

	return 0;

}

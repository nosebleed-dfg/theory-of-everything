#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

static const uint32_t KK[64]={0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
#define RR(x,n) (((x)>>(n))|((x)<<(32-(n))))
static void sha256_block(const uint8_t *bl,uint32_t *st){uint32_t W[64],a,b,c,d,e,f,g,h,T1,T2;int i;for(i=0;i<16;i++)W[i]=((uint32_t)bl[i*4]<<24)|((uint32_t)bl[i*4+1]<<16)|((uint32_t)bl[i*4+2]<<8)|bl[i*4+3];for(i=16;i<64;i++)W[i]=(RR(W[i-2],17)^RR(W[i-2],19)^(W[i-2]>>10))+W[i-7]+(RR(W[i-15],7)^RR(W[i-15],18)^(W[i-15]>>3))+W[i-16];a=st[0];b=st[1];c=st[2];d=st[3];e=st[4];f=st[5];g=st[6];h=st[7];for(i=0;i<64;i++){T1=h+(RR(e,6)^RR(e,11)^RR(e,25))+((e&f)^((~e)&g))+KK[i]+W[i];T2=(RR(a,2)^RR(a,13)^RR(a,22))+((a&b)^(a&c)^(b&c));h=g;g=f;f=e;e=d+T1;d=c;c=b;b=a;a=T1+T2;}st[0]+=a;st[1]+=b;st[2]+=c;st[3]+=d;st[4]+=e;st[5]+=f;st[6]+=g;st[7]+=h;}
static void sha256(const uint8_t *m,size_t l,uint8_t *o){uint32_t st[8]={0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};uint8_t bl[64];size_t i,n=l/64;for(i=0;i<n;i++)sha256_block(m+i*64,st);size_t r=l%64;memset(bl,0,64);memcpy(bl,m+n*64,r);bl[r]=0x80;if(r>=56){sha256_block(bl,st);memset(bl,0,64);}uint64_t bits=l*8;for(i=0;i<8;i++)bl[63-i]=(bits>>(i*8))&0xff;sha256_block(bl,st);for(i=0;i<8;i++){o[i*4]=(st[i]>>24)&0xff;o[i*4+1]=(st[i]>>16)&0xff;o[i*4+2]=(st[i]>>8)&0xff;o[i*4+3]=st[i]&0xff;}}
static void sha256d(const uint8_t *m,size_t l,uint8_t *o){uint8_t t[32];sha256(m,l,t);sha256(t,32,o);}

int main(int argc, char **argv) {
    if (argc < 2) { printf("Usage: %s <header76hex> [start_nonce]\n", argv[0]); return 1; }
    
    uint8_t hdr[80];
    for (int i=0;i<76;i++) sscanf(argv[1]+2*i,"%2hhx",&hdr[i]);
    
    uint32_t start = 0;
    if (argc >= 3) start = (uint32_t)atoll(argv[2]);
    
    /* Pool diff 10000 target */
    uint8_t target_prefix[6] = {0, 0, 0, 0, 0, 0}; /* 45+ zeros = first 5.5 bytes zero */
    
    printf("sequential from nonce=%u\n", start);
    fflush(stdout);
    
    clock_t t0 = clock();
    
    for (uint64_t n = start; n < (uint64_t)start + 100000000ULL; n++) {
        uint32_t nonce = (uint32_t)(n & 0xFFFFFFFF);
        hdr[76]=nonce&0xff; hdr[77]=(nonce>>8)&0xff;
        hdr[78]=(nonce>>16)&0xff; hdr[79]=(nonce>>24)&0xff;
        
        uint8_t hash[32];
        sha256d(hdr, 80, hash);
        
        /* Quick check: bytes 31,30,29,28,27 must be 0 for 40+ zeros */
        if (hash[31]==0 && hash[30]==0 && hash[29]==0 && hash[28]==0 && hash[27]==0) {
            /* Full check */
            int z = 0;
            for (int i=31;i>=0;i--) {
                if(hash[i]==0) z+=8;
                else {uint8_t b=hash[i];while(!(b&0x80)){z++;b<<=1;}break;}
            }
            if (z >= 40) {
                double elapsed = (double)(clock()-t0)/CLOCKS_PER_SEC;
                printf("FOUND nonce=%u zeros=%d time=%.2fs rate=%.0f MH/s\n",
                       nonce, z, elapsed, (n-start)/elapsed/1000000);
                printf("hash=");
                for(int i=31;i>=0;i--) printf("%02x",hash[i]);
                printf("\n");
                fflush(stdout);
                if (z >= 45) {
                    printf("*** VALID SHARE ***\n");
                    return 0;
                }
            }
        }
        
        if ((n-start) % 10000000 == 0 && n > start) {
            double elapsed = (double)(clock()-t0)/CLOCKS_PER_SEC;
            printf("%lluM (%.1fs %.1f MH/s)\n", 
                   (unsigned long long)(n-start)/1000000, elapsed, (n-start)/elapsed/1000000);
            fflush(stdout);
        }
    }
    
    double elapsed = (double)(clock()-t0)/CLOCKS_PER_SEC;
    printf("done. 100M in %.1fs (%.1f MH/s)\n", elapsed, 100000000.0/elapsed);
    return 1;
}

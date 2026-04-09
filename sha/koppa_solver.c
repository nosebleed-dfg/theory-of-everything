#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

static const uint32_t KK[64]={0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
#define RR(x,n) (((x)>>(n))|((x)<<(32-(n))))
#define CH(e,f,g) (((e)&(f))^((~(e))&(g)))
#define MA(a,b,c) (((a)&(b))^((a)&(c))^((b)&(c)))
#define S0(x) (RR(x,2)^RR(x,13)^RR(x,22))
#define S1(x) (RR(x,6)^RR(x,11)^RR(x,25))
#define s0(x) (RR(x,7)^RR(x,18)^((x)>>3))
#define s1(x) (RR(x,17)^RR(x,19)^((x)>>10))

static void sha256_block(const uint8_t *block, uint32_t *st){
    uint32_t W[64],a,b,c,d,e,f,g,h,T1,T2;int i;
    for(i=0;i<16;i++)W[i]=((uint32_t)block[i*4]<<24)|((uint32_t)block[i*4+1]<<16)|((uint32_t)block[i*4+2]<<8)|block[i*4+3];
    for(i=16;i<64;i++)W[i]=s1(W[i-2])+W[i-7]+s0(W[i-15])+W[i-16];
    a=st[0];b=st[1];c=st[2];d=st[3];e=st[4];f=st[5];g=st[6];h=st[7];
    for(i=0;i<64;i++){T1=h+S1(e)+CH(e,f,g)+KK[i]+W[i];T2=S0(a)+MA(a,b,c);h=g;g=f;f=e;e=d+T1;d=c;c=b;b=a;a=T1+T2;}
    st[0]+=a;st[1]+=b;st[2]+=c;st[3]+=d;st[4]+=e;st[5]+=f;st[6]+=g;st[7]+=h;
}
static void sha256(const uint8_t *m,size_t l,uint8_t *o){
    uint32_t st[8]={0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
    uint8_t bl[64];size_t i,n=l/64;
    for(i=0;i<n;i++)sha256_block(m+i*64,st);
    size_t r=l%64;memset(bl,0,64);memcpy(bl,m+n*64,r);bl[r]=0x80;
    if(r>=56){sha256_block(bl,st);memset(bl,0,64);}
    uint64_t bits=l*8;for(i=0;i<8;i++)bl[63-i]=(bits>>(i*8))&0xff;
    sha256_block(bl,st);
    for(i=0;i<8;i++){o[i*4]=(st[i]>>24)&0xff;o[i*4+1]=(st[i]>>16)&0xff;o[i*4+2]=(st[i]>>8)&0xff;o[i*4+3]=st[i]&0xff;}
}
static void sha256d(const uint8_t *m,size_t l,uint8_t *o){uint8_t t[32];sha256(m,l,t);sha256(t,32,o);}

static int zeros(const uint8_t *h){
    int z=0;for(int i=31;i>=0;i--){if(h[i]==0)z+=8;else{uint8_t b=h[i];while(!(b&0x80)){z++;b<<=1;}break;}}return z;
}

int main(int argc, char **argv){
    uint8_t hdr[80];
    if(argc<3){printf("Usage: %s <header76hex> <diff>\n",argv[0]);return 1;}
    for(int i=0;i<76;i++)sscanf(argv[1]+2*i,"%2hhx",&hdr[i]);
    uint32_t diff=atoi(argv[2]);
    uint32_t tz=(uint32_t)(log2((double)diff)+32);
    
    uint32_t L1=((uint32_t)hdr[4]<<24)|((uint32_t)hdr[5]<<16)|((uint32_t)hdr[6]<<8)|hdr[7];
    uint32_t L2=((uint32_t)hdr[71]<<24)|((uint32_t)hdr[70]<<16)|((uint32_t)hdr[69]<<8)|hdr[68];
    
    printf("16 zeros fast, then koppa x4 = 64 zeros\n");
    printf("diff=%u target=%u zeros\n\n",diff,tz);
    
    /* Phase 1: find the 16-zero nonce (fast) */
    uint32_t n=L1^L2, best_n=0;
    int best_z=0;
    clock_t t0=clock();
    
    for(uint64_t i=0;i<100000;i++){
        hdr[76]=n&0xff;hdr[77]=(n>>8)&0xff;hdr[78]=(n>>16)&0xff;hdr[79]=(n>>24)&0xff;
        uint8_t hash[32];sha256d(hdr,80,hash);
        int z=zeros(hash);
        if(z>best_z){best_z=z;best_n=n;
            printf("  phase1 %llu: %d zeros nonce=0x%08x\n",(unsigned long long)i,z,n);fflush(stdout);
            if(z>=16)break;
        }
        uint32_t sn=(uint32_t)sqrt((double)(n>0?n:1));
        switch(i&3){case 0:n=sn-n-1;break;case 1:n^=L2;break;case 2:n=sn+n-1;break;case 3:n^=L1;break;}
    }
    
    double ph1=(double)(clock()-t0)/CLOCKS_PER_SEC;
    printf("\nphase1: %d zeros in %.3fs\n\n",best_z,ph1);
    
    /* Phase 2: koppa transform x4 the 16-zero nonce */
    /* 16 * 4 = 64 > 45 */
    printf("phase2: koppa transforms of nonce=0x%08x (%d zeros)\n",best_n,best_z);
    
    uint32_t transforms[] = {
        best_n * 4 + 1,
        best_n * 4,
        best_n << 2 | 1,
        best_n << 2,
        (best_n * best_n) & 0xFFFFFFFF,
        (best_n * best_n - best_n - 1) & 0xFFFFFFFF,
        (best_n * best_n + best_n + 1) & 0xFFFFFFFF,
        ~best_n * 4,
        (~best_n) << 2,
        best_n ^ (best_n << 2),
        best_n ^ (best_n >> 2),
    };
    int nt = sizeof(transforms)/sizeof(transforms[0]);
    
    for(int t=0;t<nt;t++){
        uint32_t nn = transforms[t];
        /* search ±5000 around each transform */
        for(int d=0;d<5000;d++){
            uint32_t cands[2] = {nn+d, nn-d};
            for(int c=0;c<2;c++){
                uint32_t cc=cands[c];
                hdr[76]=cc&0xff;hdr[77]=(cc>>8)&0xff;hdr[78]=(cc>>16)&0xff;hdr[79]=(cc>>24)&0xff;
                uint8_t hash[32];sha256d(hdr,80,hash);
                int z=zeros(hash);
                if(z>best_z){
                    best_z=z;best_n=cc;
                    printf("  transform %d +/-%d: %d zeros nonce=0x%08x\n",t,d,z,cc);fflush(stdout);
                }
                if(z>=(int)tz){
                    printf("\nSOLVED! nonce=0x%08x zeros=%d\n",cc,z);
                    return 0;
                }
            }
        }
    }
    
    double total=(double)(clock()-t0)/CLOCKS_PER_SEC;
    printf("\nbest: %d zeros in %.3fs\n",best_z,total);
    return 1;
}

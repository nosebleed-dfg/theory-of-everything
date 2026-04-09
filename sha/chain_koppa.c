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
static int zeros(const uint8_t *h){int z=0;for(int i=31;i>=0;i--){if(h[i]==0)z+=8;else{uint8_t b=h[i];while(!(b&0x80)){z++;b<<=1;}break;}}return z;}

static int test_nonce(uint8_t *hdr, uint32_t n, int tz) {
    hdr[76]=n&0xff;hdr[77]=(n>>8)&0xff;hdr[78]=(n>>16)&0xff;hdr[79]=(n>>24)&0xff;
    uint8_t hash[32]; sha256d(hdr,80,hash);
    return zeros(hash);
}

static uint32_t find_best(uint8_t *hdr, uint32_t center, int radius, int *best_z) {
    uint32_t best_n = center;
    for (int d = 0; d < radius; d++) {
        uint32_t cands[2] = {center+d, center-d};
        for (int c = 0; c < 2; c++) {
            int z = test_nonce(hdr, cands[c], 0);
            if (z > *best_z) { *best_z = z; best_n = cands[c]; }
        }
    }
    return best_n;
}

int main(int argc, char **argv) {
    if (argc < 3) { printf("Usage: %s <hdr76hex> <diff>\n", argv[0]); return 1; }
    uint8_t hdr[80];
    for (int i=0;i<76;i++) sscanf(argv[1]+2*i,"%2hhx",&hdr[i]);
    uint32_t diff = atoi(argv[2]);
    int tz = (int)(log2((double)diff) + 32);
    
    uint32_t L1=((uint32_t)hdr[4]<<24)|((uint32_t)hdr[5]<<16)|((uint32_t)hdr[6]<<8)|hdr[7];
    uint32_t L2=((uint32_t)hdr[71]<<24)|((uint32_t)hdr[70]<<16)|((uint32_t)hdr[69]<<8)|hdr[68];
    
    printf("CHAINED KOPPA SOLVER\ndiff=%u target=%d zeros\n\n", diff, tz);
    clock_t t0 = clock();
    
    /* Phase 1: gravity fall to 16 zeros */
    uint32_t n = L1 ^ L2;
    int best_z = 0;
    for (uint64_t i = 0; i < 100000; i++) {
        int z = test_nonce(hdr, n, tz);
        if (z > best_z) { best_z = z; printf("  fall %llu: %d zeros 0x%08x\n",(unsigned long long)i,z,n); fflush(stdout); }
        if (z >= tz) { printf("SOLVED in fall!\n"); return 0; }
        if (z >= 16) break;
        uint32_t sn = (uint32_t)sqrt((double)(n>0?n:1));
        switch(i&3){case 0:n=sn-n-1;break;case 1:n^=L2;break;case 2:n=sn+n-1;break;case 3:n^=L1;break;}
    }
    printf("phase1: %d zeros (%.3fs)\n\n", best_z, (double)(clock()-t0)/CLOCKS_PER_SEC);
    
    /* Phase 2: chain koppa transforms, each should add ~7 zeros */
    /* 16 -> 23 -> 30 -> 37 -> 44 -> 45+ */
    uint32_t current = n;
    for (int koppa_round = 0; koppa_round < 6; koppa_round++) {
        printf("koppa round %d (current=%d zeros):\n", koppa_round, best_z);
        
        /* Try all koppa transforms, keep the best */
        uint32_t transforms[8];
        transforms[0] = ~current * 4;           /* mirror * 4 (won last time) */
        transforms[1] = current * 4 + 1;        /* scale + axiom */
        transforms[2] = (current * current) - current - 1;  /* quadratic axiom */
        transforms[3] = ~current << 2;           /* mirror shift */
        transforms[4] = current ^ (current << 2); /* XOR shift */
        transforms[5] = current ^ (~current >> 2); /* XOR mirror shift */
        transforms[6] = (current + (~current << 1)); /* add mirror*2 */
        transforms[7] = current ^ L1 ^ L2;       /* XOR both L-ends */
        
        int round_best = best_z;
        uint32_t round_best_n = current;
        
        for (int t = 0; t < 8; t++) {
            int tz2 = best_z;
            uint32_t bn = find_best(hdr, transforms[t], 5000, &tz2);
            if (tz2 > round_best) {
                round_best = tz2;
                round_best_n = bn;
                printf("  transform %d: %d zeros nonce=0x%08x\n", t, tz2, bn);
                fflush(stdout);
            }
            if (tz2 >= tz) {
                printf("\nSOLVED! nonce=0x%08x zeros=%d time=%.3fs\n", bn, tz2,
                       (double)(clock()-t0)/CLOCKS_PER_SEC);
                return 0;
            }
        }
        
        best_z = round_best;
        current = round_best_n;
        printf("  round %d best: %d zeros\n\n", koppa_round, best_z);
    }
    
    printf("final: %d zeros (%.3fs)\n", best_z, (double)(clock()-t0)/CLOCKS_PER_SEC);
    return 1;
}

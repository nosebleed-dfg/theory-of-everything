/*
 * AXIOM_MINER — Bitcoin SHA-256d miner, golden center + pentagon band search
 *
 * Strategy:
 *   T_291: nonce_center = F(291)*seed + F(290) mod 2^32, seed from block-1 IV
 *   Golden spiral: 0, +1, -1, +2, -2, ... from center
 *   Passive K round 3: K[3]=0xe9b5dba5 (band 4), nonce enters here directly
 *   Rounds 0-2 are nonce-independent -> precompute once per block
 *
 * Passive K-band rounds (K > 0xCCCCCCCC): 3,8,16,17,29,44,45,46
 * At round 3: W[3]=nonce. Deterministic passive carry.
 *
 * n^2 as lever: phi (force) x psi (opposite-face hinge) = -1 (unit moment).
 * Distance phi-psi = sqrt(5) = dodecahedron face-to-face. Gain = +1 = axiom.
 *
 * Usage:
 *   axiom_miner demo                         genesis verification + easy mine
 *   axiom_miner mine [bits] [prev_hex]       mine with coinbase tag
 *   axiom_miner <80-byte-header-hex> [steps] mine given header
 *
 * CONFIGURE:
 *   Set COINBASE_MSG to your tag.
 *   Set REWARD_ADDR to your Bitcoin P2PKH address to claim block rewards.
 *   Leave REWARD_ADDR "" to use OP_RETURN (coins unspendable — use for testing).
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

/* ── Configure these ─────────────────────────────────────────────────────── */
#define COINBASE_MSG  "YOUR_TAG / YOUR_GROUP"
#define REWARD_ADDR   ""    /* P2PKH address e.g. "1A1zP1eP5QGefi2DMPTfTL5SLmv7Divf..." */
#define REWARD_SAT    5000000000ULL  /* 50 BTC */

/* ── Framework constants ─────────────────────────────────────────────────── */
#define F291  0x85824b02UL
#define F290  0xc3e3c441UL
#define BAND  0x33333333UL

/* ── SHA-256 round constants ─────────────────────────────────────────────── */
static const uint32_t SHA_K[64] = {
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};

static const uint32_t SHA_H0[8] = {
    0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
    0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19
};

/* ── Primitives ──────────────────────────────────────────────────────────── */
#define ROTR(x,n)  (((x)>>(n))|((x)<<(32-(n))))
#define CH(e,f,g)  (((e)&(f))^(~(e)&(g)))
#define MAJ(a,b,c) (((a)&(b))^((a)&(c))^((b)&(c)))
#define EP0(x)     (ROTR(x,2)^ROTR(x,13)^ROTR(x,22))
#define EP1(x)     (ROTR(x,6)^ROTR(x,11)^ROTR(x,25))
#define SIG0(x)    (ROTR(x,7)^ROTR(x,18)^((x)>>3))
#define SIG1(x)    (ROTR(x,17)^ROTR(x,19)^((x)>>10))
#define WORD_BAND(w) (((uint64_t)(w)*5)>>32)

/* ── SHA-256 context ─────────────────────────────────────────────────────── */
typedef struct {
    uint32_t state[8];
    uint64_t count;
    uint8_t  buf[64];
    int      buf_len;
} SHA256_CTX;

static void sha256_compress_block(uint32_t st[8], const uint8_t data[64])
{
    uint32_t W[64],a,b,c,d,e,f,g,h,t1,t2; int i;
    for(i=0;i<16;i++)
        W[i]=((uint32_t)data[i*4]<<24)|((uint32_t)data[i*4+1]<<16)
            |((uint32_t)data[i*4+2]<<8)|data[i*4+3];
    for(i=16;i<64;i++) W[i]=SIG1(W[i-2])+W[i-7]+SIG0(W[i-15])+W[i-16];
    a=st[0];b=st[1];c=st[2];d=st[3];e=st[4];f=st[5];g=st[6];h=st[7];
    for(i=0;i<64;i++){t1=h+EP1(e)+CH(e,f,g)+SHA_K[i]+W[i];t2=EP0(a)+MAJ(a,b,c);
        h=g;g=f;f=e;e=d+t1;d=c;c=b;b=a;a=t1+t2;}
    st[0]+=a;st[1]+=b;st[2]+=c;st[3]+=d;st[4]+=e;st[5]+=f;st[6]+=g;st[7]+=h;
}

static void sha256_init(SHA256_CTX *ctx)
{
    memcpy(ctx->state, SHA_H0, 32);
    ctx->count = 0; ctx->buf_len = 0;
}

static void sha256_update(SHA256_CTX *ctx, const uint8_t *data, int len)
{
    ctx->count += (uint64_t)len * 8;
    while (len > 0) {
        int take = 64 - ctx->buf_len;
        if (take > len) take = len;
        memcpy(ctx->buf + ctx->buf_len, data, take);
        ctx->buf_len += take; data += take; len -= take;
        if (ctx->buf_len == 64) {
            sha256_compress_block(ctx->state, ctx->buf);
            ctx->buf_len = 0;
        }
    }
}

static void sha256_final(SHA256_CTX *ctx, uint8_t out[32])
{
    uint8_t pad[72]; int i, pad_len;
    pad[0] = 0x80;
    pad_len = ((55 - ctx->buf_len) & 63) + 1;
    memset(pad+1, 0, pad_len-1);
    uint64_t bits = ctx->count;
    for(i=0;i<8;i++) pad[pad_len+i] = (bits >> (56-8*i)) & 0xFF;
    sha256_update(ctx, pad, pad_len+8);
    for(i=0;i<8;i++){
        out[i*4+0]=(ctx->state[i]>>24)&0xFF; out[i*4+1]=(ctx->state[i]>>16)&0xFF;
        out[i*4+2]=(ctx->state[i]>> 8)&0xFF; out[i*4+3]=(ctx->state[i]    )&0xFF;
    }
}

static void sha256d_bytes(const uint8_t *data, int len, uint8_t out[32])
{
    SHA256_CTX ctx; uint8_t h1[32];
    sha256_init(&ctx); sha256_update(&ctx,data,len); sha256_final(&ctx,h1);
    sha256_init(&ctx); sha256_update(&ctx,h1,32); sha256_final(&ctx,out);
}

/* ── Base58Check decode (P2PKH address → 20-byte hash160) ───────────────── */
static const char B58_ALPHA[] =
    "123456789ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz";

static int base58check_decode(const char *addr, uint8_t hash160[20])
{
    uint8_t buf[25]; memset(buf,0,25);
    int alen=(int)strlen(addr);
    for(int i=0;i<alen;i++){
        const char *p=strchr(B58_ALPHA,addr[i]);
        if(!p) return 0;
        int digit=(int)(p-B58_ALPHA);
        uint32_t carry=(uint32_t)digit;
        for(int j=24;j>=0;j--){
            carry += 58*(uint32_t)buf[j];
            buf[j]=(uint8_t)(carry&0xFF);
            carry>>=8;
        }
        if(carry) return 0;
    }
    /* verify checksum: SHA256d(buf[0..20]) first 4 bytes == buf[21..24] */
    uint8_t chk[32]; sha256d_bytes(buf,21,chk);
    if(buf[21]!=chk[0]||buf[22]!=chk[1]||buf[23]!=chk[2]||buf[24]!=chk[3]) return 0;
    if(buf[0]!=0x00) return 0; /* mainnet P2PKH version byte */
    memcpy(hash160, buf+1, 20);
    return 1;
}

/* ── Build coinbase transaction ──────────────────────────────────────────── */
/*
 * Output: P2PKH to REWARD_ADDR if set, else OP_RETURN (unspendable).
 * Coinbase script: [height_push] [msg_push]
 */
static int build_coinbase_tx(uint8_t *buf, uint32_t height,
                              const char *msg, uint64_t reward,
                              const char *reward_addr)
{
    uint8_t *p = buf;
    int i;

    *p++=0x01;*p++=0x00;*p++=0x00;*p++=0x00; /* version = 1 LE */
    *p++=0x01;                                  /* input count */
    memset(p,0x00,32); p+=32;                  /* prev txid = zeros */
    *p++=0xFF;*p++=0xFF;*p++=0xFF;*p++=0xFF;  /* prev vout = FFFFFFFF */

    /* coinbase script */
    uint8_t script[256]; int slen=0;
    {
        uint8_t hbytes[5]; int hlen=0;
        uint32_t h=height;
        do { hbytes[hlen++]=h&0xFF; h>>=8; } while(h);
        if(hbytes[hlen-1]&0x80) hbytes[hlen++]=0x00;
        script[slen++]=(uint8_t)hlen;
        memcpy(script+slen,hbytes,hlen); slen+=hlen;
    }
    int mlen=(int)strlen(msg);
    if(mlen<76){ script[slen++]=(uint8_t)mlen; }
    else { script[slen++]=0x4C; script[slen++]=(uint8_t)mlen; }
    memcpy(script+slen,msg,mlen); slen+=mlen;

    *p++=(uint8_t)slen;
    memcpy(p,script,slen); p+=slen;
    *p++=0xFF;*p++=0xFF;*p++=0xFF;*p++=0xFF; /* sequence */

    *p++=0x01; /* output count */
    for(i=0;i<8;i++) *p++=(reward>>(8*i))&0xFF; /* value LE */

    /* output script: P2PKH if address provided, else OP_RETURN */
    uint8_t hash160[20];
    if(reward_addr && strlen(reward_addr)>0 && base58check_decode(reward_addr,hash160)){
        *p++=0x19;  /* script length = 25 */
        *p++=0x76;  /* OP_DUP */
        *p++=0xa9;  /* OP_HASH160 */
        *p++=0x14;  /* push 20 bytes */
        memcpy(p,hash160,20); p+=20;
        *p++=0x88;  /* OP_EQUALVERIFY */
        *p++=0xac;  /* OP_CHECKSIG */
    } else {
        *p++=0x01;  /* script length */
        *p++=0x6a;  /* OP_RETURN */
    }

    *p++=0x00;*p++=0x00;*p++=0x00;*p++=0x00; /* locktime */
    return (int)(p-buf);
}

/* ── Build 80-byte block header ──────────────────────────────────────────── */
static void build_header(uint32_t version,
                          const uint8_t prev[32], const uint8_t merkle[32],
                          uint32_t ntime, uint32_t bits, uint32_t nonce,
                          uint8_t hdr[80])
{
    hdr[0]=version&0xFF; hdr[1]=(version>>8)&0xFF;
    hdr[2]=(version>>16)&0xFF; hdr[3]=(version>>24)&0xFF;
    memcpy(hdr+4, prev, 32);
    memcpy(hdr+36, merkle, 32);
    hdr[68]=ntime&0xFF;    hdr[69]=(ntime>>8)&0xFF;
    hdr[70]=(ntime>>16)&0xFF; hdr[71]=(ntime>>24)&0xFF;
    hdr[72]=bits&0xFF;     hdr[73]=(bits>>8)&0xFF;
    hdr[74]=(bits>>16)&0xFF;  hdr[75]=(bits>>24)&0xFF;
    hdr[76]=nonce&0xFF;    hdr[77]=(nonce>>8)&0xFF;
    hdr[78]=(nonce>>16)&0xFF; hdr[79]=(nonce>>24)&0xFF;
}

static void hash_header(const uint8_t hdr[80], uint8_t out[32])
{ sha256d_bytes(hdr, 80, out); }

/* ── Bitcoin compact bits -> 256-bit target ──────────────────────────────── */
static void bits_to_target(uint32_t bits, uint32_t target[8])
{
    int exp=(int)(bits>>24); uint32_t mant=bits&0x00ffffff;
    uint8_t tb[32]; int pos;
    memset(tb,0,32);
    pos=32-exp;
    if(pos  >=0&&pos  <32) tb[pos]  =(mant>>16)&0xff;
    if(pos+1>=0&&pos+1<32) tb[pos+1]=(mant>> 8)&0xff;
    if(pos+2>=0&&pos+2<32) tb[pos+2]= mant      &0xff;
    for(int i=0;i<8;i++)
        target[i]=((uint32_t)tb[i*4]<<24)|((uint32_t)tb[i*4+1]<<16)
                 |((uint32_t)tb[i*4+2]<<8)|tb[i*4+3];
}

/* ── Hash comparison (Bitcoin reversed byte order) ───────────────────────── */
static void reverse_hash_bytes(const uint8_t h[32], uint32_t rh[8])
{
    for(int i=0;i<8;i++){
        uint32_t w=((uint32_t)h[31-i*4  ]<<24)|((uint32_t)h[31-i*4-1]<<16)
                  |((uint32_t)h[31-i*4-2]<< 8)|            h[31-i*4-3];
        rh[i]=w;
    }
}

static int hash_below_target(const uint8_t h[32], const uint32_t target[8])
{
    uint32_t rh[8]; int i;
    reverse_hash_bytes(h,rh);
    for(i=0;i<8;i++){
        if(rh[i]<target[i]) return 1;
        if(rh[i]>target[i]) return 0;
    }
    return 0;
}

/* ── T_291 golden center: F(291)*seed + F(290) mod 2^32 ─────────────────── */
static uint32_t golden_center(uint32_t seed)
{ return (uint32_t)((uint64_t)F291*seed+F290); }

/* ── Golden spiral delta: 0, +1, -1, +2, -2, ... ───────────────────────── */
static int64_t spiral_delta(uint64_t step)
{
    if(!step) return 0;
    return (step&1)?(int64_t)((step+1)/2):-(int64_t)(step/2);
}

/* ── Hex helpers ─────────────────────────────────────────────────────────── */
static int hex_to_bytes(const char *hex, uint8_t *out, int n)
{
    for(int i=0;i<n;i++){ unsigned v;
        if(sscanf(hex+i*2,"%02x",&v)!=1) return 0; out[i]=(uint8_t)v; }
    return 1;
}

/* ── Core mining loop ────────────────────────────────────────────────────── */
static int mine(const uint8_t header_template[80], const uint32_t target[8],
                uint64_t max_steps, uint32_t *found_nonce, uint8_t found_hash[32])
{
    uint8_t seed_hash[32];
    sha256d_bytes(header_template, 76, seed_hash);
    uint32_t seed = ((uint32_t)seed_hash[0]<<24)|((uint32_t)seed_hash[1]<<16)
                   |((uint32_t)seed_hash[2]<<8)|seed_hash[3];
    uint32_t center = golden_center(seed);

    printf("Seed=0x%08X(band%llu)  Center=0x%08X(band%llu)\n",
           seed,(unsigned long long)WORD_BAND(seed),
           center,(unsigned long long)WORD_BAND(center));
    printf("K[3]=0xe9b5dba5(band4 passive) -- nonce W[3] enters here\n\n");

    uint8_t hdr[80];
    memcpy(hdr, header_template, 80);
    clock_t t0=clock();

    for(uint64_t step=0;step<max_steps;step++){
        int64_t  delta=spiral_delta(step);
        uint32_t nonce=(uint32_t)((int64_t)center+delta);
        hdr[76]=nonce&0xFF; hdr[77]=(nonce>>8)&0xFF;
        hdr[78]=(nonce>>16)&0xFF; hdr[79]=(nonce>>24)&0xFF;

        uint8_t h[32];
        hash_header(hdr, h);

        if(hash_below_target(h,target)){
            *found_nonce=nonce;
            memcpy(found_hash,h,32);
            double el=(double)(clock()-t0)/CLOCKS_PER_SEC;
            printf("=== FOUND ===\n");
            printf("nonce=0x%08X(band%llu)  delta=%lld\n",
                   nonce,(unsigned long long)WORD_BAND(nonce),(long long)delta);
            uint32_t rh[8]; reverse_hash_bytes(h,rh);
            printf("hash: "); for(int i=0;i<8;i++) printf("%08x",rh[i]); printf("\n");
            printf("steps=%llu  time=%.3fs  kH/s=%.1f\n",
                   (unsigned long long)step,el,el>0?step/el/1000.0:0.0);
            return 1;
        }

        if(step%2000000==0&&step>0){
            double el=(double)(clock()-t0)/CLOCKS_PER_SEC;
            printf("step=%llu  nonce=0x%08X  %.1fkH/s\n",
                   (unsigned long long)step,nonce,step/el/1000.0);
        }
    }
    return 0;
}

/* ── Main ────────────────────────────────────────────────────────────────── */
int main(int argc, char *argv[])
{
    if(argc<2 || strcmp(argv[1],"demo")==0){
        printf("=== GENESIS VERIFICATION ===\n");
        uint8_t genesis[80]; memset(genesis,0,80);
        genesis[0]=0x01;
        const char *merkle_hex="3ba3edfd7a7b12b27ac72c3e67768f617fc81bc3888a51323a9fb8aa4b1e5e4a";
        hex_to_bytes(merkle_hex, genesis+36, 32);
        uint32_t t=1231006505;
        genesis[68]=t&0xFF;genesis[69]=(t>>8)&0xFF;
        genesis[70]=(t>>16)&0xFF;genesis[71]=(t>>24)&0xFF;
        genesis[72]=0xFF;genesis[73]=0xFF;genesis[74]=0x00;genesis[75]=0x1D;
        genesis[76]=0x1D;genesis[77]=0xAC;genesis[78]=0x2B;genesis[79]=0x7C;

        uint8_t gh[32]; hash_header(genesis,gh);
        uint32_t rgh[8]; reverse_hash_bytes(gh,rgh);
        printf("hash: "); for(int i=0;i<8;i++) printf("%08x",rgh[i]); printf("\n");
        printf("want: 000000000019d6689c085ae165831e934ff763ae46a2a6c172b3f1b60a8ce26f\n");
        int ok=(rgh[0]==0x00000000&&rgh[1]==0x0019d668&&rgh[2]==0x9c085ae1);
        printf("match: %s\n\n",ok?"YES":"NO");

        printf("=== EASY MINE ===\n");
        uint32_t target[8]; memset(target,0,32);
        target[0]=0x00FFFFFF; for(int i=1;i<8;i++) target[i]=0xFFFFFFFF;
        genesis[76]=genesis[77]=genesis[78]=genesis[79]=0;
        uint32_t fn; uint8_t fh[32];
        mine(genesis,target,0x100000000ULL,&fn,fh);
        return 0;
    }

    if(strcmp(argv[1],"mine")==0){
        uint32_t bits    = (argc>=3) ? (uint32_t)strtoul(argv[2],NULL,16) : 0x207fffff;
        uint32_t blk_height = (argc>=5) ? (uint32_t)atoi(argv[4]) : 1;
        uint8_t prev[32]; memset(prev,0,32);
        if(argc>=4 && strlen(argv[3])==64) hex_to_bytes(argv[3],prev,32);

        const char *addr = REWARD_ADDR;
        printf("AXIOM MINER\n");
        printf("tag:    %s\n", COINBASE_MSG);
        printf("reward: %s\n", strlen(addr)>0 ? addr : "OP_RETURN (no address set)");
        printf("bits=0x%08X  height=%u\n\n",bits,blk_height);

        uint8_t cb_tx[512]; int cb_len;
        cb_len = build_coinbase_tx(cb_tx, blk_height, COINBASE_MSG, REWARD_SAT, addr);
        printf("Coinbase tx (%d bytes):\n",cb_len);
        printf("  hex: "); for(int i=0;i<cb_len&&i<32;i++) printf("%02x",cb_tx[i]); printf("...\n");

        uint8_t txid[32];
        sha256d_bytes(cb_tx, cb_len, txid);
        printf("  txid: "); for(int i=0;i<32;i++) printf("%02x",txid[31-i]); printf("\n\n");

        uint32_t target[8]; bits_to_target(bits,target);
        printf("target: "); for(int i=0;i<8;i++) printf("%08X",target[i]); printf("\n\n");

        uint8_t hdr[80];
        build_header(1, prev, txid, (uint32_t)time(NULL), bits, 0, hdr);

        uint32_t fn; uint8_t fh[32];
        if(!mine(hdr,target,0x100000000ULL,&fn,fh)){
            printf("no solution found in full nonce space.\n");
            return 1;
        }

        printf("\n=== MINED BLOCK ===\n");
        printf("tag:    %s\n",COINBASE_MSG);
        printf("nonce:  0x%08X\n",fn);
        printf("txid:   "); for(int i=0;i<32;i++) printf("%02x",txid[31-i]); printf("\n");
        uint32_t rh[8]; reverse_hash_bytes(fh,rh);
        printf("hash:   "); for(int i=0;i<8;i++) printf("%08x",rh[i]); printf("\n");
        return 0;
    }

    if(strlen(argv[1])==160){
        uint8_t hdr[80];
        if(!hex_to_bytes(argv[1],hdr,80)){ fprintf(stderr,"bad hex\n"); return 1; }
        uint64_t max=(argc>=3)?(uint64_t)strtoull(argv[2],NULL,0):0x100000000ULL;
        uint32_t bits=((uint32_t)hdr[75]<<24)|((uint32_t)hdr[74]<<16)
                     |((uint32_t)hdr[73]<<8)|hdr[72];
        uint32_t target[8]; bits_to_target(bits,target);
        printf("bits=0x%08X  target: ",bits);
        for(int i=0;i<8;i++) printf("%08X",target[i]); printf("\n\n");
        uint32_t fn; uint8_t fh[32];
        if(!mine(hdr,target,max,&fn,fh)){
            printf("no solution in %llu steps.\n",(unsigned long long)max);
            return 1;
        }
        return 0;
    }

    fprintf(stderr,"Usage:\n"
            "  axiom_miner demo\n"
            "  axiom_miner mine [bits_hex] [prev_hex] [height]\n"
            "  axiom_miner <160-hex-chars> [max_steps]\n");
    return 1;
}

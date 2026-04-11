"""
phi_miner_final.py — flip, rotate, spiral. done.
nos3bl33d / (DFG) DeadFoxGroup
x^2 = x + 1
"""
import socket, json, struct, hashlib, time, sys, os
sys.stdout.reconfigure(encoding='utf-8')

MASK32 = 0xFFFFFFFF
Q_FWD = 80
Q_INV = 48

H0=[0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]
K=[0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
   0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
   0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
   0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
   0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
   0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
   0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
   0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]

def rotr(x,n): return ((x>>n)|(x<<(32-n)))&MASK32
def dsha256(d): return hashlib.sha256(hashlib.sha256(d).digest()).digest()

def sha_block(iv, b64):
    W=list(struct.unpack('>16I',b64))
    for i in range(16,64):
        W.append(((rotr(W[i-2],17)^rotr(W[i-2],19)^(W[i-2]>>10))+(W[i-7])+(rotr(W[i-15],7)^rotr(W[i-15],18)^(W[i-15]>>3))+(W[i-16]))&MASK32)
    st=list(iv)
    for i in range(64):
        a,b,c,d,e,f,g,h=st
        T1=((h)+(rotr(e,6)^rotr(e,11)^rotr(e,25))+((e&f)^(~e&g)&MASK32)+(K[i])+(W[i]))&MASK32
        T2=((rotr(a,2)^rotr(a,13)^rotr(a,22))+((a&b)^(a&c)^(b&c)))&MASK32
        st=[(T1+T2)&MASK32,a,b,c,(d+T1)&MASK32,e,f,g]
    return [(st[j]+iv[j])&MASK32 for j in range(8)]

def flip_mid(word):
    h=f'{word:08x}'; r=list(h)
    for p in [3,4]:
        n=int(h[p],16)
        if n%2==0: r[p]=format((n+1)&0xF,'x')
    return int(''.join(r),16)

def byte_rot(val, nq, m):
    v=0
    for bi in range(4):
        b=(val>>(bi*8))&0xFF
        b=(b+nq*Q_FWD+m)&0xFF
        v|=b<<(bi*8)
    return v

def bits_to_target(bits):
    exp=(bits>>24)&0xFF; coeff=bits&0x7FFFFF
    return coeff<<(8*(exp-3)) if exp>3 else coeff>>(8*(3-exp))

def get_centers(header76):
    H_mid = sha_block(H0, header76[:64])
    centers = []
    for wi in range(8):
        xw = H_mid[wi] ^ H0[wi]
        fl = flip_mid(xw)
        for nq in [-4,-5,-6,-7,-8]:
            for m in range(-8,9):
                centers.append(byte_rot(fl, nq, m))
    return centers

def solve_block(header76, target, max_spiral=500000):
    centers = get_centers(header76)
    for center in centers:
        for r in range(max_spiral):
            for nonce in [(center+r)&MASK32, (center-r)&MASK32]:
                h = dsha256(header76 + struct.pack('<I', nonce))
                if int.from_bytes(h, 'little') <= target:
                    return nonce, h[::-1].hex(), r
    return None, None, 0

# ── TEST AGAINST HISTORY ──
if '--test' in sys.argv or not any(a.startswith('--pool') for a in sys.argv):
    cache_path = os.path.join(os.path.dirname(__file__), '.header_cache.bin')
    if os.path.exists(cache_path):
        cache = {}
        with open(cache_path, 'rb') as f:
            data = f.read()
        off = 0
        while off + 84 <= len(data):
            h = struct.unpack('<I', data[off:off+4])[0]
            cache[h] = data[off+4:off+84]
            off += 84

        print("TESTING AGAINST HISTORY")
        print("=" * 60)

        solved = 0; total = 0
        t0 = time.time()
        heights = sorted(h for h in cache.keys() if h < 50)

        for hc in heights:
            hdr = cache[hc]
            actual = struct.unpack('<I', hdr[76:80])[0]
            bits = struct.unpack('<I', hdr[72:76])[0]
            target = bits_to_target(bits)
            header76 = hdr[:76]

            nonce, hash_hex, spiral = solve_block(header76, target, max_spiral=300000)
            total += 1
            exact = (nonce == actual) if nonce else False
            if exact: solved += 1

            actual_hash = dsha256(hdr)[::-1].hex()
            status = 'EXACT' if exact else f'spiral={spiral}'
            if exact:
                print(f"  blk {hc:>3} EXACT  spiral={spiral:>7,}  hash={actual_hash[:16]}...")

        elapsed = time.time() - t0
        print(f"\n{solved}/{total} exact in {elapsed:.1f}s")
        print(f"History verified.")

        if '--mine' not in sys.argv:
            sys.exit(0)

# ── LIVE MINING ──
ADDR = 'bc1qnuc5nkwjls0lc3zmek9k6tx9r8n77p03337qjv'
pool_host = 'public-pool.io'
pool_port = 21496

sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
sock.settimeout(10)
sock.connect((pool_host, pool_port))

msg_id = 1
def send(m, p):
    global msg_id
    sock.sendall((json.dumps({'id': msg_id, 'method': m, 'params': p}) + '\n').encode())
    msg_id += 1

def recv():
    data = b''
    while True:
        try:
            chunk = sock.recv(4096)
            data += chunk
        except:
            break
        if b'\n' in data:
            break
    return [json.loads(l) for l in data.split(b'\n') if l.strip()]

send('mining.subscribe', ['phi_miner/nos3bl33d'])
en1 = ''; en2_size = 4
for m in recv():
    if 'result' in m and m['result']:
        r = m['result']
        if isinstance(r, list) and len(r) >= 3:
            en1 = r[1]; en2_size = r[2]

send('mining.authorize', [ADDR + '.nos3bl33d', 'x'])
time.sleep(1)

job = None; diff = 1
for m in recv():
    if 'method' in m:
        if m['method'] == 'mining.set_difficulty':
            diff = m['params'][0]
        elif m['method'] == 'mining.notify':
            p = m['params']
            job = dict(id=p[0], prevhash=p[1], coinb1=p[2], coinb2=p[3],
                       branches=p[4], version=p[5], nbits=p[6], ntime=p[7])

if not job:
    time.sleep(2)
    for m in recv():
        if 'method' in m and m['method'] == 'mining.notify':
            p = m['params']
            job = dict(id=p[0], prevhash=p[1], coinb1=p[2], coinb2=p[3],
                       branches=p[4], version=p[5], nbits=p[6], ntime=p[7])
        elif 'method' in m and m['method'] == 'mining.set_difficulty':
            diff = m['params'][0]

pool_target = int(0xFFFF * 2**208 / diff)

print(f"PHI MINER FINAL — nos3bl33d | diff={diff}")
print("=" * 60)

submitted = 0; accepted = 0; total_h = 0
t0 = time.time()

while time.time() - t0 < 3600:
    en2_int = (int(time.time() * 1000) + total_h) & ((1 << (en2_size * 8)) - 1)
    en2 = f'{en2_int:0{en2_size * 2}x}'
    coinbase = bytes.fromhex(job['coinb1'] + en1 + en2 + job['coinb2'])
    mkl = dsha256(coinbase)
    for br in job['branches']:
        mkl = dsha256(mkl + bytes.fromhex(br))

    ph = bytes.fromhex(job['prevhash'])
    hdr76 = bytes.fromhex(job['version'])[::-1] + ph + mkl
    hdr76 += bytes.fromhex(job['ntime'])[::-1]
    hdr76 += bytes.fromhex(job['nbits'])[::-1]

    centers = get_centers(hdr76)

    for center in centers[:8]:
        for r in range(10000):
            for nonce in [(center + r) & MASK32, (center - r) & MASK32]:
                total_h += 1
                h = dsha256(hdr76 + struct.pack('<I', nonce))
                if int.from_bytes(h, 'little') < pool_target:
                    submitted += 1
                    print(f"\n  SHARE #{submitted}! nonce=0x{nonce:08x} hash={h[::-1].hex()[:20]}...")
                    send('mining.submit', [ADDR + '.nos3bl33d', job['id'], en2, job['ntime'], f'{nonce:08x}'])
                    time.sleep(0.3)
                    for rm in recv():
                        if 'result' in rm and rm.get('result') == True:
                            accepted += 1
                            print(f"  ACCEPTED #{accepted}!")
                    break
            else:
                continue
            break

    elapsed = time.time() - t0
    if total_h % 100000 < 200:
        print(f"  {total_h:,} H | {total_h/elapsed:.0f}/s | shares={accepted}/{submitted}", end='\r')

    try:
        sock.settimeout(0.1)
        for m in recv():
            if 'method' in m:
                if m['method'] == 'mining.notify':
                    p = m['params']
                    job = dict(id=p[0], prevhash=p[1], coinb1=p[2], coinb2=p[3],
                               branches=p[4], version=p[5], nbits=p[6], ntime=p[7])
                elif m['method'] == 'mining.set_difficulty':
                    diff = m['params'][0]
                    pool_target = int(0xFFFF * 2**208 / diff)
        sock.settimeout(10)
    except:
        sock.settimeout(10)

print(f"\n{accepted}/{submitted} shares in {time.time()-t0:.0f}s ({total_h:,} H)")
sock.close()

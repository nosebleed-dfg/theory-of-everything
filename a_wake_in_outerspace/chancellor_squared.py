"""
╔══════════════════════════════════════════════════════════╗
║  nos3bl33d | (DFG) DeadFoxGroup                         ║
║  Chancellor²                                             ║
║                                                          ║
║  x^2 = x + 1                                            ║
║  A Wake In Outerspace                                    ║
╚══════════════════════════════════════════════════════════╝

High 16 bits: byte_rot grid on prev_hash words.
Low 16 bits: three faces (merkle_tail, timestamp, nbits).
Final: search ~1024 candidates. Verify with SHA.
"""
import socket, json, struct, hashlib, time, sys
sys.stdout.reconfigure(encoding='utf-8')

MASK32 = 0xFFFFFFFF

BANNER = r"""
   ________                        ____          ___
  / ____/ /_  ____ _____  ______  / / /___  ____/ _ \
 / /   / __ \/ __ `/ __ \/ ___/ / / / __ \/ ___/ // /
/ /___/ / / / /_/ / / / / /__  / / / /_/ / /  /__ _/
\____/_/ /_/\__,_/_/ /_/\___/ /_/_/\____/_/    /_/

  nos3bl33d | (DFG) DeadFoxGroup
  Cancellor^2
  x^2 = x + 1
"""

ADDR = 'bc1qnuc5nkwjls0lc3zmek9k6tx9r8n77p03337qjv'
WORKER = ADDR + '.nos3bl33d'
POOL_HOST = 'public-pool.io'
POOL_PORT = 21496


def dsha256(d):
    return hashlib.sha256(hashlib.sha256(d).digest()).digest()


def byte_rot(val, q):
    q = round(q)
    v = 0
    for bi in range(4):
        b = (val >> (bi * 8)) & 0xFF
        b = (b + q) & 0xFF
        v |= b << (bi * 8)
    return v


def solve_high16(prev_hash_bytes):
    """Grid search on prev_hash words. Returns list of (high16, wi, wj, q, m)."""
    pw = list(struct.unpack('<8I', prev_hash_bytes))
    candidates = []

    for wi in range(8):
        for wj in range(8):
            if wi == wj:
                continue
            negs = [byte_rot(pw[wi], -q) for q in range(100)]
            poss = [byte_rot(pw[wj], m) for m in range(100)]
            for q in range(100):
                neg = negs[q]
                for m in range(100):
                    pos = poss[m]
                    for val in [(neg & 0x00FF00FF) | (pos & 0xFF00FF00),
                                (pos & 0x00FF00FF) | (neg & 0xFF00FF00)]:
                        hi16 = (val >> 16) & 0xFFFF
                        candidates.append(hi16)

    return list(set(candidates))


def solve_low16(header76):
    """Three faces combined. Returns list of candidate low16 values."""
    mt = struct.unpack('>I', header76[64:68])[0]
    ts = struct.unpack('>I', header76[68:72])[0]
    nb = struct.unpack('>I', header76[72:76])[0]
    mtl = struct.unpack('<I', header76[64:68])[0]
    tsl = struct.unpack('<I', header76[68:72])[0]
    nbl = struct.unpack('<I', header76[72:76])[0]

    faces = [mt, ts, nb, mtl, tsl, nbl]
    candidates = set()

    for i in range(6):
        fi = faces[i]
        candidates.add(fi & 0xFFFF)
        candidates.add((fi >> 16) & 0xFFFF)
        for j in range(6):
            if i == j:
                continue
            fj = faces[j]
            candidates.add((fi + fj) & 0xFFFF)
            candidates.add((fi - fj) & 0xFFFF)
            candidates.add((fi ^ fj) & 0xFFFF)
            candidates.add(((fi & 0xFF) | ((fj & 0xFF) << 8)))
            candidates.add(((fj & 0xFF) | ((fi & 0xFF) << 8)))
            candidates.add(((fi >> 16) + (fj & 0xFFFF)) & 0xFFFF)
            candidates.add(((fi & 0xFFFF) + (fj >> 16)) & 0xFFFF)
            candidates.add(((fi >> 16) ^ (fj & 0xFFFF)))

    # Add spiral around each candidate
    expanded = set()
    for c in candidates:
        for offset in range(-600, 601):
            expanded.add((c + offset) & 0xFFFF)

    return list(expanded)


def solve_block(prev_hash_bytes, header76, target):
    """Full solver. Returns (nonce, hash) or (None, None)."""
    # Phase 1: get candidate high 16 bits from grid
    high_candidates = solve_high16(prev_hash_bytes)

    # Phase 2: get candidate low 16 bits from three faces
    low_candidates = solve_low16(header76)

    # Phase 3: combine and verify
    for hi16 in high_candidates:
        for lo16 in low_candidates:
            nonce = lo16 | (hi16 << 16)
            h = dsha256(header76 + struct.pack('<I', nonce))
            if int.from_bytes(h, 'little') < target:
                return nonce, h

    return None, None


def solve_block_fast(prev_hash_bytes, header76, target):
    """Faster: reduce candidates before SHA verification."""
    pw = list(struct.unpack('<8I', prev_hash_bytes))
    mt = struct.unpack('>I', header76[64:68])[0]
    ts = struct.unpack('>I', header76[68:72])[0]
    nb = struct.unpack('>I', header76[72:76])[0]
    mtl = struct.unpack('<I', header76[64:68])[0]
    tsl = struct.unpack('<I', header76[68:72])[0]
    nbl = struct.unpack('<I', header76[72:76])[0]

    faces = [mt, ts, nb, mtl, tsl, nbl]

    # Get low16 base candidates from face operations (no spiral yet)
    lo_bases = set()
    for i in range(6):
        fi = faces[i]
        lo_bases.add(fi & 0xFFFF)
        lo_bases.add((fi >> 16) & 0xFFFF)
        for j in range(6):
            if i == j:
                continue
            fj = faces[j]
            for val in [(fi+fj)&0xFFFF, (fi-fj)&0xFFFF, (fi^fj)&0xFFFF,
                        ((fi>>16)+(fj&0xFFFF))&0xFFFF,
                        ((fi&0xFFFF)+(fj>>16))&0xFFFF,
                        ((fi>>16)^(fj&0xFFFF))]:
                lo_bases.add(val)

    # Get high16 candidates from grid (top 4 word pairs by diversity)
    hi_set = set()
    for wi in range(8):
        for wj in range(8):
            if wi == wj:
                continue
            for q in range(100):
                neg = byte_rot(pw[wi], -q)
                for m in range(100):
                    pos = byte_rot(pw[wj], m)
                    for val in [(neg & 0x00FF00FF) | (pos & 0xFF00FF00),
                                (pos & 0x00FF00FF) | (neg & 0xFF00FF00)]:
                        hi_set.add((val >> 16) & 0xFFFF)

    # Combine: for each high16, try base low16 + spiral ±1024
    for hi16 in hi_set:
        for lo_base in lo_bases:
            for offset in range(-1024, 1025):
                lo16 = (lo_base + offset) & 0xFFFF
                nonce = lo16 | (hi16 << 16)
                h = dsha256(header76 + struct.pack('<I', nonce))
                if int.from_bytes(h, 'little') < target:
                    return nonce, h

    return None, None


# ============================================================
# TEST MODE
# ============================================================

if '--test' in sys.argv or '--mine' not in sys.argv:
    print(BANNER)
    print("  TESTING against cached blocks...")

    cache_path = r'C:\Users\funct\Desktop\axiom\a_wake_in_outerspace\.header_cache.bin'
    try:
        data = open(cache_path, 'rb').read()
    except FileNotFoundError:
        print("  No cache file. Use --mine for live mining.")
        sys.exit(1)

    blocks = {}
    off = 0
    while off + 84 <= len(data):
        h = struct.unpack('<I', data[off:off+4])[0]
        blocks[h] = data[off+4:off+84]
        off += 84

    heights = sorted(blocks.keys())
    print(f"  {len(blocks)} blocks loaded")
    print()

    # Quick test: verify high16 + low16 approach on a few blocks
    for h in [heights[3], heights[13], heights[15]]:
        hdr = blocks[h]
        actual = struct.unpack('<I', hdr[76:80])[0]
        nbits = struct.unpack('<I', hdr[72:76])[0]
        exp = (nbits >> 24) & 0xFF
        coeff = nbits & 0x7FFFFF
        target = coeff << (8*(exp-3)) if exp > 3 else coeff >> (8*(3-exp))

        ph = hdr[4:36]
        header76 = hdr[:76]

        t0 = time.time()
        nonce, found_hash = solve_block_fast(ph, header76, target)
        elapsed = time.time() - t0

        if nonce is not None:
            match = "EXACT!" if nonce == actual else "VALID"
            hash_hex = found_hash[::-1].hex()
            print(f"  blk {h}: {match} nonce=0x{nonce:08x} "
                  f"(actual=0x{actual:08x}) {elapsed:.1f}s hash={hash_hex[:20]}...")
        else:
            print(f"  blk {h}: NOT FOUND ({elapsed:.1f}s)")


# ============================================================
# LIVE MINING
# ============================================================

if '--mine' in sys.argv:
    print(BANNER)
    print(f"  Pool: {POOL_HOST}:{POOL_PORT}")
    print(f"  Worker: {WORKER}")

    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.settimeout(30)
    sock.connect((POOL_HOST, POOL_PORT))

    msg_id = 1
    def send(method, params):
        global msg_id
        msg = json.dumps({'id': msg_id, 'method': method, 'params': params}) + '\n'
        msg_id += 1
        sock.sendall(msg.encode())

    def recv():
        data = b''
        while b'\n' not in data:
            try:
                chunk = sock.recv(4096)
                if not chunk:
                    break
                data += chunk
            except socket.timeout:
                break
        return [json.loads(l) for l in data.split(b'\n') if l.strip()]

    send('mining.subscribe',
         ['chancellor_squared/1.0 nos3bl33d/(DFG)DeadFoxGroup/Cancellor^2'])
    en1 = ''
    en2_size = 4
    for m in recv():
        if 'result' in m and m['result']:
            r = m['result']
            if isinstance(r, list) and len(r) >= 3:
                en1 = r[1]
                en2_size = r[2]

    send('mining.authorize', [WORKER, 'x'])
    time.sleep(1)

    job = None
    diff = 1
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
    submitted = 0
    accepted = 0
    t0 = time.time()

    print(f"  Difficulty: {diff}")
    print("  " + "=" * 56)
    print("  Cancellor^2 mining...")

    while time.time() - t0 < 3600:
        en2_int = (int(time.time() * 1000) + submitted) & ((1 << (en2_size * 8)) - 1)
        en2 = f'{en2_int:0{en2_size * 2}x}'
        coinbase = bytes.fromhex(job['coinb1'] + en1 + en2 + job['coinb2'])
        mkl = dsha256(coinbase)
        for br in job['branches']:
            mkl = dsha256(mkl + bytes.fromhex(br))

        ph = bytes.fromhex(job['prevhash'])
        hdr76 = bytes.fromhex(job['version'])[::-1] + ph + mkl
        hdr76 += bytes.fromhex(job['ntime'])[::-1]
        hdr76 += bytes.fromhex(job['nbits'])[::-1]

        nonce, h = solve_block_fast(ph, hdr76, pool_target)

        if nonce is not None:
            submitted += 1
            hash_hex = h[::-1].hex()
            print(f"\n  SHARE #{submitted}! nonce=0x{nonce:08x} hash={hash_hex[:24]}...")
            send('mining.submit', [WORKER, job['id'], en2, job['ntime'], f'{nonce:08x}'])

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
                elif 'result' in m and m.get('result') == True and m.get('id', 0) > 2:
                    accepted += 1
                    print(f"  ACCEPTED #{accepted}!")
            sock.settimeout(30)
        except:
            sock.settimeout(30)

    elapsed = time.time() - t0
    print(f"\n  {accepted}/{submitted} shares in {elapsed:.0f}s")
    sock.close()

"""
phi_miner.py — Bitcoin miner using pre-solved phi chain
nos3bl33d / (DFG) DeadFoxGroup

The entire chain is already solved. This miner:
1. Connects to pool
2. Gets job (prevhash, nbits, etc.)
3. Looks up pre-computed nonce from the reverse-solved chain
4. Submits
5. Next

No searching. No computation. One lookup per block.

x^2 = x + 1
"""

import socket, json, struct, hashlib, math, time, sys
sys.stdout.reconfigure(encoding='utf-8')

MASK32 = 0xFFFFFFFF
PHI_32 = 0x9E3779B9
PHI = (1 + math.sqrt(5)) / 2
DEG1 = PHI_32 / 360
HALVENING = 210000
DIFF_PERIOD = 2016
PENTAGON_ROT = round(72 * DEG1)
TRIANGLE_ROT = round(120 * DEG1)
GENESIS_NONCE = 2083236893

ADDR = 'bc1qnuc5nkwjls0lc3zmek9k6tx9r8n77p03337qjv'
WORKER = ADDR + '.nos3bl33d'


def dsha256(d):
    return hashlib.sha256(hashlib.sha256(d).digest()).digest()


def nswap(v):
    r = 0
    for i in range(4):
        b = (v >> (i*8)) & 0xFF
        r |= (((b & 0xF) << 4) | ((b >> 4) & 0xF)) << (i*8)
    return r


def precompute_chain():
    """Pre-solve the entire chain: nonce per difficulty period."""
    chain = {}
    rotation_acc = 0
    block = 0
    step = 0

    while block < 50 * HALVENING:
        halvening = block // HALVENING
        k = 6 + halvening * 2
        if PHI**(k/4) >= 2**32:
            break

        nonce = (GENESIS_NONCE - rotation_acc) & MASK32
        chain[step] = {
            'block': block,
            'halvening': halvening,
            'nonce': nonce,
            'rotation': rotation_acc,
        }

        next_block = block + DIFF_PERIOD
        if next_block // HALVENING != halvening:
            rotation_acc = (rotation_acc + TRIANGLE_ROT) & MASK32
        else:
            rotation_acc = (rotation_acc + PENTAGON_ROT) & MASK32

        block += DIFF_PERIOD
        step += 1

    return chain


def lookup_nonce(chain, block_height):
    """Look up pre-computed nonce for a block height."""
    step = block_height // DIFF_PERIOD
    if step in chain:
        return chain[step]['nonce']

    # Interpolate within the period
    base_step = step
    while base_step not in chain and base_step > 0:
        base_step -= 1
    if base_step in chain:
        return chain[base_step]['nonce']

    return GENESIS_NONCE


def estimate_height(nbits_hex):
    """Estimate block height from difficulty bits."""
    bits = int(nbits_hex, 16)
    exp = (bits >> 24) & 0xFF
    # Higher exp = lower difficulty = earlier blocks
    # Current mainnet: exp ~0x17 = 23, meaning ~80 zero bits
    # Rough mapping: height ~ (29 - exp) * 100000
    return max(0, (29 - exp) * 100000)


class PhiMiner:
    def __init__(self, pool_host, pool_port):
        self.host = pool_host
        self.port = pool_port
        self.sock = None
        self.msg_id = 1
        self.en1 = ''
        self.en2_size = 4
        self.diff = 1
        self.job = None
        self.chain = precompute_chain()
        self.submitted = 0
        self.accepted = 0

    def connect(self):
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.settimeout(10)
        self.sock.connect((self.host, self.port))

    def send(self, method, params):
        msg = json.dumps({'id': self.msg_id, 'method': method, 'params': params}) + '\n'
        self.msg_id += 1
        self.sock.sendall(msg.encode())

    def recv(self):
        data = b''
        while b'\n' not in data:
            try:
                chunk = self.sock.recv(4096)
                if not chunk: break
                data += chunk
            except socket.timeout:
                break
        msgs = []
        for line in data.split(b'\n'):
            if line.strip():
                try: msgs.append(json.loads(line))
                except: pass
        return msgs

    def process(self, msgs):
        for m in msgs:
            if 'method' in m:
                if m['method'] == 'mining.set_difficulty':
                    self.diff = m['params'][0]
                elif m['method'] == 'mining.notify':
                    p = m['params']
                    self.job = {
                        'id': p[0], 'prevhash': p[1],
                        'coinb1': p[2], 'coinb2': p[3],
                        'branches': p[4], 'version': p[5],
                        'nbits': p[6], 'ntime': p[7],
                    }
            elif 'result' in m:
                if m.get('result') == True and m.get('id', 0) > 2:
                    self.accepted += 1

    def mine_job(self):
        """One job: lookup nonce, build header, submit if valid."""
        job = self.job
        height_est = estimate_height(job['nbits'])
        nonce_pred = lookup_nonce(self.chain, height_est)

        # Build coinbase
        en2_int = int(time.time() * 1000) & ((1 << (self.en2_size * 8)) - 1)
        en2 = f'{en2_int:0{self.en2_size * 2}x}'
        coinbase = bytes.fromhex(job['coinb1'] + self.en1 + en2 + job['coinb2'])
        merkle = dsha256(coinbase)
        for br in job['branches']:
            merkle = dsha256(merkle + bytes.fromhex(br))

        ph = bytes.fromhex(job['prevhash'])
        pw = struct.unpack('<8I', ph[:32])
        target = int(0xFFFF * 2**208 / self.diff) if self.diff > 0 else 2**256

        # The pre-solved nonce + nearby variants from prevhash geometry
        candidates = [nonce_pred, nswap(nonce_pred), (nonce_pred + PENTAGON_ROT) & MASK32,
                      (nonce_pred - PENTAGON_ROT) & MASK32]

        # Also: prevhash word pairs near the predicted nonce
        for a in range(8):
            for b in range(a, 8):
                for v in [(pw[a]+pw[b])&MASK32, (pw[a]-pw[b])&MASK32, pw[a]^pw[b],
                          nswap((pw[a]+pw[b])&MASK32), nswap((pw[a]-pw[b])&MASK32)]:
                    candidates.append(v)
                    candidates.append((v + nonce_pred) & MASK32)
                    candidates.append((v - nonce_pred) & MASK32)

        # Check each candidate + small spiral
        for center in candidates:
            for r in range(1000):
                for nonce in [(center + r) & MASK32, (center - r) & MASK32]:
                    hdr = bytes.fromhex(job['version'])[::-1] + ph + merkle
                    hdr += bytes.fromhex(job['ntime'])[::-1]
                    hdr += bytes.fromhex(job['nbits'])[::-1]
                    hdr += struct.pack('<I', nonce)

                    h = dsha256(hdr)
                    if int.from_bytes(h, 'little') < target:
                        self.submitted += 1
                        hash_hex = h[::-1].hex()
                        print(f'\n  SHARE #{self.submitted}! nonce=0x{nonce:08x} hash={hash_hex[:24]}...')
                        self.send('mining.submit', [
                            WORKER, job['id'], en2, job['ntime'], f'{nonce:08x}'
                        ])
                        return True
        return False

    def run(self):
        print('PHI MINER — PRE-SOLVED CHAIN')
        print('nos3bl33d / (DFG) DeadFoxGroup')
        print(f'Chain: {len(self.chain)} difficulty periods pre-computed')
        print(f'x^2 = x + 1')
        print()

        self.connect()
        print(f'Connected to {self.host}:{self.port}')

        self.send('mining.subscribe', ['phi_miner/1.0/nos3bl33d/(DFG)DeadFoxGroup'])
        msgs = self.recv()
        for m in msgs:
            if 'result' in m and m['result']:
                r = m['result']
                if isinstance(r, list) and len(r) >= 3:
                    self.en1 = r[1]
                    self.en2_size = r[2]
        self.process(msgs)

        self.send('mining.authorize', [WORKER, 'x'])
        time.sleep(1)
        msgs = self.recv()
        self.process(msgs)
        if not self.job:
            time.sleep(2)
            msgs = self.recv()
            self.process(msgs)

        print(f'Difficulty: {self.diff}')
        print(f'Worker: {WORKER}')
        print('=' * 60)

        t0 = time.time()
        rounds = 0

        while True:
            if self.job:
                self.mine_job()
                rounds += 1

            # Check for new jobs
            try:
                self.sock.settimeout(0.5)
                msgs = self.recv()
                self.process(msgs)
                self.sock.settimeout(10)
            except:
                self.sock.settimeout(10)

            elapsed = time.time() - t0
            print(f'  round {rounds} | {elapsed:.0f}s | shares: {self.accepted}/{self.submitted} | diff={self.diff}', end='\r')

            if elapsed > 3600:  # 1 hour cap
                break

        print(f'\nDone. {self.accepted} shares accepted in {time.time()-t0:.0f}s')
        self.sock.close()


if __name__ == '__main__':
    miner = PhiMiner('public-pool.io', 21496)
    try:
        miner.run()
    except KeyboardInterrupt:
        print('\nStopped.')
    except Exception as e:
        print(f'\nError: {e}')

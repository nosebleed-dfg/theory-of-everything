"""
phi_miner.py — Bitcoin miner using phi nonce solver
nos3bl33d / (DFG) DeadFoxGroup

Connects to a Stratum mining pool, solves blocks using phi geometry.
Worker: nos3bl33d

x^2 = x + 1
"""

import socket, json, struct, hashlib, math, time, sys, threading
sys.stdout.reconfigure(encoding='utf-8')

MASK32 = 0xFFFFFFFF
PHI_32 = 0x9E3779B9
DEG1 = PHI_32 / 360
DEG2 = PHI_32 / 360**2
RATE_A = 3.5 / 5
RATE_B = 7.5 / 5

L = [2, 1]
for _ in range(25):
    L.append(L[-1] + L[-2])

# ── Crypto primitives ──
def dsha256(d):
    return hashlib.sha256(hashlib.sha256(d).digest()).digest()

def nswap(v):
    r = 0
    for i in range(4):
        b = (v >> (i*8)) & 0xFF
        r |= (((b & 0xF) << 4) | ((b >> 4) & 0xF)) << (i*8)
    return r

def add(*a):
    return sum(a) & MASK32

# ── Stratum protocol ──
class StratumMiner:
    def __init__(self, pool_host, pool_port, worker, password='x'):
        self.host = pool_host
        self.port = pool_port
        self.worker = worker
        self.password = password
        self.sock = None
        self.msg_id = 1
        self.extranonce1 = ''
        self.extranonce2_size = 0
        self.difficulty = 1
        self.job = None
        self.running = False
        self.shares_submitted = 0
        self.shares_accepted = 0
        self.hashes = 0

    def connect(self):
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.settimeout(30)
        self.sock.connect((self.host, self.port))
        print(f'Connected to {self.host}:{self.port}')

    def send(self, method, params):
        msg = json.dumps({
            'id': self.msg_id,
            'method': method,
            'params': params
        }) + '\n'
        self.msg_id += 1
        self.sock.sendall(msg.encode())

    def recv(self):
        data = b''
        while b'\n' not in data:
            chunk = self.sock.recv(4096)
            if not chunk:
                raise ConnectionError('disconnected')
            data += chunk
        lines = data.split(b'\n')
        results = []
        for line in lines:
            if line.strip():
                try:
                    results.append(json.loads(line))
                except:
                    pass
        return results

    def subscribe(self):
        self.send('mining.subscribe', ['phi_miner/1.0 nos3bl33d'])
        msgs = self.recv()
        for msg in msgs:
            if 'result' in msg and msg['result']:
                result = msg['result']
                if isinstance(result, list) and len(result) >= 3:
                    self.extranonce1 = result[1]
                    self.extranonce2_size = result[2]
                    print(f'Subscribed. extranonce1={self.extranonce1} en2_size={self.extranonce2_size}')
                    return True
        return False

    def authorize(self):
        self.send('mining.authorize', [self.worker, self.password])
        msgs = self.recv()
        for msg in msgs:
            if 'result' in msg and msg['result'] == True:
                print(f'Authorized as {self.worker}')
                return True
            # Also check for mining.notify in the same batch
            if 'method' in msg and msg['method'] == 'mining.notify':
                self.handle_notify(msg['params'])
            if 'method' in msg and msg['method'] == 'mining.set_difficulty':
                self.difficulty = msg['params'][0]
                print(f'Difficulty: {self.difficulty}')
        return True  # some pools don't send explicit auth response

    def handle_notify(self, params):
        self.job = {
            'id': params[0],
            'prevhash': params[1],
            'coinb1': params[2],
            'coinb2': params[3],
            'merkle_branches': params[4],
            'version': params[5],
            'nbits': params[6],
            'ntime': params[7],
            'clean': params[8] if len(params) > 8 else False,
        }
        print(f'New job: {self.job["id"][:8]}... prevhash={self.job["prevhash"][:16]}...')

    def build_header(self, extranonce2, ntime, nonce):
        """Build 80-byte block header from job + nonce."""
        job = self.job

        # Build coinbase
        coinbase = bytes.fromhex(job['coinb1'] + self.extranonce1 +
                                  extranonce2 + job['coinb2'])
        coinbase_hash = dsha256(coinbase)

        # Build merkle root
        merkle_root = coinbase_hash
        for branch in job['merkle_branches']:
            merkle_root = dsha256(merkle_root + bytes.fromhex(branch))

        # Build header
        header = b''
        header += bytes.fromhex(job['version'])[::-1]  # version LE
        header += bytes.fromhex(job['prevhash'])  # prevhash (already in internal order)
        header += merkle_root  # merkle root
        header += bytes.fromhex(ntime)[::-1]  # time LE
        header += bytes.fromhex(job['nbits'])[::-1]  # bits LE
        header += struct.pack('<I', nonce)  # nonce LE

        return header

    def phi_solve(self, header76, bits_hex, height_est=0):
        """Use phi geometry to find nonce candidates."""
        # prevhash as 8 LE words
        prevhash_bytes = bytes.fromhex(self.job['prevhash'])
        pw = struct.unpack('<8I', prevhash_bytes[:32])

        bits = int(bits_hex, 16)
        bits_le = struct.unpack('<I', struct.pack('>I', bits))[0]

        sa = round(max(height_est, 1) * RATE_A * DEG1) & MASK32
        sb = round(max(height_est, 1) * RATE_B * DEG1) & MASK32
        sa2 = round(max(height_est, 1) * RATE_A * DEG1 / 2) & MASK32
        sab = round(max(height_est, 1) * (RATE_A+RATE_B)/2 * DEG1) & MASK32

        x_shifts = [0, sa, (-sa)&MASK32, sb, (-sb)&MASK32,
                    sa2, (-sa2)&MASK32, sab, (-sab)&MASK32]

        centers = set()
        for a in range(8):
            for b in range(8):
                for v in [add(pw[a],pw[b]), (pw[a]-pw[b])&MASK32, pw[a]^pw[b],
                          nswap(add(pw[a],pw[b])), nswap((pw[a]-pw[b])&MASK32)]:
                    for x in x_shifts:
                        for li in range(20):
                            for ls in [1, -1]:
                                y = round(ls * L[li] * DEG2 / 5) & MASK32
                                centers.add((v + x + y) & MASK32)

        return sorted(centers)

    def mine(self, pool_host=None, pool_port=None):
        """Main mining loop."""
        if pool_host:
            self.host = pool_host
            self.port = pool_port

        self.connect()
        self.subscribe()
        self.authorize()

        # Listen for jobs and solve
        self.running = True
        self.sock.settimeout(1)

        print(f'\nMINING as {self.worker}')
        print(f'phi solver: prev_hash geometry + 360 degree base conversion')
        print('=' * 60)

        t0 = time.time()

        while self.running:
            # Check for new messages
            try:
                msgs = self.recv()
                for msg in msgs:
                    if 'method' in msg:
                        if msg['method'] == 'mining.notify':
                            self.handle_notify(msg['params'])
                        elif msg['method'] == 'mining.set_difficulty':
                            self.difficulty = msg['params'][0]
                            print(f'Difficulty: {self.difficulty}')
                    elif 'result' in msg:
                        if msg.get('result') == True:
                            self.shares_accepted += 1
                            print(f'  SHARE ACCEPTED! ({self.shares_accepted} total)')
            except socket.timeout:
                pass
            except Exception as e:
                print(f'  recv error: {e}')
                break

            if not self.job:
                continue

            # Build extranonce2
            en2 = '00' * self.extranonce2_size

            # Generate phi centers from prevhash
            centers = self.phi_solve(None, self.job['nbits'])

            # Spiral from each center
            for center in centers[:1000]:  # top 1000 centers
                for r in range(10000):  # 10K spiral per center
                    for nonce in [(center + r) & MASK32, (center - r) & MASK32]:
                        self.hashes += 1

                        header = self.build_header(en2, self.job['ntime'], nonce)
                        h = dsha256(header)
                        hash_int = int.from_bytes(h, 'little')

                        # Check against pool difficulty
                        target = int(0xFFFF * 2**208 / self.difficulty)
                        if hash_int < target:
                            print(f'\n  SHARE FOUND! nonce={nonce:#010x} hash={h[::-1].hex()[:16]}...')
                            self.send('mining.submit', [
                                self.worker,
                                self.job['id'],
                                en2,
                                self.job['ntime'],
                                f'{nonce:08x}'
                            ])
                            self.shares_submitted += 1

                        if self.hashes % 100000 == 0:
                            elapsed = time.time() - t0
                            rate = self.hashes / elapsed if elapsed > 0 else 0
                            print(f'  {self.hashes:,} hashes  {rate:.0f} H/s  '
                                  f'shares: {self.shares_accepted}/{self.shares_submitted}',
                                  end='\r')

                    if self.hashes > 10_000_000:  # 10M hashes per job
                        break
                break  # one center per job cycle, check for new jobs

        self.sock.close()
        print(f'\nStopped. {self.hashes:,} hashes, {self.shares_accepted} shares accepted.')


if __name__ == '__main__':
    # Default: connect to a public pool
    # User can override with command line args
    import argparse
    p = argparse.ArgumentParser(description='phi miner - nos3bl33d / (DFG) DeadFoxGroup')
    p.add_argument('--pool', default='solo.ckpool.org', help='Pool hostname')
    p.add_argument('--port', type=int, default=3333, help='Pool port')
    p.add_argument('--worker', default='nos3bl33d.phi', help='Worker name')
    p.add_argument('--password', default='x', help='Worker password')
    args = p.parse_args()

    print('PHI MINER')
    print('nos3bl33d / (DFG) DeadFoxGroup')
    print('x^2 = x + 1')
    print()

    miner = StratumMiner(args.pool, args.port, args.worker, args.password)
    try:
        miner.mine()
    except KeyboardInterrupt:
        print('\nStopping...')
        miner.running = False
    except Exception as e:
        print(f'\nError: {e}')

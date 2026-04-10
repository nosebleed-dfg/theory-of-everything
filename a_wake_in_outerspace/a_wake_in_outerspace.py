"""
A WAKE IN OUTERSPACE
(DFG) DeadFoxGroup

Real Bitcoin miner. Solves blocks sequentially. Propagates normally.
Each coinbase embeds your tag — visible to every blockchain analyzer.
Like Satoshi's Times headline in the genesis block, but in every block.

Connects to Bitcoin node via JSON-RPC (getblocktemplate / submitblock).
Default: regtest mode for local testing. Switch to testnet/mainnet via config.

Flow:
  1. getblocktemplate — get current work from node
  2. Build coinbase tx with your tag in scriptSig
  3. Build merkle tree
  4. Construct block header
  5. Solve (find nonce via brute force)
  6. submitblock — propagate to network
  7. Wait for confirmation
  8. Next block. Repeat. Forever.

The big bang in digital space.
"""

import struct
import hashlib
import json
import time
import sys
import os
import binascii
from datetime import datetime, timezone

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

# ============================================================
# CONFIG
# ============================================================

DEFAULT_CONFIG = {
    'rpc_url': 'http://127.0.0.1:18443',  # regtest default
    'rpc_user': 'bitcoin',
    'rpc_pass': 'bitcoin',
    'network': 'regtest',                   # regtest / testnet / mainnet

    # Coinbase signature — embedded in every block we mine
    'coinbase_message': 'YOUR_TAG / YOUR_GROUP',

    # Wallet — P2PKH address to receive block reward
    # This is a test address for regtest. Replace for real networks.
    'reward_address': None,  # will be generated or fetched from node

    # Mining
    'poll_interval': 1.0,     # seconds between getblocktemplate calls
    'max_nonce': 0xFFFFFFFF,
    'log_interval': 500000,
}

# ============================================================
# BITCOIN PRIMITIVES
# ============================================================

def double_sha256(data):
    return hashlib.sha256(hashlib.sha256(data).digest()).digest()

def hash256(data):
    return double_sha256(data)

def sha256(data):
    return hashlib.sha256(data).digest()

def ripemd160(data):
    h = hashlib.new('ripemd160')
    h.update(data)
    return h.digest()

def hash160(data):
    return ripemd160(sha256(data))

def merkle_root(hashes):
    """Compute merkle root from list of tx hashes (bytes, internal order)."""
    if len(hashes) == 0:
        return b'\x00' * 32
    if len(hashes) == 1:
        return hashes[0]
    while len(hashes) > 1:
        if len(hashes) % 2 == 1:
            hashes.append(hashes[-1])  # duplicate last if odd
        new_level = []
        for i in range(0, len(hashes), 2):
            new_level.append(hash256(hashes[i] + hashes[i+1]))
        hashes = new_level
    return hashes[0]

def compact_size(n):
    """Bitcoin variable-length integer encoding."""
    if n < 0xfd:
        return struct.pack('<B', n)
    elif n <= 0xffff:
        return b'\xfd' + struct.pack('<H', n)
    elif n <= 0xffffffff:
        return b'\xfe' + struct.pack('<I', n)
    else:
        return b'\xff' + struct.pack('<Q', n)

def int_to_le_hex(n, length):
    """Integer to little-endian hex string."""
    return n.to_bytes(length, 'little').hex()

def hex_to_bytes(h):
    return bytes.fromhex(h)

def bytes_to_hex(b):
    return b.hex()

def reverse_hex(h):
    """Reverse byte order of hex string."""
    return bytes.fromhex(h)[::-1].hex()

# ============================================================
# COINBASE TRANSACTION
# ============================================================

def build_coinbase_tx(height, reward_sats, reward_script, message, extra_nonce=0, witness_commitment=None):
    """
    Build a coinbase transaction.

    The coinbase scriptSig contains:
      - Block height (BIP34)
      - Custom coinbase message
      - Extra nonce for additional search space

    The output pays the block reward to reward_script (P2PKH or P2SH).
    """
    # Version 2 (segwit-era)
    tx = struct.pack('<I', 2)

    # Segwit marker + flag
    tx += b'\x00\x01'

    # Input count
    tx += compact_size(1)

    # Coinbase input:
    #   prev_hash = 0x00...00 (32 bytes)
    #   prev_index = 0xFFFFFFFF
    #   scriptSig = height + message + extra_nonce
    #   sequence = 0xFFFFFFFF
    tx += b'\x00' * 32                        # prev hash (null for coinbase)
    tx += struct.pack('<I', 0xFFFFFFFF)        # prev index

    # ScriptSig: BIP34 height encoding
    # Bitcoin Core uses OP_N for heights 0-16, push-data for larger
    if height == 0:
        height_push = b'\x00'  # OP_0
    elif 1 <= height <= 16:
        height_push = bytes([0x50 + height])  # OP_1=0x51 through OP_16=0x60
    elif height <= 0x7f:
        height_push = b'\x01' + bytes([height])
    elif height <= 0x7fff:
        height_push = b'\x02' + struct.pack('<H', height)
    elif height <= 0x7fffff:
        height_push = b'\x03' + height.to_bytes(3, 'little')
    else:
        height_push = b'\x04' + struct.pack('<I', height)

    msg_bytes = message.encode('utf-8')
    msg_push = bytes([len(msg_bytes)]) + msg_bytes if len(msg_bytes) < 76 else b''

    extra_bytes = struct.pack('<Q', extra_nonce)
    extra_push = bytes([len(extra_bytes)]) + extra_bytes

    script_sig = height_push + msg_push + extra_push
    tx += compact_size(len(script_sig))
    tx += script_sig
    tx += struct.pack('<I', 0xFFFFFFFF)        # sequence

    # Outputs
    n_outputs = 1
    if witness_commitment:
        n_outputs = 2
    tx += compact_size(n_outputs)

    # Output 0: reward to address
    tx += struct.pack('<Q', reward_sats)       # value in satoshis
    tx += compact_size(len(reward_script))
    tx += reward_script

    # Output 1 (optional): witness commitment (segwit)
    if witness_commitment:
        wc_script = bytes.fromhex(witness_commitment)
        tx += struct.pack('<Q', 0)             # 0 value
        tx += compact_size(len(wc_script))
        tx += wc_script

    # Witness (coinbase must have a witness stack with one 32-byte zero element)
    tx += b'\x01'                              # 1 witness item
    tx += b'\x20'                              # 32 bytes
    tx += b'\x00' * 32                         # all zeros (coinbase witness nonce)

    # Locktime
    tx += struct.pack('<I', 0)

    return tx

def p2pkh_script(address_hash160):
    """Build P2PKH output script: OP_DUP OP_HASH160 <hash> OP_EQUALVERIFY OP_CHECKSIG"""
    return b'\x76\xa9\x14' + address_hash160 + b'\x88\xac'

# ============================================================
# BLOCK CONSTRUCTION
# ============================================================

def build_block_header(version, prev_hash, merkle, timestamp, bits, nonce):
    """Construct 80-byte block header."""
    return struct.pack('<I', version) + \
           prev_hash + \
           merkle + \
           struct.pack('<I', timestamp) + \
           struct.pack('<I', bits) + \
           struct.pack('<I', nonce)

def bits_to_target(bits):
    exp = bits >> 24
    coeff = bits & 0x7fffff
    if exp <= 3:
        return coeff >> (8 * (3 - exp))
    return coeff << (8 * (exp - 3))

def serialize_block(header, transactions):
    """Full block serialization: header + tx_count + transactions."""
    block = header
    block += compact_size(len(transactions))
    for tx in transactions:
        block += tx
    return block

# ============================================================
# RPC CLIENT
# ============================================================

class BitcoinRPC:
    def __init__(self, url, user, password):
        self.url = url
        self.user = user
        self.password = password
        self.id_counter = 0

    def call(self, method, params=None):
        if not HAS_REQUESTS:
            raise RuntimeError("requests library not installed. Run: pip install requests")
        self.id_counter += 1
        payload = {
            'jsonrpc': '2.0',
            'id': self.id_counter,
            'method': method,
            'params': params or [],
        }
        resp = requests.post(
            self.url,
            json=payload,
            auth=(self.user, self.password),
            timeout=30,
        )
        result = resp.json()
        if result.get('error'):
            raise RuntimeError(f"RPC error: {result['error']}")
        return result.get('result')

    def getblocktemplate(self):
        return self.call('getblocktemplate', [{'rules': ['segwit']}])

    def submitblock(self, block_hex):
        return self.call('submitblock', [block_hex])

    def getblockcount(self):
        return self.call('getblockcount')

    def getbestblockhash(self):
        return self.call('getbestblockhash')

    def getnewaddress(self):
        return self.call('getnewaddress')

    def getaddressinfo(self, addr):
        return self.call('getaddressinfo', [addr])

    def generatetoaddress(self, n, addr):
        """Regtest only: generate n blocks to address."""
        return self.call('generatetoaddress', [n, addr])

# ============================================================
# THE MINER
# ============================================================

class AxiomBitcoinMiner:
    def __init__(self, config=None):
        self.config = {**DEFAULT_CONFIG, **(config or {})}
        self.rpc = BitcoinRPC(
            self.config['rpc_url'],
            self.config['rpc_user'],
            self.config['rpc_pass'],
        )
        self.blocks_mined = 0
        self.total_reward = 0
        self.start_time = time.time()

    @staticmethod
    def _strip_witness(tx_bytes):
        """Remove segwit witness data for txid computation (BIP141)."""
        version = tx_bytes[:4]
        pos = 4
        has_witness = tx_bytes[pos] == 0x00 and tx_bytes[pos+1] == 0x01
        if has_witness:
            pos += 2
        # Scan through inputs
        inp_start = pos
        n_in = tx_bytes[pos]; pos += 1
        for _ in range(n_in):
            pos += 36
            ss_len = tx_bytes[pos]; pos += 1
            pos += ss_len + 4
        # Scan through outputs
        n_out = tx_bytes[pos]; pos += 1
        for _ in range(n_out):
            pos += 8
            sc_len = tx_bytes[pos]; pos += 1
            pos += sc_len
        out_end = pos
        locktime = tx_bytes[-4:]
        # Non-witness = version + inputs + outputs + locktime (no marker/flag/witness)
        return version + tx_bytes[inp_start:out_end] + locktime

    def get_reward_script(self):
        """Get or create a P2PKH script for block rewards."""
        if self.config['reward_address']:
            addr = self.config['reward_address']
        else:
            try:
                addr = self.rpc.getnewaddress()
                self.config['reward_address'] = addr
                print(f"  Wallet: {addr}")
            except Exception as e:
                # Fallback: use a test hash160
                print(f"  No wallet from node ({e}), using test address")
                test_hash160 = hash160(b"awake-dfg-test-wallet")
                return p2pkh_script(test_hash160)

        try:
            info = self.rpc.getaddressinfo(addr)
            script_hex = info.get('scriptPubKey', '')
            if script_hex:
                return bytes.fromhex(script_hex)
        except Exception:
            pass

        # Fallback P2PKH from address hash
        test_hash160 = hash160(b"awake-dfg-test-wallet")
        return p2pkh_script(test_hash160)

    def solve_block(self, header_template, target, nonce_start=0):
        """
        Find nonce such that hash(header) < target.
        Returns (nonce, hash, tries, elapsed) or (None, None, tries, elapsed).
        """
        t0 = time.time()
        max_nonce = self.config['max_nonce']
        log_interval = self.config['log_interval']

        for nonce in range(nonce_start, max_nonce + 1):
            header = header_template[:76] + struct.pack('<I', nonce)
            h = double_sha256(header)
            h_int = int.from_bytes(h[::-1], 'big')

            if h_int < target:
                return nonce, h, nonce - nonce_start + 1, time.time() - t0

            if nonce > nonce_start and (nonce - nonce_start) % log_interval == 0:
                elapsed = time.time() - t0
                rate = (nonce - nonce_start) / elapsed
                sys.stdout.write(f"\r    {nonce - nonce_start:,} hashes ({rate:,.0f} H/s)")
                sys.stdout.flush()

        return None, None, max_nonce - nonce_start + 1, time.time() - t0

    def mine_from_template(self, template):
        """Mine a block from getblocktemplate response."""
        version = template['version']
        prev_hash = bytes.fromhex(template['previousblockhash'])[::-1]  # internal byte order
        height = template['height']
        bits = int(template['bits'], 16)
        target = bits_to_target(bits)
        cur_time = template.get('curtime', int(time.time()))
        reward_sats = template.get('coinbasevalue', 5000000000)  # default 50 BTC regtest

        # Build coinbase with tag signature
        msg = self.config['coinbase_message']
        reward_script = self.get_reward_script()
        witness_commitment = template.get('default_witness_commitment')
        coinbase_tx = build_coinbase_tx(height, reward_sats, reward_script, msg,
                                         witness_commitment=witness_commitment)
        # txid = hash of tx WITHOUT witness (BIP141)
        # Strip: version(4) + marker(1) + flag(1) ... witness ... locktime(4)
        # Non-witness serialization: version + inputs + outputs + locktime
        coinbase_no_witness = self._strip_witness(coinbase_tx)
        coinbase_hash = hash256(coinbase_no_witness)

        # Gather transaction hashes (coinbase + mempool txs)
        tx_hashes = [coinbase_hash]
        txs = [coinbase_tx]

        for tx_data in template.get('transactions', []):
            tx_bytes = bytes.fromhex(tx_data['data'])
            txs.append(tx_bytes)
            # txid is hash of tx, but template gives it directly
            if 'txid' in tx_data:
                tx_hashes.append(bytes.fromhex(tx_data['txid'])[::-1])
            else:
                tx_hashes.append(hash256(tx_bytes))

        # Merkle root
        mk_root = merkle_root(tx_hashes)

        # Build header (nonce=0 placeholder)
        header = build_block_header(version, prev_hash, mk_root, cur_time, bits, 0)

        print(f"\n  Block #{height}")
        print(f"  Previous: {template['previousblockhash'][:16]}...")
        print(f"  Transactions: {len(txs)} (coinbase + {len(txs)-1} mempool)")
        print(f"  Reward: {reward_sats / 1e8:.8f} BTC")
        print(f"  Target: {hex(target)[:20]}...")
        print(f"  Coinbase: \"{msg}\"")

        # SOLVE
        print(f"  Mining...", end='', flush=True)
        nonce, block_hash, tries, elapsed = self.solve_block(header, target)

        if nonce is None:
            print(f"\r  Nonce space exhausted. No solution in {tries:,} tries.")
            return None

        print(f"\r  SOLVED! nonce={nonce} ({tries:,} hashes, {elapsed:.3f}s, {tries/elapsed:,.0f} H/s)")
        print(f"  Hash: {block_hash[::-1].hex()}")

        # Build final block
        final_header = build_block_header(version, prev_hash, mk_root, cur_time, bits, nonce)
        full_block = serialize_block(final_header, txs)

        return {
            'header': final_header,
            'block': full_block,
            'block_hex': full_block.hex(),
            'hash': block_hash[::-1].hex(),
            'height': height,
            'nonce': nonce,
            'tries': tries,
            'elapsed': elapsed,
            'reward_btc': reward_sats / 1e8,
            'coinbase_msg': msg,
            'tx_count': len(txs),
        }

    def submit_block(self, block_result):
        """Submit solved block to network."""
        print(f"  Submitting block #{block_result['height']}...")
        try:
            result = self.rpc.submitblock(block_result['block_hex'])
            if result is None:
                print(f"  ACCEPTED! Block #{block_result['height']} on chain.")
                return True
            else:
                print(f"  Rejected: {result}")
                return False
        except Exception as e:
            print(f"  Submit error: {e}")
            return False

    def wait_confirmation(self, expected_height, timeout=30):
        """Wait for block to appear on chain."""
        t0 = time.time()
        while time.time() - t0 < timeout:
            try:
                height = self.rpc.getblockcount()
                if height >= expected_height:
                    print(f"  Confirmed at height {height}.")
                    return True
            except Exception:
                pass
            time.sleep(0.5)
        print(f"  Confirmation timeout after {timeout}s")
        return False

    def run(self, n_blocks=None):
        """Main mining loop. Mine n_blocks or forever if None."""
        banner()
        print(f"  Network: {self.config['network']}")
        print(f"  RPC: {self.config['rpc_url']}")
        print(f"  Signature: {self.config['coinbase_message']}")

        # Test RPC connection
        try:
            height = self.rpc.getblockcount()
            best = self.rpc.getbestblockhash()
            print(f"  Connected! Height: {height}, tip: {best[:16]}...")
        except Exception as e:
            print(f"\n  Cannot connect to Bitcoin node: {e}")
            print(f"  For local testing, start Bitcoin Core in regtest mode:")
            print(f"    bitcoind -regtest -rpcuser=bitcoin -rpcpassword=bitcoin -daemon")
            print(f"\n  Or run in offline demo mode:")
            self.run_offline_demo(n_blocks or 10)
            return

        count = 0
        while n_blocks is None or count < n_blocks:
            try:
                template = self.rpc.getblocktemplate()
            except Exception as e:
                print(f"\n  getblocktemplate failed: {e}")
                print(f"  Retrying in {self.config['poll_interval']}s...")
                time.sleep(self.config['poll_interval'])
                continue

            result = self.mine_from_template(template)
            if result is None:
                continue

            accepted = self.submit_block(result)
            if accepted:
                self.blocks_mined += 1
                self.total_reward += result['reward_btc']
                count += 1

                self.wait_confirmation(result['height'])

                print(f"\n  === STATS ===")
                print(f"  Blocks mined: {self.blocks_mined}")
                print(f"  Total reward: {self.total_reward:.8f} BTC")
                print(f"  Uptime: {time.time() - self.start_time:.1f}s")
                print(f"  Tag embedded in every block")
                print()

            time.sleep(0.1)  # brief pause before next block

        print(f"\n  Mining complete. {self.blocks_mined} blocks mined.")
        print(f"  Total reward: {self.total_reward:.8f} BTC")
        print(f"  Signature in every coinbase: {self.config['coinbase_message']}")

    def run_offline_demo(self, n_blocks=10):
        """Offline demo: simulate mining without a node."""
        print(f"\n  OFFLINE DEMO MODE")
        print(f"  Simulating {n_blocks} blocks (no network)")
        print(f"  Difficulty: regtest (easy)")
        print()

        prev_hash = b'\x00' * 32
        height = 0
        reward_sats = 5000000000  # 50 BTC regtest
        bits = 0x207fffff  # regtest minimum difficulty
        target = bits_to_target(bits)
        reward_script = p2pkh_script(hash160(b"awake-dfg-test"))

        snapshots = []

        for i in range(n_blocks):
            msg = self.config['coinbase_message']
            coinbase = build_coinbase_tx(height, reward_sats, reward_script, msg, extra_nonce=i)
            coinbase_hash = hash256(coinbase)
            mk_root = merkle_root([coinbase_hash])

            timestamp = int(time.time())
            header = build_block_header(1, prev_hash, mk_root, timestamp, bits, 0)

            sys.stdout.write(f"  Block #{height}: mining...")
            sys.stdout.flush()

            nonce, block_hash, tries, elapsed = self.solve_block(header, target)

            if nonce is not None:
                final_header = build_block_header(1, prev_hash, mk_root, timestamp, bits, nonce)
                full_block = serialize_block(final_header, [coinbase])

                sys.stdout.write(f"\r  Block #{height}: {block_hash[::-1].hex()[:16]}...")
                print(f" nonce={nonce} ({tries:,} H, {elapsed:.3f}s)")
                print(f"    Coinbase: \"{msg}\"")
                print(f"    Reward: {reward_sats/1e8:.2f} BTC")

                snapshots.append({
                    'height': height,
                    'hash': block_hash[::-1].hex(),
                    'nonce': nonce,
                    'tries': tries,
                    'time': elapsed,
                    'coinbase': msg,
                    'reward_btc': reward_sats / 1e8,
                    'timestamp': datetime.fromtimestamp(timestamp, tz=timezone.utc).isoformat(),
                    'block_size': len(full_block),
                })

                prev_hash = block_hash
                height += 1
                self.blocks_mined += 1
                self.total_reward += reward_sats / 1e8
            else:
                print(f"\r  Block #{height}: FAILED")

        # Save snapshot
        snap_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'snapshots')
        os.makedirs(snap_dir, exist_ok=True)
        snap_file = os.path.join(snap_dir, f'wake_{int(time.time())}.json')
        snapshot = {
            'event': 'A Wake In Outerspace',
            'author': '(DFG) DeadFoxGroup',
            'timestamp': datetime.now(timezone.utc).isoformat(),
            'blocks_mined': self.blocks_mined,
            'total_reward_btc': self.total_reward,
            'signature_in_every_block': True,
            'signature': self.config['coinbase_message'],
            'blocks': snapshots,
        }
        with open(snap_file, 'w') as f:
            json.dump(snapshot, f, indent=2)

        print(f"\n  ============================================")
        print(f"  THE WAKE IS COMPLETE")
        print(f"  ============================================")
        print(f"  Blocks: {self.blocks_mined}")
        print(f"  Reward: {self.total_reward:.2f} BTC")
        print(f"  Signature: \"{self.config['coinbase_message']}\"")
        print(f"  Present in: EVERY block")
        print(f"  Snapshot: {os.path.basename(snap_file)}")
        print()
        print(f"  Every block stamped. Every nonce solved.")
        print(f"  Signature in every coinbase: {self.config['coinbase_message']}")
        print()
        print(f"  All is number.")
        print()

# ============================================================
# BANNER
# ============================================================

def banner():
    print()
    print("  ============================================")
    print("  A   W A K E   I N   O U T E R S P A C E")
    print("  ============================================")
    print("  (DFG) DeadFoxGroup")
    print()
    print("  Your tag in every coinbase.")
    print("  The big bang in digital space.")
    print()

# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    n_blocks = 10
    config = {}

    # Parse simple args
    args = sys.argv[1:]
    for i, arg in enumerate(args):
        if arg.isdigit():
            n_blocks = int(arg)
        elif arg == '--testnet':
            config['network'] = 'testnet'
            config['rpc_url'] = 'http://127.0.0.1:18332'
        elif arg == '--mainnet':
            config['network'] = 'mainnet'
            config['rpc_url'] = 'http://127.0.0.1:8332'
        elif arg == '--regtest':
            config['network'] = 'regtest'
            config['rpc_url'] = 'http://127.0.0.1:18443'
        elif arg.startswith('--rpc='):
            config['rpc_url'] = arg.split('=', 1)[1]
        elif arg.startswith('--user='):
            config['rpc_user'] = arg.split('=', 1)[1]
        elif arg.startswith('--pass='):
            config['rpc_pass'] = arg.split('=', 1)[1]
        elif arg.startswith('--msg='):
            config['coinbase_message'] = arg.split('=', 1)[1]

    miner = AxiomBitcoinMiner(config)
    miner.run(n_blocks=n_blocks)

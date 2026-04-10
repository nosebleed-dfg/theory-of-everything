"""
AXIOM MINER — Solo Bitcoin Mining via Stratum
(DFG) DeadFoxGroup | nos3bl33d | April 2026

x^2 = x + 1 | The nonce is forced, not found.

Standalone CLI miner. Connects to solo.ckpool.org via Stratum v1.
Dual-layer solving: axiom-guided seed search + sequential sweep.
Midstate optimization via hashlib.copy() for maximum throughput.

Usage: python miner.py
"""

import struct
import hashlib
import socket
import json
import time
import threading
import sys

# ============================================================
# CONSTANTS
# ============================================================

PHI = (1 + 5 ** 0.5) / 2
MASK32 = 0xFFFFFFFF

WALLET = "bc1qnuc5nkwjls0lc3zmek9k6tx9r8n77p03337qjv"
WORKER = "axiom"
POOL_HOST = "solo.ckpool.org"
POOL_PORT = 3333
USER_AGENT = "AxiomMiner/1.0"

# Pool difficulty 1 target
DIFF1_TARGET = 0x00000000ffff0000000000000000000000000000000000000000000000000000


# ============================================================
# STRATUM CLIENT
# ============================================================

class StratumClient:
    """Stratum v1 client for Bitcoin mining pools."""

    def __init__(self, host, port, wallet, worker):
        self.host = host
        self.port = port
        self.wallet = wallet
        self.worker = worker

        self.sock = None
        self.buf = ""
        self.msg_id = 0
        self.lock = threading.Lock()

        self.extranonce1 = ""
        self.extranonce2_size = 4

        self.difficulty = 1
        self.target = DIFF1_TARGET
        self.current_job = None
        self.job_lock = threading.Lock()
        self.job_counter = 0  # increments on every new job

        self.accepted = 0
        self.rejected = 0
        self.connected = False

    def connect(self):
        """Connect and perform subscribe + authorize handshake."""
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.settimeout(30)
        self.sock.connect((self.host, self.port))
        self.buf = ""
        self.connected = True

        # mining.subscribe
        result = self._call("mining.subscribe", [USER_AGENT])
        if not isinstance(result, list) or len(result) < 3:
            raise ConnectionError(f"subscribe failed: {result}")
        self.extranonce1 = result[1]
        self.extranonce2_size = result[2]

        # Request lower difficulty (pool may ignore this)
        self._send_msg("mining.suggest_difficulty", [256])

        # mining.authorize (d=256 in password for pools that support it)
        self._call("mining.authorize", [f"{self.wallet}.{self.worker}", "d=256"])
        return True

    def _send_msg(self, method, params):
        """Send JSON-RPC message. Returns message id."""
        with self.lock:
            self.msg_id += 1
            mid = self.msg_id
        line = json.dumps({"id": mid, "method": method, "params": params}) + "\n"
        try:
            self.sock.sendall(line.encode())
        except (OSError, BrokenPipeError) as exc:
            self.connected = False
            raise ConnectionError(f"Send failed: {exc}") from exc
        return mid

    def _call(self, method, params, timeout=15):
        """Send request and wait for matching response."""
        mid = self._send_msg(method, params)
        deadline = time.time() + timeout
        while time.time() < deadline:
            for msg in self._recv_messages(timeout=1.0):
                self._dispatch(msg)
                if isinstance(msg, dict) and msg.get("id") == mid:
                    err = msg.get("error")
                    if err is not None:
                        raise ConnectionError(f"RPC error: {err}")
                    return msg.get("result")
        return None

    def _recv_messages(self, timeout=0.05):
        """Receive available data and parse complete JSON lines."""
        messages = []
        self.sock.settimeout(timeout)
        try:
            while True:
                try:
                    data = self.sock.recv(4096)
                    if not data:
                        self.connected = False
                        break
                    self.buf += data.decode('utf-8', errors='replace')
                except socket.timeout:
                    break
                except OSError:
                    self.connected = False
                    break
        finally:
            self.sock.settimeout(30)

        while '\n' in self.buf:
            line, self.buf = self.buf.split('\n', 1)
            line = line.strip()
            if not line:
                continue
            try:
                messages.append(json.loads(line))
            except json.JSONDecodeError:
                continue
        return messages

    def _dispatch(self, msg):
        """Handle incoming stratum message."""
        if not isinstance(msg, dict):
            return

        method = msg.get("method")

        if method == "mining.notify":
            p = msg["params"]
            clean = p[8] if len(p) > 8 else False
            job = {
                "id": p[0], "prevhash": p[1], "coinb1": p[2], "coinb2": p[3],
                "merkle_branches": p[4], "version": p[5], "nbits": p[6],
                "ntime": p[7], "clean": clean,
            }
            with self.job_lock:
                old_id = self.current_job["id"] if self.current_job else None
                old_prev = self.current_job["prevhash"] if self.current_job else None
                self.current_job = job
                # Only signal new work if prevhash changed (new block) or clean_jobs
                if job["prevhash"] != old_prev or clean:
                    self.job_counter += 1

        elif method == "mining.set_difficulty":
            self.difficulty = msg["params"][0]
            self.target = int(DIFF1_TARGET / max(1, self.difficulty))

    def poll(self):
        """Poll for new messages."""
        try:
            for msg in self._recv_messages(timeout=0.05):
                self._dispatch(msg)
        except Exception:
            pass

    def get_job(self, timeout=30):
        """Wait for and return a job."""
        with self.job_lock:
            if self.current_job is not None:
                return self.current_job
        deadline = time.time() + timeout
        while time.time() < deadline:
            self.poll()
            with self.job_lock:
                if self.current_job is not None:
                    return self.current_job
            time.sleep(0.1)
        return None

    def submit_share(self, job_id, en2_hex, ntime_hex, nonce_hex):
        """
        Submit a share. Returns True if accepted, (False, msg) if rejected.
        Nonce hex is the hex of the 4 LE header bytes (stratum convention).
        """
        params = [f"{self.wallet}.{self.worker}", job_id, en2_hex, ntime_hex, nonce_hex]
        mid = self._send_msg("mining.submit", params)
        deadline = time.time() + 10
        while time.time() < deadline:
            for msg in self._recv_messages(timeout=1.0):
                self._dispatch(msg)
                if isinstance(msg, dict) and msg.get("id") == mid:
                    result = msg.get("result")
                    error = msg.get("error")
                    if result is True:
                        self.accepted += 1
                        return True
                    else:
                        self.rejected += 1
                        err_msg = ""
                        if isinstance(error, list) and len(error) >= 2:
                            err_msg = str(error[1])
                        elif error:
                            err_msg = str(error)
                        return (False, err_msg)
        return None

    def close(self):
        """Close the connection."""
        self.connected = False
        try:
            self.sock.shutdown(socket.SHUT_RDWR)
        except Exception:
            pass
        try:
            self.sock.close()
        except Exception:
            pass


# ============================================================
# BLOCK HEADER CONSTRUCTION
# ============================================================

def build_header(job, extranonce1, extranonce2_hex):
    """
    Build an 80-byte Bitcoin block header from stratum job data.

    Stratum prevhash: 8 groups of 4-byte words, each group needs byte-reversal
    to produce the correct internal byte order.

    Returns: bytearray(80) with nonce = 0
    """
    # Coinbase = coinb1 + extranonce1 + extranonce2 + coinb2
    coinbase = bytes.fromhex(job["coinb1"] + extranonce1 + extranonce2_hex + job["coinb2"])
    cb_hash = hashlib.sha256(hashlib.sha256(coinbase).digest()).digest()

    # Merkle root
    merkle = cb_hash
    for branch in job["merkle_branches"]:
        merkle = hashlib.sha256(hashlib.sha256(merkle + bytes.fromhex(branch)).digest()).digest()

    # Version (LE)
    version = struct.pack('<I', int(job["version"], 16))

    # Prevhash: swap each 4-byte group
    ph = job["prevhash"]
    prevhash = b""
    for i in range(0, 64, 8):
        prevhash += bytes.fromhex(ph[i:i+8])[::-1]

    # ntime + nbits (LE)
    ntime = struct.pack('<I', int(job["ntime"], 16))
    nbits = struct.pack('<I', int(job["nbits"], 16))

    header = bytearray(80)
    header[0:4] = version
    header[4:36] = prevhash
    header[36:68] = merkle
    header[68:72] = ntime
    header[72:76] = nbits
    # header[76:80] = nonce (zeros)
    return header


# ============================================================
# MINING ENGINE
# ============================================================

class MiningEngine:
    """
    High-performance mining engine with midstate optimization.

    Two-phase solving:
    1. Axiom phase: Planck seed + golden spiral candidates
    2. Sequential sweep: brute force nonce 0 -> 2^32
    """

    def __init__(self, stratum):
        self.stratum = stratum
        self.running = False
        self.total_hashes = 0
        self.hash_lock = threading.Lock()

    def start(self):
        """Start mining (blocks until stopped or disconnected)."""
        self.running = True
        threading.Thread(target=self._receiver_loop, daemon=True).start()
        threading.Thread(target=self._report_loop, daemon=True).start()
        self._mine_loop()

    def stop(self):
        self.running = False

    def _receiver_loop(self):
        while self.running:
            try:
                self.stratum.poll()
            except Exception:
                pass
            time.sleep(0.05)

    def _report_loop(self):
        last_time = time.time()
        last_hashes = 0
        while self.running:
            time.sleep(5)
            now = time.time()
            with self.hash_lock:
                cur = self.total_hashes
            dt = now - last_time
            dh = cur - last_hashes
            if dt > 0:
                rate = dh / dt
                unit, disp = "H/s", rate
                if rate >= 1_000_000:
                    unit, disp = "MH/s", rate / 1_000_000
                elif rate >= 1_000:
                    unit, disp = "KH/s", rate / 1_000
                print(f"  [{time.strftime('%H:%M:%S')}] {disp:>8.2f} {unit} | "
                      f"total: {cur:>12,} | diff: {self.stratum.difficulty} | "
                      f"shares: {self.stratum.accepted}A/{self.stratum.rejected}R")
                sys.stdout.flush()
            last_time = now
            last_hashes = cur

    def _mine_loop(self):
        en2_counter = 0
        pack_into = struct.pack_into  # local reference for speed
        sha256 = hashlib.sha256      # local reference for speed

        while self.running:
            job = self.stratum.get_job(timeout=30)
            if job is None:
                print("  [!] No job received, retrying...")
                sys.stdout.flush()
                continue

            job_id = job["id"]
            difficulty = self.stratum.difficulty
            target = self.stratum.target
            ntime_hex = job["ntime"]

            # Record job counter to detect new jobs
            with self.stratum.job_lock:
                start_job_counter = self.stratum.job_counter

            # Target as bytes for comparison
            target_bytes = target.to_bytes(32, 'big')

            # Count leading zero bytes for fast rejection
            zero_bytes = 0
            for b in target_bytes:
                if b == 0:
                    zero_bytes += 1
                else:
                    break

            # Generate extranonce2
            en2_hex = format(en2_counter, '0' + str(self.stratum.extranonce2_size * 2) + 'x')
            en2_counter = (en2_counter + 1) & ((1 << (self.stratum.extranonce2_size * 8)) - 1)

            # Build header
            try:
                header = build_header(job, self.stratum.extranonce1, en2_hex)
            except Exception as exc:
                print(f"  [!] Header build error: {exc}")
                sys.stdout.flush()
                time.sleep(1)
                continue

            print(f"  [{time.strftime('%H:%M:%S')}] Job {job_id[:8]}... | "
                  f"diff: {difficulty} | en2: {en2_hex[:8]}...")
            sys.stdout.flush()

            # Midstate: hash first 64 bytes once
            midstate_sha = sha256(bytes(header[:64]))
            tail = bytearray(header[64:80])

            # ============================================================
            # PHASE 1: AXIOM CANDIDATES
            # ============================================================
            prevhash_int = int.from_bytes(header[4:36], 'little')
            planck_seed = ((prevhash_int * 3) // (50 * (1 << 224))) & MASK32
            diff_int = max(1, int(difficulty))

            # Build axiom candidate list
            candidates = []

            # Planck neighborhood
            radius = min(diff_int * 2, 200_000)
            for delta in range(radius):
                candidates.append((planck_seed + delta) & MASK32)
                if delta > 0:
                    candidates.append((planck_seed - delta) & MASK32)

            # Golden spiral
            for k in range(64):
                sp = int(planck_seed * (PHI ** k)) & MASK32
                candidates.append(sp)
                for d in range(1, min(diff_int, 200)):
                    candidates.append((sp + d) & MASK32)
                    candidates.append((sp - d) & MASK32)

            # Mine axiom candidates
            found = False
            axiom_hashed = 0

            for nonce in candidates:
                pack_into('<I', tail, 12, nonce)
                inner = midstate_sha.copy()
                inner.update(tail)
                result = sha256(inner.digest()).digest()

                axiom_hashed += 1

                # Fast reject: check high bytes of hash (Bitcoin LE comparison)
                if result[31] == 0 and result[30] == 0 and result[29] == 0 and result[28] == 0:
                    if result[27] == 0 and result[::-1] < target_bytes:
                        with self.hash_lock:
                            self.total_hashes += axiom_hashed
                        # Nonce for stratum: LE header bytes as hex
                        nonce_hex = struct.pack('<I', nonce).hex()
                        self._submit(job_id, en2_hex, ntime_hex, nonce_hex,
                                     result, nonce, "axiom")
                        found = True
                        break

                # Check for new job every 50K hashes
                if axiom_hashed % 50_000 == 0:
                    with self.hash_lock:
                        self.total_hashes += 50_000
                    axiom_hashed -= 50_000
                    with self.stratum.job_lock:
                        if self.stratum.job_counter != start_job_counter:
                            break
                    if not self.running:
                        return

            # Flush remaining axiom hash count
            with self.hash_lock:
                self.total_hashes += axiom_hashed

            if found:
                with self.stratum.job_lock:
                    if self.stratum.current_job and self.stratum.current_job["id"] == job_id:
                        self.stratum.current_job = None
                continue

            # Check if job was superseded during axiom phase
            with self.stratum.job_lock:
                if self.stratum.job_counter != start_job_counter:
                    continue

            # ============================================================
            # PHASE 2: SEQUENTIAL SWEEP
            # ============================================================
            BATCH = 100_000
            nonce = 0
            found = False

            while nonce <= MASK32 and self.running:
                batch_end = min(nonce + BATCH, MASK32 + 1)
                local_count = 0

                for n in range(nonce, batch_end):
                    pack_into('<I', tail, 12, n)
                    inner = midstate_sha.copy()
                    inner.update(tail)
                    result = sha256(inner.digest()).digest()

                    # Fast reject: for diff 10000, need ~5 leading zero bytes
                    if result[31] == 0 and result[30] == 0 and result[29] == 0 and result[28] == 0:
                        if result[27] == 0 and result[::-1] < target_bytes:
                            local_count = n - nonce + 1
                            with self.hash_lock:
                                self.total_hashes += local_count
                            nonce_hex = struct.pack('<I', n).hex()
                            self._submit(job_id, en2_hex, ntime_hex, nonce_hex,
                                         result, n, "sweep")
                            found = True
                            break

                if found:
                    break

                local_count = batch_end - nonce
                with self.hash_lock:
                    self.total_hashes += local_count
                nonce = batch_end

                # Check for new job between batches
                with self.stratum.job_lock:
                    if self.stratum.job_counter != start_job_counter:
                        break

            if found:
                with self.stratum.job_lock:
                    if self.stratum.current_job and self.stratum.current_job["id"] == job_id:
                        self.stratum.current_job = None

    def _submit(self, job_id, en2_hex, ntime_hex, nonce_hex, result_hash, nonce_int, source):
        """Submit a found share to the pool."""
        disp = result_hash[::-1].hex()
        print(f"\n  *** SHARE FOUND ({source}) ***")
        print(f"  Nonce:  {nonce_int} (0x{nonce_int:08x})")
        print(f"  Hash:   {disp}")
        print(f"  Submit: job={job_id[:8]} en2={en2_hex} ntime={ntime_hex} nonce={nonce_hex}")
        sys.stdout.flush()

        result = self.stratum.submit_share(job_id, en2_hex, ntime_hex, nonce_hex)

        if result is True:
            print(f"  >>> ACCEPTED! ({self.stratum.accepted} total) <<<\n")
        elif isinstance(result, tuple) and len(result) == 2:
            print(f"  >>> REJECTED: {result[1]} ({self.stratum.rejected} total) <<<\n")
        else:
            print(f"  >>> Submit response: {result} <<<\n")
        sys.stdout.flush()


# ============================================================
# MAIN
# ============================================================

def main():
    print()
    print("  ============================================")
    print("  AXIOM MINER v1.0")
    print("  (DFG) DeadFoxGroup | nos3bl33d")
    print("  x^2 = x + 1 | The nonce is forced.")
    print("  ============================================")
    print()
    print(f"  Pool:   {POOL_HOST}:{POOL_PORT}")
    print(f"  Wallet: {WALLET}")
    print(f"  Worker: {WORKER}")
    print()
    sys.stdout.flush()

    while True:
        stratum = StratumClient(POOL_HOST, POOL_PORT, WALLET, WORKER)

        try:
            print(f"  [{time.strftime('%H:%M:%S')}] Connecting to {POOL_HOST}:{POOL_PORT}...")
            sys.stdout.flush()
            stratum.connect()
            print(f"  [{time.strftime('%H:%M:%S')}] Connected!")
            print(f"  Extranonce1: {stratum.extranonce1}")
            print(f"  Extranonce2 size: {stratum.extranonce2_size}")
            print(f"  Difficulty: {stratum.difficulty}")
            print()
            sys.stdout.flush()

            # Let initial difficulty + job arrive
            time.sleep(1)
            stratum.poll()

            print(f"  [{time.strftime('%H:%M:%S')}] Difficulty: {stratum.difficulty}")
            tgt = stratum.target
            print(f"  [{time.strftime('%H:%M:%S')}] Target: {tgt.to_bytes(32,'big')[:8].hex()}...")
            print(f"  [{time.strftime('%H:%M:%S')}] Mining started!")
            print()
            sys.stdout.flush()

            engine = MiningEngine(stratum)
            engine.start()

        except KeyboardInterrupt:
            print("\n  [*] Shutting down...")
            sys.stdout.flush()
            stratum.close()
            break

        except Exception as exc:
            print(f"  [!] Error: {exc}")
            sys.stdout.flush()
            stratum.close()
            print(f"  [{time.strftime('%H:%M:%S')}] Reconnecting in 5s...")
            sys.stdout.flush()
            time.sleep(5)


if __name__ == "__main__":
    main()

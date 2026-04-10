"""
phi_chain_local.py — The entire Bitcoin blockchain from n0. One calc per block.
nos3bl33d

No mining. No searching. No internet.
Each block: nonce = phi(halvening), hash = dsha256(header). Next.
The whole chain in seconds.

x^2 = x + 1.
"""

import struct, hashlib, math, time, sys, os
sys.stdout.reconfigure(encoding='utf-8')

PHI    = (1 + math.sqrt(5)) / 2
MASK32 = 0xFFFFFFFF
MOD32  = 2**32

GENESIS_HEADER = bytes.fromhex(
    '0100000000000000000000000000000000000000000000000000000000000000'
    '000000003ba3edfd7a7b12b27ac72c3e67768f617fc81bc3888a51323a9fb8aa'
    '4b1e5e4a29ab5f49ffff001d1dac2b7c'
)
GENESIS_NONCE = 2083236893
GENESIS_BITS  = 0x1d00ffff
HALVENING     = 210_000

OUTPUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "btc_genesis_chain.txt")

def dsha256(data):
    return hashlib.sha256(hashlib.sha256(data).digest()).digest()

def quarter_step(k):
    v = MOD32 / PHI**(k / 4.0)
    return round(v) if v >= 1 else 0

def mirror(n):
    return struct.unpack('>I', struct.pack('<I', n & MASK32))[0]

def phi_nonce(halvening):
    """ONE CALC. The nonce for this halvening period."""
    k = 6 + halvening * 2
    return quarter_step(k)


def run(total_blocks):
    t0 = time.time()

    print(f"PHI CHAIN: {total_blocks:,} blocks from n0 = {GENESIS_NONCE}")

    # Genesis
    prev_hash = dsha256(GENESIS_HEADER)

    lines = []
    lines.append(f"0\t{prev_hash[::-1].hex()}\t{GENESIS_NONCE}")

    current_halvening = 0
    current_nonce = GENESIS_NONCE

    for height in range(1, total_blocks):
        halvening = height // HALVENING

        if halvening != current_halvening:
            current_halvening = halvening
            current_nonce = phi_nonce(halvening)
            if current_nonce < 1:
                print(f"  ZERO SPACE at block {height:,}")
                break

        # Build header: version + prev_hash + merkle(=hash of height+prev) + ts + bits + nonce
        merkle = dsha256(struct.pack('<I', height) + prev_hash)
        header = struct.pack('<I', 1)           # version
        header += prev_hash                      # prev hash (already LE from dsha256)
        header += merkle                         # merkle
        header += struct.pack('<I', height)      # timestamp = height (compact)
        header += struct.pack('<I', GENESIS_BITS)# bits (constant for phi chain)
        header += struct.pack('<I', current_nonce & MASK32)

        block_hash = dsha256(header)
        prev_hash = block_hash

        lines.append(f"{height}\t{block_hash[::-1].hex()}\t{current_nonce}")

        if height % 500_000 == 0:
            elapsed = time.time() - t0
            rate = height / elapsed
            print(f"  {height:>9,}  halv {halvening:>2}  nonce {current_nonce:>14,}  "
                  f"{rate:,.0f} blk/s  {elapsed:.1f}s")

    elapsed = time.time() - t0
    total = len(lines)
    rate = total / elapsed if elapsed > 0 else 0

    print(f"\n  {total:,} blocks in {elapsed:.1f}s ({rate:,.0f} blocks/sec)")

    # Write
    print(f"  Writing {OUTPUT}...")

    with open(OUTPUT, 'w', encoding='utf-8') as f:
        f.write("# BTC PHI CHAIN — from n0, no internet, one calc per block\n")
        f.write(f"# nos3bl33d | x^2 = x + 1 | {total:,} blocks\n")
        f.write(f"# Generated {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}\n")
        f.write(f"# n0 = {GENESIS_NONCE} | grid = 2^32/phi^(k/4)\n")
        f.write("# HEIGHT\tHASH\tNONCE\n")
        f.write('\n'.join(lines))
        f.write('\n')

    fsize = os.path.getsize(OUTPUT)
    print(f"  {fsize/1024/1024:.1f} MB written")
    print(f"  done.")


if __name__ == "__main__":
    # All 33 halvenings = ~6.93M blocks = the whole thing to zero space
    run(33 * HALVENING)

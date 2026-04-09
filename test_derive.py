"""
Test: does phi^291 actually DERIVE the Bitcoin nonce?

F291 = Fibonacci(291) mod 2^32 = phi^291 integer coefficient
F290 = Fibonacci(290) mod 2^32 = phi^291 constant term

derivation: nonce = F291 * sha256d(h76)[:4] + F290 mod 2^32
            (i.e., phi^291 * seed mod 2^32)

We test this against REAL difficulty targets, not trivially easy ones.
With BITS = 1e0fffff (~1M expected hashes), if the center is always
near-zero steps away, the math works. If steps ~ 500K, it's random.

nos3bl33d
"""

import struct, hashlib, time, urllib.request

F291 = 0x85824b02
F290 = 0xc3e3c441

def sha256d(d): return hashlib.sha256(hashlib.sha256(d).digest()).digest()

def fetch_tip():
    try:
        def g(p):
            r = urllib.request.Request("https://mempool.space/api/" + p,
                headers={"User-Agent": "AWake/1.0"})
            with urllib.request.urlopen(r, timeout=10) as x:
                return x.read().decode().strip()
        return g("blocks/tip/hash"), int(g("blocks/tip/height"))
    except Exception as e:
        print(f"  API fail: {e}")
        return "0"*64, 0

def make_cb(height, tag):
    hb = height.to_bytes((height.bit_length()+8)//8, "little")
    tb = tag.encode()
    ss  = bytes([len(hb)]) + hb + bytes([len(tb)]) + tb
    os_ = b"\x6a" + bytes([len(tb)]) + tb
    tx  = struct.pack("<I", 1) + b"\x01"
    tx += b"\x00"*32 + b"\xff\xff\xff\xff"
    tx += bytes([len(ss)]) + ss + b"\xff\xff\xff\xff"
    tx += b"\x01" + struct.pack("<q", 0) + bytes([len(os_)]) + os_
    tx += struct.pack("<I", 0)
    return tx

def bits_to_target(bits_hex):
    bi = int(bits_hex, 16)
    exp = bi >> 24
    coeff = bi & 0x7fffff
    return coeff << (8*(exp-3))

TAG = "nos3bl33d / DeadFoxGroup"

print("=" * 60)
print("  PHI^291 DERIVATION TEST")
print("  nos3bl33d")
print("=" * 60)
print()

# Test three difficulty levels
difficulties = [
    ("trivial  BITS=207fffff", "207fffff", "any nonce works — no info"),
    ("easy     BITS=1f3fffff", "1f3fffff", "~16K expected hashes"),
    ("medium   BITS=1e0fffff", "1e0fffff", "~1M expected hashes"),
]

tip_hash, tip_height = fetch_tip()
print(f"  tip #{tip_height}  {tip_hash[:20]}...")
print()

prev_bytes = bytes.fromhex(tip_hash)[::-1]

for label, bits_hex, note in difficulties:
    tgt = bits_to_target(bits_hex)
    print(f"  --- {label}  ({note}) ---")

    results = []
    h_prev = prev_bytes

    for i in range(5):
        height = tip_height + i + 1
        cb  = make_cb(height, TAG)
        mk  = sha256d(cb)
        hdr76 = (struct.pack("<I", 0x20000000) +
                 h_prev +
                 mk +
                 struct.pack("<I", int(time.time())) +
                 bytes.fromhex(bits_hex)[::-1])

        # PURE DERIVATION: no spiral, just phi^291 * seed
        seed = struct.unpack(">I", sha256d(hdr76)[:4])[0]
        nonce_derived = (F291 * seed + F290) & 0xFFFFFFFF
        h = sha256d(hdr76 + struct.pack("<I", nonce_derived))
        h_int = int.from_bytes(h[::-1], "big")
        derived_hit = h_int < tgt

        if derived_hit:
            results.append(("DERIVED", 0, nonce_derived))
            h_prev = h
        else:
            # how far is the center from a valid nonce?
            step = 0
            best = h_int
            for s in range(1, 100_001):
                off = (s+1)//2 if s & 1 else -(s//2)
                nc = (nonce_derived + off) & 0xFFFFFFFF
                hc = sha256d(hdr76 + struct.pack("<I", nc))
                hi = int.from_bytes(hc[::-1], "big")
                if hi < tgt:
                    step = s
                    h_prev = hc
                    results.append(("spiral", step, nc))
                    break
            else:
                results.append(("timeout", 100_001, 0))
                h_prev = bytes(32)

    print(f"  {'block':>6} {'result':>10} {'steps':>8}")
    for j, (res, steps, nc) in enumerate(results):
        print(f"  {tip_height+j+1:>6} {res:>10} {steps:>8}")

    direct_hits = sum(1 for r,_,_ in results if r == "DERIVED")
    avg_steps = sum(s for _,s,_ in results) / len(results) if results else 0
    expected = (2**256 // tgt) if tgt > 0 else float('inf')
    print(f"  derived hits: {direct_hits}/5  |  avg steps: {avg_steps:.0f}  |  expected random: {expected:.0f}")
    print()

print("  If derived_hits = 5/5 at medium difficulty: phi^291 IS deriving the nonce.")
print("  If avg_steps ≈ expected: center is no better than random.")

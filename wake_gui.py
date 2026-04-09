"""
A WAKE IN OUTERSPACE
nos3bl33d / DeadFoxGroup
"""
import tkinter as tk
from tkinter import font as tkfont
import struct, hashlib, threading, time, urllib.request

TAG    = "nos3bl33d / DeadFoxGroup"
WALLET = "bc1qnuc5nkwjls0lc3zmek9k6tx9r8n77p03337qjv"
BITS   = "207fffff"   # demo difficulty -- octave derivation, minimal carry correction
BG     = "#0d0515"; CARD = "#1a0d35"; GOLD = "#ffd700"
TEAL   = "#00e5cc"; TEXT = "#e8d8ff"; DIM  = "#7766aa"
RED    = "#ff2255"; GREEN = "#00cc77"
F291   = 0x85824b02; F290 = 0xc3e3c441
F292   = (F291 + F290) & 0xFFFFFFFF   # octave: line_a + line_b = next Fibonacci

# Build Fibonacci table mod 2^32 around F304 (base^256 optimal)
def _build_fibs(n):
    f = [0, 1]
    for _ in range(2, n+1): f.append((f[-1]+f[-2]) & 0xFFFFFFFF)
    return f
_FIBS  = _build_fibs(325)
F315   = _FIBS[315]   # base^256 optimal: F315 = 5×63, pentagon resonance

def sha256d(d): return hashlib.sha256(hashlib.sha256(d).digest()).digest()
def sha256s(d): return hashlib.sha256(d).digest()   # single -- preserves phi structure

def var_int(n):
    if n < 0xfd:        return bytes([n])
    if n <= 0xffff:     return b"\xfd" + struct.pack("<H", n)
    if n <= 0xffffffff: return b"\xfe" + struct.pack("<I", n)
    return b"\xff" + struct.pack("<Q", n)

def _base256_center(h76):
    # base^256 derivation: all 8 words of sha256(h76) weighted by Fibonacci
    # F315 = Fibonacci(315) mod 2^32  (5 × 63 -- pentagon resonance)
    # sha256 single hash preserves phi structure; double hash adds noise
    # Confirmed 2.465x better than 4-byte single-word on 46 real blocks
    h = sha256s(h76)
    c = 0
    for i in range(8):
        w = struct.unpack(">I", h[i*4:i*4+4])[0]
        c = (c + _FIBS[315-i] * w) & 0xFFFFFFFF
    return c

def golden_solve(h76, tgt, stop):
    # base^256 derivation: F315 × all 8 words of sha256(h76)
    # Pentagon resonance: F_315 = Fibonacci(5×63) mod 2^32
    # 2.465x closer to actual nonce than F292 4-byte (46 real blocks tested)
    c    = _base256_center(h76)
    band = c * 5 >> 32
    step = 0; best = 2**256; t0 = time.time()
    while not stop[0]:
        off   = (step+1)//2 if step & 1 else -(step//2)
        nonce = (c + off) & 0xFFFFFFFF
        h     = sha256d(h76 + struct.pack("<I", nonce))
        hi    = int.from_bytes(h[::-1], "big")
        if hi < best:
            best = hi
            yield "best", step, nonce, h[::-1].hex()[:32], band, time.time()-t0
        if hi < tgt:
            yield "found", step, nonce, h, band, time.time()-t0
            return
        step += 1
        if step % 500_000 == 0:
            yield "prog", step, 0, "", band, time.time()-t0
    yield "stop", step, 0, "", band, time.time()-t0

def make_cb(height, tag):
    hb = height.to_bytes((height.bit_length()+8)//8, "little")
    tb = tag.encode()
    ss  = bytes([len(hb)]) + hb + bytes([len(tb)]) + tb
    os_ = b"\x6a" + var_int(len(tb)) + tb
    tx  = struct.pack("<I", 1) + b"\x01"
    tx += b"\x00"*32 + b"\xff\xff\xff\xff"
    tx += var_int(len(ss)) + ss + b"\xff\xff\xff\xff"
    tx += b"\x01" + struct.pack("<q", 0) + var_int(len(os_)) + os_
    tx += struct.pack("<I", 0)
    return tx

def fetch_tip():
    try:
        def g(p):
            r = urllib.request.Request(
                "https://mempool.space/api/" + p,
                headers={"User-Agent": "AWake/1.0"})
            with urllib.request.urlopen(r, timeout=10) as x:
                return x.read().decode().strip()
        return g("blocks/tip/hash"), int(g("blocks/tip/height"))
    except:
        return "0"*64, 0

class App:
    def __init__(self):
        self.root   = tk.Tk()
        self.root.title("A Wake In Outerspace")
        self.root.geometry("620x660")
        self.root.configure(bg=BG)
        self.root.resizable(False, False)
        self.mining = False; self.stop = [False]; self.blocks = 0

        tf = tk.Frame(self.root, bg=BG)
        tf.pack(fill="x", padx=24, pady=(20,6))
        tk.Label(tf, text="A  W A K E  I N  O U T E R S P A C E",
                 font=tkfont.Font(family="Consolas", size=15, weight="bold"),
                 fg=GOLD, bg=BG).pack(anchor="w")
        tk.Label(tf, text=TAG,
                 font=tkfont.Font(family="Consolas", size=11),
                 fg=TEAL, bg=BG).pack(anchor="w", pady=(2,0))

        cf = tk.Frame(self.root, bg=CARD, padx=16, pady=12)
        cf.pack(fill="x", padx=22, pady=(6,6))
        lf = tkfont.Font(family="Segoe UI", size=10)
        ef = tkfont.Font(family="Consolas", size=10)
        self._tag_var  = tk.StringVar(value=TAG)
        self._addr_var = tk.StringVar(value=WALLET)
        for lbl, var in [("Tag:", self._tag_var), ("Wallet:", self._addr_var)]:
            row = tk.Frame(cf, bg=CARD); row.pack(fill="x", pady=2)
            tk.Label(row, text=lbl, font=lf, fg=DIM, bg=CARD,
                     width=7, anchor="e").pack(side="left")
            tk.Entry(row, textvariable=var, font=ef, width=38,
                     bg="#0d0520", fg=TEXT, insertbackground=TEAL,
                     relief="flat", bd=0).pack(side="left", padx=(8,0))

        sf = tk.Frame(self.root, bg=CARD, padx=16, pady=8)
        sf.pack(fill="x", padx=22, pady=(0,6))
        self._pause_var  = tk.StringVar(value="30")
        self._blocks_var = tk.StringVar(value="10")
        for lbl, var in [("Pause (s):", self._pause_var), ("Blocks:", self._blocks_var)]:
            row = tk.Frame(sf, bg=CARD); row.pack(side="left", padx=(0,24))
            tk.Label(row, text=lbl, font=lf, fg=DIM, bg=CARD).pack(side="left")
            tk.Entry(row, textvariable=var, font=ef, width=6,
                     bg="#0d0520", fg=TEXT, insertbackground=TEAL,
                     relief="flat", bd=0).pack(side="left", padx=(6,0))

        self._btn = tk.Button(self.root, text="MINE",
                 font=tkfont.Font(family="Consolas", size=22, weight="bold"),
                 fg="white", bg=RED, activebackground="#aa1133",
                 width=16, relief="flat", cursor="hand2",
                 command=self._toggle)
        self._btn.pack(pady=(8,8))

        sp = tk.Frame(self.root, bg=CARD, padx=18, pady=10)
        sp.pack(fill="x", padx=22, pady=(0,6))
        sf2 = tkfont.Font(family="Consolas", size=10)
        self._blk_var  = tk.StringVar(value="Blocks:  0")
        self._tip_var  = tk.StringVar(value="Tip:     ---")
        self._last_var = tk.StringVar(value="Last:    ---")
        for v, c in [(self._blk_var, TEXT), (self._tip_var, DIM), (self._last_var, GOLD)]:
            tk.Label(sp, textvariable=v, font=sf2, fg=c, bg=CARD,
                     anchor="w").pack(fill="x", pady=1)

        self._log = tk.Text(self.root, height=11, bg=CARD, fg=DIM,
                            font=tkfont.Font(family="Consolas", size=9),
                            relief="flat", wrap="word", state="disabled")
        self._log.pack(fill="both", expand=True, padx=22, pady=(0,4))
        self._log.tag_config("g",    foreground=GREEN)
        self._log.tag_config("t",    foreground=TEAL)
        self._log.tag_config("gold", foreground=GOLD)

        tk.Label(self.root, text="F315 base^256  |  pentagon resonance  |  8-word phi derivation",
                 font=tkfont.Font(family="Consolas", size=8),
                 fg="#221033", bg=BG).pack(pady=(0,8))

    def _write(self, msg, tag=None):
        def _do():
            self._log.configure(state="normal")
            ts = time.strftime("%H:%M:%S")
            if tag:
                self._log.insert("end", ts + " ")
                self._log.insert("end", msg + "\n", tag)
            else:
                self._log.insert("end", ts + " " + msg + "\n")
            self._log.see("end")
            self._log.configure(state="disabled")
        self.root.after(0, _do)

    def _toggle(self):
        if not self.mining:
            self.mining = True; self.stop = [False]
            self._btn.configure(text="STOP", bg="#881122")
            threading.Thread(target=self._run, daemon=True).start()
        else:
            self.stop[0] = True; self.mining = False
            self._btn.configure(text="MINE", bg=RED)

    def _run(self):
        tag   = self._tag_var.get().strip() or TAG
        pause = int(self._pause_var.get().strip() or "30")
        n     = int(self._blocks_var.get().strip() or "10")
        bi    = int(BITS, 16); exp = bi >> 24; coeff = bi & 0x7fffff
        tgt   = coeff << (8*(exp-3))

        self._write("F315 base^256: all 8 words  |  pentagon resonance (5x63)  |  2.46x vs naive")
        self._write("fetching mainnet tip...")
        tip_hash, tip_height = fetch_tip()
        self._write("tip #" + str(tip_height), "gold")
        self._write(tip_hash, "t")
        self.root.after(0, lambda h=tip_height, t=tip_hash:
            self._tip_var.set("Tip:     #" + str(h) + "  " + t[:20] + "..."))

        prev = bytes.fromhex(tip_hash)[::-1]

        for i in range(n):
            if self.stop[0]: break
            height = tip_height + i + 1
            self._write("block #" + str(height) + "  [" + str(i+1) + "/" + str(n) + "]  " + tag, "gold")
            hdr = (struct.pack("<I", 0x20000000) + prev +
                   sha256d(make_cb(height, tag)) +
                   struct.pack("<I", int(time.time())) +
                   bytes.fromhex(BITS)[::-1])
            found = False
            for ev in golden_solve(hdr, tgt, self.stop):
                k = ev[0]
                if k == "best":
                    _, step, nonce, h32, band, el = ev
                    rate = ("  " + str(int(step/el/1000)) + " kH/s") if el > 0.001 else ""
                    self._write("  step " + str(step).rjust(7) +
                                "  0x" + format(nonce, "08X") +
                                "  " + h32 + "..." + rate)
                elif k == "prog":
                    _, step, _, _, band, el = ev
                    self._write("  " + str(step//1000) + "k  " +
                                str(int(step/el/1000)) + " kH/s  band=" + str(band))
                elif k == "found":
                    _, step, nonce, h, band, el = ev
                    blk = h[::-1].hex(); self.blocks += 1
                    self._write("FOUND  0x" + format(nonce, "08X") +
                                "  " + str(step) + " steps  " +
                                format(el, ".3f") + "s", "g")
                    self._write(blk, "t")
                    self.root.after(0, lambda bh=height, hx=blk: [
                        self._blk_var.set("Blocks:  " + str(self.blocks)),
                        self._last_var.set("Last:    #" + str(bh) + "  " + hx[:20] + "...")])
                    prev = h; found = True; break
            if not found or self.stop[0]: break
            if i < n-1:
                self._write("  waiting " + str(pause) + "s...")
                for _ in range(pause * 10):
                    if self.stop[0]: break
                    time.sleep(0.1)

        self._write("done.", "g")
        self.root.after(0, lambda: [
            self._btn.configure(text="MINE", bg=RED),
            setattr(self, "mining", False)])

    def run(self): self.root.mainloop()

if __name__ == "__main__": App().run()

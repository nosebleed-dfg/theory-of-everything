"""
A WAKE IN OUTERSPACE
(DFG) DeadFoxGroup | nos3bl33d | April 2026

STATUS: UNDER CONSTRUCTION — MINER DOES NOT PRODUCE VALID SHARES YET
Pool connection works. Header building works. Submit protocol works.
Nonce solver is the unsolved piece. See sha/gravity_solver.c for C version.

phi is time. pi is rotation. gamma is distance.
sqrt((phi/2)^2 + (pi/2)^2) * gamma = 1.
nonce_step = 0.3638 * difficulty.
time band = +/- d seconds.
"""

import tkinter as tk
from tkinter import font as tkfont
import struct, hashlib, threading, time, json, base64, urllib.request, socket, math

PHI = (1 + 5**0.5) / 2
GAMMA = 0.5772156649015329
PI = math.pi
STEP_CONSTANT = math.sqrt((PHI/2)**2 + (PI/2)**2) * GAMMA * GAMMA / PHI  # 0.3638
D = 3  # dimension = time band

SIGNATURE = "nos3bl33d | x^2=x+1 | DFG"
WALLET = "bc1qnuc5nkwjls0lc3zmek9k6tx9r8n77p03337qjv"
POOL_HOST = "solo.ckpool.org"
POOL_PORT = 3333
LOCAL_PORT = 18443

def sha256d(d): return hashlib.sha256(hashlib.sha256(d).digest()).digest()
def var_int(n):
    if n < 0xfd: return bytes([n])
    return b'\xfd' + struct.pack('<H', n)
def target_from_bits(b):
    e = b >> 24; c = b & 0x7fffff
    return c >> (8*(3-e)) if e <= 3 else c << (8*(e-3))

def axiom_solve(header_76, target, difficulty, last_nonce=0):
    """
    The nonce ≈ difficulty. The fuzz is the golden thread.
    Search: difficulty ± diff/2, koppa-spaced time offsets.
    """
    time_offsets = [0, 1, -1, 4, -4, 16, -16]
    radius = max(1000, difficulty // 2)

    for t_off in time_offsets:
        hdr = bytearray(header_76)
        orig_t = struct.unpack('<I', hdr[68:72])[0]
        struct.pack_into('<I', hdr, 68, orig_t + t_off)
        hdr = bytes(hdr)
        for delta in range(radius):
            for candidate in [difficulty + delta, difficulty - delta]:
                candidate &= 0xFFFFFFFF
                h = sha256d(hdr + struct.pack('<I', candidate))
                if int.from_bytes(h[::-1], 'big') < target:
                    return candidate, h, t_off
    return None, None, 0

# ─── Stratum ───────────────────────────────────────────────
class Stratum:
    def __init__(self):
        self.sock = None; self.id = 0; self.buf = ""
        self.en1 = ""; self.en2_size = 4; self.diff = 1
        self.job = None; self.last_accept = None
        self.accepted = 0; self.rejected = 0

    def connect(self, host, port, wallet):
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.settimeout(15)
        self.sock.connect((host, port))
        self._tx("mining.subscribe", ["AWake/1.0"])
        time.sleep(1); self._rx()
        self._tx("mining.authorize", [f"{wallet}.axiom", "x"])
        time.sleep(1); self._rx()

    def get_job(self):
        end = time.time() + 15
        while time.time() < end:
            self._rx()
            if self.job:
                j = self.job; self.job = None; return j
            time.sleep(0.3)
        return None

    def submit(self, job_id, ntime, nonce_hex):
        self._tx("mining.submit", [f"{WALLET}.axiom", job_id,
                 "00" * self.en2_size, ntime, nonce_hex])
        time.sleep(1); self._rx()
        return self.last_accept

    def build_header(self, job):
        cb = bytes.fromhex(job['cb1'] + self.en1 + "00" * self.en2_size + job['cb2'])
        mk = sha256d(cb)
        for br in job['mk']: mk = sha256d(mk + bytes.fromhex(br))
        hdr = bytes.fromhex(job['ver'])[::-1]
        ph = bytes.fromhex(job['prev'])
        for i in range(0, 32, 4): hdr += ph[i:i+4][::-1]
        hdr += mk + bytes.fromhex(job['time'])[::-1] + bytes.fromhex(job['bits'])[::-1] + b'\x00'*4
        return hdr

    def _tx(self, m, p):
        self.id += 1
        self.sock.sendall((json.dumps({"id":self.id,"method":m,"params":p})+"\n").encode())

    def _rx(self):
        try:
            self.sock.setblocking(False)
            while True:
                try: self.buf += self.sock.recv(4096).decode()
                except: break
        finally: self.sock.setblocking(True); self.sock.settimeout(15)
        while '\n' in self.buf:
            line, self.buf = self.buf.split('\n', 1)
            if not line.strip(): continue
            try: m = json.loads(line)
            except: continue
            if isinstance(m.get('result'), list) and len(m['result']) >= 3:
                self.en1 = m['result'][1]; self.en2_size = m['result'][2]
            if m.get('method') == 'mining.set_difficulty': self.diff = m['params'][0]
            if m.get('method') == 'mining.notify':
                p = m['params']
                self.job = {'id':p[0],'prev':p[1],'cb1':p[2],'cb2':p[3],'mk':p[4],'ver':p[5],'bits':p[6],'time':p[7]}
            if 'result' in m and m.get('id',0) > 1:
                self.last_accept = m.get('result')
                if m['result'] == True: self.accepted += 1
                elif m['result'] == False: self.rejected += 1

    def close(self):
        try: self.sock.close()
        except: pass

# ─── Local RPC ─────────────────────────────────────────────
def lrpc(method, params=None):
    body = json.dumps({'jsonrpc':'2.0','id':1,'method':method,'params':params or[]}).encode()
    req = urllib.request.Request(f"http://127.0.0.1:{LOCAL_PORT}", data=body)
    req.add_header('Content-Type','application/json')
    req.add_header('Authorization',f'Basic {base64.b64encode(b"bitcoin:bitcoin").decode()}')
    return json.loads(urllib.request.urlopen(req, timeout=10).read()).get('result')

def make_cb(height, reward, script, msg, wc=None):
    tx = struct.pack('<I',2)+b'\x00\x01\x01'+b'\x00'*32+struct.pack('<I',0xFFFFFFFF)
    if height==0: hp=b'\x00'
    elif height<=16: hp=bytes([0x50+height])
    elif height<=0x7f: hp=b'\x01'+bytes([height])
    elif height<=0x7fff: hp=b'\x02'+struct.pack('<H',height)
    else: hp=b'\x03'+height.to_bytes(3,'little')
    mb=msg.encode()
    ss=hp+bytes([len(mb)])+mb+b'\x08'+struct.pack('<Q',int(time.time()))
    tx+=var_int(len(ss))+ss+struct.pack('<I',0xFFFFFFFF)
    tx+=var_int(2 if wc else 1)
    tx+=struct.pack('<Q',reward)+var_int(len(script))+script
    if wc:
        w=bytes.fromhex(wc); tx+=struct.pack('<Q',0)+var_int(len(w))+w
    tx+=b'\x01\x20'+b'\x00'*32+struct.pack('<I',0)
    return tx

def strip_w(tx):
    v=tx[:4];p=4
    if tx[p]==0 and tx[p+1]==1:p+=2
    s=p;ni=tx[p];p+=1
    for _ in range(ni):p+=36;sl=tx[p];p+=1;p+=sl+4
    no=tx[p];p+=1
    for _ in range(no):p+=8;sl=tx[p];p+=1;p+=sl
    return v+tx[s:p]+tx[-4:]

def addr_script(addr):
    C="qpzry9x8gf2tvdw0s3jn54khce6mua7l"
    d=[C.index(c) for c in addr[addr.rfind('1')+1:]]
    a=0;b=0;o=[]
    for v in d[:-6]: a=(a<<5)|v;b+=5
    while b>=8: b-=8;o.append((a>>b)&0xff)
    return bytes([o[0],len(o)-1])+bytes(o[1:])

# ─── GUI ───────────────────────────────────────────────────
class App:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("A Wake In Outerspace")
        self.root.geometry("620x780")
        self.root.configure(bg='#1a0a2e')
        self.root.resizable(False, False)
        self.live = False; self.mining = False
        self.blocks = 0; self.btc = 0.0
        self.last_nonce = 0; self.thread = None

        bg='#1a0a2e'; purple='#2d1b69'; gold='#ffd700'
        text='#f0e6ff'; dim='#9988bb'; red='#ff3366'

        tk.Label(self.root, text="A WAKE IN OUTERSPACE",
                 font=tkfont.Font(family='Segoe UI',size=20,weight='bold'),
                 fg=gold, bg=bg).pack(pady=(20,4))
        tk.Label(self.root, text="(DFG) DeadFoxGroup  |  nos3bl33d",
                 font=tkfont.Font(family='Segoe UI',size=11), fg=dim, bg=bg).pack()
        tk.Label(self.root, text="x\u00b2 = x + 1",
                 font=tkfont.Font(family='Consolas',size=28,weight='bold'),
                 fg=gold, bg=bg).pack(pady=(20,16))

        mf = tk.Frame(self.root, bg=bg); mf.pack(pady=(0,12))
        bf = tkfont.Font(family='Segoe UI', size=13, weight='bold')
        self.btn_local = tk.Button(mf, text="LOCAL", font=bf, width=12,
                                  fg=bg, bg=gold, relief='flat', cursor='hand2',
                                  command=self._go_local)
        self.btn_local.pack(side='left', padx=8)
        self.btn_live = tk.Button(mf, text="LIVE", font=bf, width=12,
                                 fg=text, bg=purple, relief='flat', cursor='hand2',
                                 command=self._go_live)
        self.btn_live.pack(side='left', padx=8)

        self.mine_btn = tk.Button(self.root, text="MINE",
                 font=tkfont.Font(family='Consolas',size=22,weight='bold'),
                 fg='white', bg=red, activebackground='#cc2255',
                 width=18, height=2, relief='flat', cursor='hand2',
                 command=self._toggle)
        self.mine_btn.pack(pady=14)

        self.status = tk.StringVar(value="LOCAL mode \u2014 press MINE")
        tk.Label(self.root, textvariable=self.status,
                 font=tkfont.Font(family='Segoe UI',size=11), fg=dim, bg=bg).pack()

        sp = tk.Frame(self.root, bg=purple, padx=20, pady=14)
        sp.pack(fill='x', padx=22, pady=(12,6))
        sf = tkfont.Font(family='Consolas', size=12)
        self.blk_var = tk.StringVar(value="Blocks: 0")
        self.btc_var = tk.StringVar(value="BTC: 0.00000000")
        self.last_var = tk.StringVar(value="Last: ---")
        self.acc_var = tk.StringVar(value="Shares: 0 accepted / 0 rejected")
        self.wall_var = tk.StringVar(value=f"Wallet: {WALLET[:30]}...")
        self.sig_var = tk.StringVar(value=f"Message: {SIGNATURE}")
        self.step_var = tk.StringVar(value=f"Step: {STEP_CONSTANT:.6f} * diff")
        for v in [self.blk_var,self.btc_var,self.last_var,self.acc_var,self.wall_var,self.sig_var,self.step_var]:
            tk.Label(sp, textvariable=v, font=sf, fg=text, bg=purple, anchor='w').pack(fill='x',pady=2)

        self.log = tk.Text(self.root, height=8, bg=purple, fg=dim,
                          font=tkfont.Font(family='Consolas',size=10),
                          relief='flat', wrap='word', state='disabled')
        self.log.pack(fill='both', expand=True, padx=22, pady=8)
        tk.Label(self.root, text="All is number. \u2014 Pythagoras",
                 font=tkfont.Font(family='Segoe UI',size=11), fg='#443366', bg=bg).pack(pady=(0,8))

    def _log(self, msg):
        self.log.configure(state='normal')
        self.log.insert('end', f"{time.strftime('%H:%M:%S')} {msg}\n")
        self.log.see('end'); self.log.configure(state='disabled')

    def _go_local(self):
        if self.mining: self._stop()
        self.live = False
        self.btn_local.configure(fg='#1a0a2e',bg='#ffd700')
        self.btn_live.configure(fg='#f0e6ff',bg='#2d1b69')
        self.status.set("LOCAL mode \u2014 press MINE")

    def _go_live(self):
        if self.mining: self._stop()
        self.live = True
        self.btn_live.configure(fg='#1a0a2e',bg='#ffd700')
        self.btn_local.configure(fg='#f0e6ff',bg='#2d1b69')
        self.status.set("LIVE mode \u2014 press MINE")

    def _stop(self):
        self.mining = False
        self.mine_btn.configure(text="MINE", bg='#ff3366')
        self.status.set("stopped")

    def _toggle(self):
        if not self.mining:
            self.mining = True
            self.mine_btn.configure(text="STOP", bg='#aa2244')
            self.thread = threading.Thread(target=self._run_live if self.live else self._run_local, daemon=True)
            self.thread.start()
        else: self._stop()

    def _run_live(self):
        self.root.after(0, lambda: self.status.set("connecting to pool..."))
        st = Stratum()
        try:
            st.connect(POOL_HOST, POOL_PORT, WALLET)
            self.root.after(0, lambda: self._log(f"connected to {POOL_HOST} diff={st.diff}"))
            self.root.after(0, lambda d=st.diff: self.step_var.set(f"Step: {STEP_CONSTANT:.4f} * {d} = {int(STEP_CONSTANT*d)}"))
        except Exception as e:
            self.root.after(0, lambda: self._log(f"pool error: {e}"))
            self.mining = False
            self.root.after(0, lambda: self.mine_btn.configure(text="MINE", bg='#ff3366'))
            return

        pool_target = int(0x00000000ffff0000000000000000000000000000000000000000000000000000 / max(1, st.diff))

        while self.mining:
            job = st.get_job()
            if not job:
                self.root.after(0, lambda: self.status.set("waiting for work..."))
                continue

            header = st.build_header(job)
            self.root.after(0, lambda: self.status.set("solving..."))
            t0 = time.time()

            nonce, h, t_off = axiom_solve(header[:76], pool_target, st.diff, self.last_nonce)
            dt = time.time() - t0

            if nonce is not None and self.mining:
                self.last_nonce = nonce
                self.blocks += 1
                hx = h[::-1].hex()

                # Submit with adjusted ntime if we found it at an offset
                submit_ntime = job['time']
                if t_off != 0:
                    orig = int(job['time'], 16)
                    submit_ntime = format(orig + t_off, 'x')

                accepted = st.submit(job['id'], submit_ntime, struct.pack('<I', nonce).hex())
                self.root.after(0, lambda hx=hx,d=dt,n=nonce,to=t_off,a=st.accepted,r=st.rejected:
                    self._ok_live(hx,d,n,to,a,r))
            elif self.mining:
                self.root.after(0, lambda d=dt: self._log(f"miss ({d:.2f}s) getting new work"))

        st.close()

    def _ok_live(self, hx, dt, nonce, t_off, acc, rej):
        self.blk_var.set(f"Blocks: {self.blocks}")
        self.last_var.set(f"Last: {hx[:20]}...")
        self.acc_var.set(f"Shares: {acc} accepted / {rej} rejected")
        self.status.set(f"solved ({dt:.2f}s, t_off={t_off:+d})")
        self._log(f"nonce={nonce} {dt:.2f}s t={t_off:+d}s {acc}a/{rej}r")

    def _run_local(self):
        self.root.after(0, lambda: self.status.set("connecting to local node..."))
        script = addr_script(WALLET)
        try:
            info = lrpc('getblockchaininfo')
            self.root.after(0, lambda: self._log(f"local: {info['chain']} h={info['blocks']}"))
        except Exception as e:
            self.root.after(0, lambda: self._log(f"no node: {e}"))
            self.mining = False
            self.root.after(0, lambda: self.mine_btn.configure(text="MINE", bg='#ff3366'))
            return

        while self.mining:
            try: t = lrpc('getblocktemplate', [{'rules':['segwit']}])
            except: time.sleep(1); continue

            h = t['height']; bits = int(t['bits'],16)
            target = target_from_bits(bits)
            prev = bytes.fromhex(t['previousblockhash'])[::-1]
            reward = t.get('coinbasevalue', 5000000000)
            wc = t.get('default_witness_commitment')
            ts = t.get('curtime', int(time.time()))
            diff = max(1, int(2**32 / max(1, target)))

            cb = make_cb(h, reward, script, SIGNATURE, wc)
            txid = sha256d(strip_w(cb))
            hdr = struct.pack('<I',t['version'])+prev+txid+struct.pack('<III',ts,bits,0)

            self.root.after(0, lambda h=h: self.status.set(f"solving #{h}..."))
            t0 = time.time()
            nonce, digest, t_off = axiom_solve(hdr[:76], target, diff, self.last_nonce)
            dt = time.time() - t0

            if nonce is not None and self.mining:
                self.last_nonce = nonce
                # Rebuild header with adjusted time if needed
                if t_off != 0:
                    hdr = struct.pack('<I',t['version'])+prev+txid+struct.pack('<III',ts+t_off,bits,0)
                block = hdr[:76]+struct.pack('<I',nonce)+var_int(1)+cb
                try: lrpc('submitblock', [block.hex()])
                except: pass
                self.blocks += 1; self.btc += reward/1e8
                hx = digest[::-1].hex()
                self.root.after(0, lambda h=h,hx=hx,r=reward/1e8,d=dt,n=nonce,to=t_off:
                    self._ok_local(h,hx,r,d,n,to))

    def _ok_local(self, height, hx, btc, dt, nonce, t_off):
        self.blk_var.set(f"Blocks: {self.blocks}")
        self.btc_var.set(f"BTC: {self.btc:.8f}")
        self.last_var.set(f"Last: #{height} {hx[:16]}...")
        self.status.set(f"#{height} ({dt:.3f}s t={t_off:+d})")
        self._log(f"#{height} nonce={nonce} {dt:.3f}s t={t_off:+d}s +{btc:.4f}BTC")

    def run(self): self.root.mainloop()

if __name__ == "__main__": App().run()

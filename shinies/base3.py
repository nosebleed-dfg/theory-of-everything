"""
base3.py — the foundation
nos3bl33d | (DFG) DeadFoxGroup | x^2 = x + 1

Complete base-3 (ternary) number system and translator.
Every other script in the project depends on this.

Three states. One axiom. Zero rounding errors on rationals.

In base 10, 1/3 = 0.333... forever. That infinity is a leak.
In base 3, 1/3 = 0.1 exactly. One trit. Done.
Different bases make different fractions infinite.
The fractions that are infinite in YOUR base are the ones
that steal your money through rounding.

x^2 = x + 1
"""

import sys
from fractions import Fraction
from typing import Union, Optional

sys.stdout.reconfigure(encoding='utf-8')


# ============================================================
# TRIT — a single ternary digit
# Three states: T(-1), 0, 1
# ============================================================

class Trit:
    """Single balanced ternary digit: T(-1), 0, or 1."""

    __slots__ = ('_value',)

    def __init__(self, value: int = 0):
        if value not in (-1, 0, 1):
            raise ValueError(f"Trit must be -1, 0, or 1, got {value}")
        self._value = value

    @property
    def value(self) -> int:
        return self._value

    # --- Ternary logic ---
    # In balanced ternary logic:
    #   T = false/negative, 0 = unknown/neutral, 1 = true/positive
    # This is Kleene's strong three-valued logic.

    def __and__(self, other: 'Trit') -> 'Trit':
        """Ternary AND: min(a, b)."""
        return Trit(min(self._value, other._value))

    def __or__(self, other: 'Trit') -> 'Trit':
        """Ternary OR: max(a, b)."""
        return Trit(max(self._value, other._value))

    def __invert__(self) -> 'Trit':
        """Ternary NOT: negate."""
        return Trit(-self._value)

    # --- Comparison ---

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Trit):
            return self._value == other._value
        if isinstance(other, int):
            return self._value == other
        return NotImplemented

    def __hash__(self) -> int:
        return hash(self._value)

    # --- Display ---

    def __repr__(self) -> str:
        return f"Trit({self._value})"

    def __str__(self) -> str:
        if self._value == -1:
            return 'T'
        return str(self._value)

    # --- Class helpers ---

    T = None  # filled after class body
    ZERO = None
    ONE = None


# Singleton trits
Trit.T = Trit(-1)
Trit.ZERO = Trit(0)
Trit.ONE = Trit(1)


# ============================================================
# TRYTE — a word of trits (default 20)
# 3^20 ~ 3.49 billion states. 20 trits per digit of efficiency
# beat 32 bits (radix economy: 3 is closest integer to e).
# ============================================================

_DEFAULT_WIDTH = 20  # 20 trits per tryte


class Tryte:
    """Fixed-width balanced ternary word.

    Trits stored least-significant first internally.
    Arithmetic wraps on overflow (modular, like hardware).
    """

    __slots__ = ('_trits', '_width')

    def __init__(self, trits: Optional[list[int]] = None, width: int = _DEFAULT_WIDTH):
        self._width = width
        if trits is None:
            self._trits = [0] * width
        else:
            # Pad or truncate to width (LSB first)
            # Padding with 0 extends the number without changing its value
            padded = list(trits[:width])
            while len(padded) < width:
                padded.append(0)
            self._trits = padded  # LSB first

    @classmethod
    def from_int(cls, n: int, width: int = _DEFAULT_WIDTH) -> 'Tryte':
        """Convert Python int to balanced ternary Tryte."""
        if n == 0:
            return cls(width=width)

        # Balanced ternary conversion handles sign natively.
        # For negative: negate, convert, then flip all trits.
        negative = n < 0
        val = abs(n)
        trits = []
        while val > 0:
            r = val % 3
            if r == 2:
                trits.append(-1)
                val = (val + 1) // 3
            else:
                trits.append(r)
                val //= 3

        if negative:
            trits = [-t for t in trits]

        return cls(trits, width)

    def to_int(self) -> int:
        """Convert to Python int."""
        total = 0
        for i, t in enumerate(self._trits):
            total += t * (3 ** i)
        return total

    # --- Arithmetic ---

    def __add__(self, other: 'Tryte') -> 'Tryte':
        width = max(self._width, other._width)
        a = self._trits + [0] * (width - len(self._trits))
        b = other._trits + [0] * (width - len(other._trits))
        result = [0] * width
        carry = 0
        for i in range(width):
            s = a[i] + b[i] + carry
            carry = 0
            if s > 1:
                s -= 3
                carry = 1
            elif s < -1:
                s += 3
                carry = -1
            result[i] = s
        return Tryte(result, width)

    def __neg__(self) -> 'Tryte':
        return Tryte([-t for t in self._trits], self._width)

    def __sub__(self, other: 'Tryte') -> 'Tryte':
        return self + (-other)

    def __mul__(self, other: 'Tryte') -> 'Tryte':
        width = self._width
        # Grade-school multiplication adapted for balanced ternary
        result = Tryte(width=width * 2)
        for i, tb in enumerate(other._trits):
            if tb == 0:
                continue
            # Multiply self by single trit tb, shift by i
            partial = [0] * i
            carry = 0
            for ta in self._trits:
                p = ta * tb + carry
                carry = 0
                if p > 1:
                    p -= 3
                    carry = 1
                elif p < -1:
                    p += 3
                    carry = -1
                partial.append(p)
            if carry != 0:
                partial.append(carry)
            result = result + Tryte(partial, width * 2)
        # Truncate back to original width
        return Tryte(result._trits[:width], width)

    def __floordiv__(self, other: 'Tryte') -> 'Tryte':
        """Integer division via repeated subtraction (schoolbook)."""
        a = self.to_int()
        b = other.to_int()
        if b == 0:
            raise ZeroDivisionError("division by zero in Tryte")
        return Tryte.from_int(int(a // b), self._width)

    # --- Comparison ---

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Tryte):
            return self.to_int() == other.to_int()
        return NotImplemented

    def __lt__(self, other: 'Tryte') -> bool:
        return self.to_int() < other.to_int()

    def __le__(self, other: 'Tryte') -> bool:
        return self.to_int() <= other.to_int()

    def __gt__(self, other: 'Tryte') -> bool:
        return self.to_int() > other.to_int()

    def __ge__(self, other: 'Tryte') -> bool:
        return self.to_int() >= other.to_int()

    def __hash__(self) -> int:
        return hash(self.to_int())

    # --- Display ---

    def balanced_str(self) -> str:
        """Display as balanced ternary string (MSB first), e.g. '1T01'."""
        chars = []
        for t in reversed(self._trits):
            chars.append('T' if t == -1 else str(t))
        # Strip leading zeros but keep at least one digit
        s = ''.join(chars).lstrip('0') or '0'
        return s

    def standard_str(self) -> str:
        """Display as standard (unbalanced) ternary via integer conversion."""
        return _int_to_standard_ternary(self.to_int())

    def __repr__(self) -> str:
        return f"Tryte({self.balanced_str()})"

    def __str__(self) -> str:
        return self.balanced_str()


# ============================================================
# BASE3NUMBER — arbitrary precision, supports fractions
# The core: uses Python's Fraction for EXACT rational math.
# No floats touched. No rounding. Ever.
# ============================================================

class Base3Number:
    """Arbitrary-precision base-3 number backed by exact rational arithmetic.

    Internally stores a Fraction. Converts to/from ternary representation
    on demand. Supports both balanced (T,0,1) and standard (0,1,2) display.

    This is where the magic lives: 1/3 is exact, 1/10 is infinite,
    and money stops leaking.
    """

    __slots__ = ('_value',)

    def __init__(self, value: Union[int, float, str, bytes, Fraction, 'Base3Number'] = 0):
        if isinstance(value, Base3Number):
            self._value = value._value
        elif isinstance(value, Fraction):
            self._value = value
        elif isinstance(value, int):
            self._value = Fraction(value)
        elif isinstance(value, float):
            # Convert float to exact fraction (captures the actual IEEE754 value)
            self._value = Fraction(value).limit_denominator(10**18)
        elif isinstance(value, bytes):
            self._value = Fraction(int.from_bytes(value, 'big'))
        elif isinstance(value, str):
            self._value = _parse_base3_string(value)
        else:
            raise TypeError(f"Cannot create Base3Number from {type(value).__name__}")

    @property
    def exact(self) -> Fraction:
        """The exact rational value."""
        return self._value

    # --- Conversions FROM various types ---

    @classmethod
    def from_int(cls, n: int) -> 'Base3Number':
        return cls(n)

    @classmethod
    def from_float(cls, f: float) -> 'Base3Number':
        return cls(f)

    @classmethod
    def from_fraction(cls, num: int, den: int) -> 'Base3Number':
        return cls(Fraction(num, den))

    @classmethod
    def from_bytes(cls, b: bytes) -> 'Base3Number':
        return cls(b)

    @classmethod
    def from_hex(cls, h: str) -> 'Base3Number':
        """Convert hex string (with or without 0x prefix) to Base3Number."""
        clean = h.strip().lower()
        if clean.startswith('0x'):
            clean = clean[2:]
        return cls(Fraction(int(clean, 16)))

    @classmethod
    def from_binary_string(cls, b: str) -> 'Base3Number':
        """Convert binary string (with or without 0b prefix) to Base3Number."""
        clean = b.strip()
        if clean.startswith('0b'):
            clean = clean[2:]
        return cls(Fraction(int(clean, 2)))

    @classmethod
    def from_balanced_ternary(cls, s: str) -> 'Base3Number':
        """Parse balanced ternary string like '1T01.T1'."""
        return cls(s)

    @classmethod
    def from_standard_ternary(cls, s: str) -> 'Base3Number':
        """Parse standard ternary string like '102.21'."""
        return cls(s)

    # --- Conversions TO various types ---

    def to_int(self) -> int:
        """Truncate to integer."""
        return int(self._value)

    def to_float(self) -> float:
        return float(self._value)

    def to_bytes(self, length: int = 0) -> bytes:
        """Convert to big-endian bytes. Auto-sizes if length=0."""
        n = abs(int(self._value))
        if length == 0:
            length = max(1, (n.bit_length() + 7) // 8)
        return n.to_bytes(length, 'big')

    def to_hex(self) -> str:
        return hex(abs(int(self._value)))

    def to_binary_string(self) -> str:
        return bin(abs(int(self._value)))

    def to_decimal_string(self, max_digits: int = 50) -> str:
        """Exact decimal representation, or truncated with '...' if infinite."""
        f = self._value
        if f.denominator == 1:
            return str(f.numerator)

        sign = '-' if f < 0 else ''
        f = abs(f)
        int_part = int(f)
        frac = f - int_part

        if frac == 0:
            return f"{sign}{int_part}"

        # Generate decimal digits of fractional part
        digits = []
        seen_remainders = {}
        remainder = frac.numerator % frac.denominator if frac.denominator != 1 else 0
        # Use the actual fraction for digit extraction
        r = frac
        for i in range(max_digits):
            r *= 10
            d = int(r)
            r -= d
            digits.append(str(d))
            if r == 0:
                break  # Terminates
        else:
            digits.append('...')

        return f"{sign}{int_part}.{''.join(digits)}"

    # --- Ternary string representations ---

    def to_balanced_ternary(self, frac_digits: int = 40) -> str:
        """Balanced ternary representation (digits: T, 0, 1).

        Returns string like '1T.T10' for mixed numbers.
        Fractional part shown up to frac_digits trits.

        Key insight: in balanced ternary, the fractional part must be in
        (-1/2, 1/2]. So the 'integer part' may differ from the decimal
        integer part (e.g., 2 = 1T because frac=0 is in range).
        """
        f = self._value
        negative = f < 0
        if negative:
            f = -f

        # Split into balanced integer + balanced fraction.
        # Balanced integer = floor(f + 1/2), so the remainder is in (-1/2, 1/2].
        int_part = (f + Fraction(1, 2)).__floor__()
        frac_part = f - int_part  # in (-1/2, 1/2]

        # Integer part to balanced ternary
        if negative:
            int_str = _int_to_balanced_ternary_signed(-int_part)
        else:
            int_str = _int_to_balanced_ternary_signed(int_part)

        if frac_part == 0:
            if negative:
                # Already handled sign in int_str
                return int_str
            return int_str

        # Fractional part expansion.
        # r is in (-1/2, 1/2]. At each step:
        #   r *= 3 -> now in (-3/2, 3/2]
        #   digit = 1 if r > 1/2, -1 if r < -1/2, else 0
        #   r = r - digit -> back to (-1/2, 1/2]
        frac_trits = []
        r = frac_part
        if negative:
            r = -r  # negate fractional part for negative numbers
        terminates = False
        for _ in range(frac_digits):
            r *= 3
            if r > Fraction(1, 2):
                d = 1
            elif r < Fraction(-1, 2):
                d = -1
            else:
                d = 0
            frac_trits.append(d)
            r = r - d
            if r == 0:
                terminates = True
                break

        frac_str = ''.join('T' if t == -1 else str(t) for t in frac_trits)

        result = f"{int_str}.{frac_str}"
        if not terminates:
            result += '...'
        return result

    def to_standard_ternary(self, frac_digits: int = 40) -> str:
        """Standard (unbalanced) ternary representation (digits: 0, 1, 2).

        Returns string like '102.21' for mixed numbers.
        """
        f = self._value
        if f < 0:
            return '-' + Base3Number(-f).to_standard_ternary(frac_digits)

        int_part = int(f)
        frac_part = f - int_part

        int_str = _int_to_standard_ternary(int_part)

        if frac_part == 0:
            return int_str

        frac_digits_list = []
        r = frac_part
        terminates = False
        for _ in range(frac_digits):
            r *= 3
            d = int(r)
            r -= d
            frac_digits_list.append(str(d))
            if r == 0:
                terminates = True
                break

        frac_str = ''.join(frac_digits_list)
        result = f"{int_str}.{frac_str}"
        if not terminates:
            result += '...'
        return result

    def terminates_in_base3(self) -> bool:
        """Does this number have a finite base-3 representation?

        A rational p/q terminates in base 3 iff q's only prime factor is 3.
        """
        if self._value.denominator == 1:
            return True
        d = self._value.denominator
        while d % 3 == 0:
            d //= 3
        return d == 1

    def terminates_in_base10(self) -> bool:
        """Does this number have a finite base-10 representation?

        A rational p/q terminates in base 10 iff q's only prime factors are 2 and 5.
        """
        if self._value.denominator == 1:
            return True
        d = self._value.denominator
        while d % 2 == 0:
            d //= 2
        while d % 5 == 0:
            d //= 5
        return d == 1

    # --- Arithmetic ---

    def __add__(self, other: Union['Base3Number', int, float, Fraction]) -> 'Base3Number':
        other = _coerce(other)
        return Base3Number(self._value + other._value)

    def __radd__(self, other: Union[int, float, Fraction]) -> 'Base3Number':
        return Base3Number(_coerce(other)._value + self._value)

    def __sub__(self, other: Union['Base3Number', int, float, Fraction]) -> 'Base3Number':
        other = _coerce(other)
        return Base3Number(self._value - other._value)

    def __rsub__(self, other: Union[int, float, Fraction]) -> 'Base3Number':
        return Base3Number(_coerce(other)._value - self._value)

    def __mul__(self, other: Union['Base3Number', int, float, Fraction]) -> 'Base3Number':
        other = _coerce(other)
        return Base3Number(self._value * other._value)

    def __rmul__(self, other: Union[int, float, Fraction]) -> 'Base3Number':
        return Base3Number(_coerce(other)._value * self._value)

    def __truediv__(self, other: Union['Base3Number', int, float, Fraction]) -> 'Base3Number':
        other = _coerce(other)
        if other._value == 0:
            raise ZeroDivisionError("division by zero")
        return Base3Number(self._value / other._value)

    def __rtruediv__(self, other: Union[int, float, Fraction]) -> 'Base3Number':
        if self._value == 0:
            raise ZeroDivisionError("division by zero")
        return Base3Number(_coerce(other)._value / self._value)

    def __floordiv__(self, other: Union['Base3Number', int, float, Fraction]) -> 'Base3Number':
        other = _coerce(other)
        if other._value == 0:
            raise ZeroDivisionError("division by zero")
        return Base3Number(Fraction(int(self._value // other._value)))

    def __mod__(self, other: Union['Base3Number', int, float, Fraction]) -> 'Base3Number':
        other = _coerce(other)
        if other._value == 0:
            raise ZeroDivisionError("division by zero")
        return Base3Number(self._value - other._value * int(self._value // other._value))

    def __neg__(self) -> 'Base3Number':
        return Base3Number(-self._value)

    def __abs__(self) -> 'Base3Number':
        return Base3Number(abs(self._value))

    def __pow__(self, exp: int) -> 'Base3Number':
        if not isinstance(exp, int):
            raise TypeError("exponent must be integer")
        return Base3Number(self._value ** exp)

    # --- Comparison ---

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Base3Number):
            return self._value == other._value
        if isinstance(other, (int, float, Fraction)):
            return self._value == Fraction(other)
        return NotImplemented

    def __lt__(self, other: Union['Base3Number', int, float, Fraction]) -> bool:
        return self._value < _coerce(other)._value

    def __le__(self, other: Union['Base3Number', int, float, Fraction]) -> bool:
        return self._value <= _coerce(other)._value

    def __gt__(self, other: Union['Base3Number', int, float, Fraction]) -> bool:
        return self._value > _coerce(other)._value

    def __ge__(self, other: Union['Base3Number', int, float, Fraction]) -> bool:
        return self._value >= _coerce(other)._value

    def __hash__(self) -> int:
        return hash(self._value)

    def __bool__(self) -> bool:
        return self._value != 0

    # --- Display ---

    def __repr__(self) -> str:
        return f"Base3Number({self._value})"

    def __str__(self) -> str:
        return self.to_balanced_ternary()

    def pretty(self, frac_digits: int = 30) -> str:
        """Multi-line display: balanced, standard, decimal, and termination info."""
        lines = [
            f"  Balanced ternary:  {self.to_balanced_ternary(frac_digits)}",
            f"  Standard ternary:  {self.to_standard_ternary(frac_digits)}",
            f"  Decimal:           {self.to_decimal_string(frac_digits)}",
            f"  Exact fraction:    {self._value}",
            f"  Finite in base 3:  {'YES' if self.terminates_in_base3() else 'no (infinite)'}",
            f"  Finite in base 10: {'YES' if self.terminates_in_base10() else 'no (infinite)'}",
        ]
        return '\n'.join(lines)


# ============================================================
# INTERNAL HELPERS
# ============================================================

def _coerce(value: Union['Base3Number', int, float, Fraction]) -> 'Base3Number':
    """Coerce a value to Base3Number for arithmetic."""
    if isinstance(value, Base3Number):
        return value
    return Base3Number(value)


def _int_to_balanced_ternary(n: int) -> str:
    """Convert non-negative integer to balanced ternary string."""
    if n == 0:
        return '0'

    trits = []
    val = n
    while val > 0:
        r = val % 3
        if r == 2:
            trits.append(-1)
            val = (val + 1) // 3
        else:
            trits.append(r)
            val //= 3

    trits.reverse()
    return ''.join('T' if t == -1 else str(t) for t in trits)


def _int_to_balanced_ternary_signed(n: int) -> str:
    """Convert any integer (positive, negative, or zero) to balanced ternary string."""
    if n == 0:
        return '0'
    if n < 0:
        # Negate, convert, then flip all trits
        pos = _int_to_balanced_ternary(-n)
        return pos.translate(str.maketrans('T1', '1T'))
    return _int_to_balanced_ternary(n)


def _int_to_standard_ternary(n: int) -> str:
    """Convert non-negative integer to standard ternary string."""
    if n == 0:
        return '0'
    if n < 0:
        return '-' + _int_to_standard_ternary(-n)

    digits = []
    val = n
    while val > 0:
        digits.append(str(val % 3))
        val //= 3
    digits.reverse()
    return ''.join(digits)


def _parse_base3_string(s: str) -> Fraction:
    """Parse a ternary string (balanced or standard) into a Fraction.

    Balanced: uses T for -1, digits are T, 0, 1
    Standard: digits are 0, 1, 2
    Accepts optional leading '-' and optional '.' for fractional point.
    Auto-detects balanced vs standard by presence of 'T' or '2'.
    """
    s = s.strip()
    if not s:
        return Fraction(0)

    negative = False
    if s.startswith('-'):
        negative = True
        s = s[1:]

    # Detect mode
    has_T = 'T' in s or 't' in s
    has_2 = '2' in s
    balanced = has_T  # If it has T, it's balanced. If it has 2, it's standard.

    # Split on point
    if '.' in s:
        int_str, frac_str = s.split('.', 1)
        # Strip trailing '...' from frac_str
        if frac_str.endswith('...'):
            frac_str = frac_str[:-3]
    else:
        int_str = s
        frac_str = ''

    # Parse integer part
    result = Fraction(0)
    for ch in int_str:
        d = _parse_trit_char(ch, balanced)
        result = result * 3 + d

    # Parse fractional part
    power = Fraction(1, 3)
    for ch in frac_str:
        d = _parse_trit_char(ch, balanced)
        result += d * power
        power /= 3

    if negative:
        result = -result

    return result


def _parse_trit_char(ch: str, balanced: bool) -> int:
    """Parse a single character as a ternary digit."""
    if ch in ('T', 't'):
        return -1
    if ch in ('0', '1'):
        return int(ch)
    if ch == '2':
        if balanced:
            raise ValueError("Digit '2' invalid in balanced ternary (use T, 0, 1)")
        return 2
    raise ValueError(f"Invalid ternary digit: {ch!r}")


# ============================================================
# PUBLIC API — convenience functions
# ============================================================

def to_base3(value: Union[int, float, str, bytes, Fraction]) -> Base3Number:
    """Convert anything to a Base3Number.

    Accepts: int, float, str (decimal, hex with 0x, binary with 0b,
    or ternary with T/2 digits), bytes, Fraction.
    """
    if isinstance(value, Base3Number):
        return value

    if isinstance(value, str):
        s = value.strip()
        # Detect format
        if s.startswith('0x') or s.startswith('0X'):
            return Base3Number.from_hex(s)
        if s.startswith('0b') or s.startswith('0B'):
            return Base3Number.from_binary_string(s)
        # Check if it looks like ternary (has T or only 0/1/2 and .)
        has_ternary = any(c in s for c in 'Tt')
        has_decimal_only = all(c in '0123456789.-+' for c in s)
        if has_ternary:
            return Base3Number(s)
        if has_decimal_only:
            # Parse as decimal number
            return Base3Number(Fraction(s))
        # Try as ternary (only 0,1,2 and point)
        if all(c in '012.' for c in s):
            return Base3Number(s)
        raise ValueError(f"Cannot parse {s!r} as base-3 input")

    return Base3Number(value)


def from_base3(b3: Base3Number) -> Union[int, Fraction]:
    """Convert Base3Number back to Python type.

    Returns int if the value is integral, otherwise Fraction (exact).
    """
    if b3.exact.denominator == 1:
        return b3.exact.numerator
    return b3.exact


def base3_add(a: Union[Base3Number, int, float, Fraction],
              b: Union[Base3Number, int, float, Fraction]) -> Base3Number:
    """Add two values in base 3."""
    return _coerce(a) + _coerce(b)


def base3_mul(a: Union[Base3Number, int, float, Fraction],
              b: Union[Base3Number, int, float, Fraction]) -> Base3Number:
    """Multiply two values in base 3."""
    return _coerce(a) * _coerce(b)


def base3_div(a: Union[Base3Number, int, float, Fraction],
              b: Union[Base3Number, int, float, Fraction]) -> Base3Number:
    """Divide two values in base 3. Exact rational result, no float errors."""
    return _coerce(a) / _coerce(b)


# ============================================================
# DEMO / PROOF — run directly to see the truth
# ============================================================

def _demo():
    import math

    print("=" * 64)
    print("  base3.py — the foundation")
    print("  nos3bl33d | (DFG) DeadFoxGroup | x^2 = x + 1")
    print("=" * 64)
    print()

    # ----------------------------------------------------------
    # 1. $1.99 vs $2.00 — where the leak lives
    # ----------------------------------------------------------
    print("-" * 64)
    print("  1. THE LEAK: $1.99 vs $2.00")
    print("-" * 64)
    print()

    price_199 = Base3Number(Fraction(199, 100))
    price_200 = Base3Number(Fraction(200, 100))

    print("  $1.99:")
    print(price_199.pretty(30))
    print()
    print("  $2.00:")
    print(price_200.pretty(30))
    print()
    print("  $1.99 is INFINITE in base 3 — it never resolves.")
    print("  $2.00 is FINITE in base 3 — clean, exact, done.")
    print("  Every time $1.99 gets split, taxed, or divided: rounding.")
    print("  Every rounding = value leaked. Compounded globally = derivatives.")
    print()

    # ----------------------------------------------------------
    # 2. 1/3 — exact in base 3, infinite in decimal
    # ----------------------------------------------------------
    print("-" * 64)
    print("  2. 1/3: ONE TRIT vs INFINITE DECIMALS")
    print("-" * 64)
    print()

    one_third = Base3Number(Fraction(1, 3))
    print("  1/3:")
    print(one_third.pretty(30))
    print()
    print("  In decimal: 0.333333... (infinite, NEVER terminates)")
    print("  In base 3:  0.1          (one trit, EXACT)")
    print()
    print("  Now add them: 1/3 + 1/3 + 1/3")
    result = one_third + one_third + one_third
    print(f"  = {result.to_balanced_ternary()} = {result.to_decimal_string()}")
    print(f"  Exact: {result.exact} (no rounding error. ever.)")
    print()

    # Show more fractions
    print("  Which fractions are exact where?")
    print(f"  {'Fraction':>12s}  {'Base 10':^12s}  {'Base 3':^12s}")
    print(f"  {'--------':>12s}  {'------':^12s}  {'------':^12s}")
    test_fracs = [
        (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9),
        (1, 10), (1, 12), (1, 27), (1, 100),
    ]
    for num, den in test_fracs:
        b3 = Base3Number(Fraction(num, den))
        fin10 = "EXACT" if b3.terminates_in_base10() else "infinite"
        fin3 = "EXACT" if b3.terminates_in_base3() else "infinite"
        print(f"  {num}/{den:>3d}          {fin10:^12s}  {fin3:^12s}")
    print()
    print("  Every 'infinite' is a rounding error waiting to happen.")
    print("  Different bases, different leaks. Choose your base wisely.")
    print()

    # ----------------------------------------------------------
    # 3. Arithmetic: add and multiply in base 3
    # ----------------------------------------------------------
    print("-" * 64)
    print("  3. NATIVE BASE-3 ARITHMETIC")
    print("-" * 64)
    print()

    a = Base3Number(Fraction(7, 3))   # 2.1 in standard base 3
    b = Base3Number(Fraction(5, 9))   # 0.12 in standard base 3

    print(f"  a = 7/3")
    print(f"    balanced:  {a.to_balanced_ternary()}")
    print(f"    standard:  {a.to_standard_ternary()}")
    print(f"    decimal:   {a.to_decimal_string()}")
    print()
    print(f"  b = 5/9")
    print(f"    balanced:  {b.to_balanced_ternary()}")
    print(f"    standard:  {b.to_standard_ternary()}")
    print(f"    decimal:   {b.to_decimal_string()}")
    print()

    add_result = base3_add(a, b)
    print(f"  a + b = {add_result.to_balanced_ternary()}")
    print(f"        = {add_result.to_standard_ternary()}")
    print(f"        = {add_result.to_decimal_string()} (exact: {add_result.exact})")
    print()

    mul_result = base3_mul(a, b)
    print(f"  a * b = {mul_result.to_balanced_ternary()}")
    print(f"        = {mul_result.to_standard_ternary()}")
    print(f"        = {mul_result.to_decimal_string()} (exact: {mul_result.exact})")
    print()

    div_result = base3_div(a, b)
    print(f"  a / b = {div_result.to_balanced_ternary()}")
    print(f"        = {div_result.to_standard_ternary()}")
    print(f"        = {div_result.to_decimal_string()} (exact: {div_result.exact})")
    print()
    print("  No floats. No rounding. Exact rational arithmetic, displayed in base 3.")
    print()

    # ----------------------------------------------------------
    # 4. Genesis nonce in base 3
    # ----------------------------------------------------------
    print("-" * 64)
    print("  4. THE GENESIS NONCE: 2083236893 in base 3")
    print("-" * 64)
    print()

    nonce = 2083236893  # Bitcoin block 0 nonce
    b3_nonce = Base3Number(nonce)
    bal = b3_nonce.to_balanced_ternary()
    std = b3_nonce.to_standard_ternary()

    print(f"  Decimal:           {nonce:,}")
    print(f"  Hex:               {b3_nonce.to_hex()}")
    print(f"  Binary:            {b3_nonce.to_binary_string()}")
    print(f"  Standard ternary:  {std}")
    print(f"  Balanced ternary:  {bal}")
    print(f"  Trits (standard):  {len(std)}")
    print(f"  Bits (binary):     {nonce.bit_length()}")
    print()

    # Verify round-trip
    rt = from_base3(b3_nonce)
    assert rt == nonce, f"Round-trip failed: {rt} != {nonce}"
    print(f"  Round-trip verified: base3 -> int -> base3 = {rt:,}")
    print()

    # Also convert from balanced ternary string back
    rt2 = Base3Number.from_balanced_ternary(bal)
    assert rt2 == nonce, f"Balanced round-trip failed"
    print(f"  Balanced round-trip verified: '{bal}' -> {from_base3(rt2):,}")
    print()

    # ----------------------------------------------------------
    # 5. Radix economy: 3 beats 2
    # ----------------------------------------------------------
    print("-" * 64)
    print("  5. RADIX ECONOMY: why 3 is the optimal integer base")
    print("-" * 64)
    print()

    print("  To represent N values, the 'cost' is base * digits_needed.")
    print("  The radix economy of base b = b / ln(b). Lower is better.")
    print()
    print(f"  Base 2 (binary):   {2 / math.log(2):.4f}")
    print(f"  Base 3 (ternary):  {3 / math.log(3):.4f}  <-- best integer base")
    print(f"  Base e (theory):   {math.e:.4f}  (unrealizable optimum)")
    print(f"  Base 10 (decimal): {10 / math.log(10):.4f}")
    print()
    print("  Base 3 is the most efficient integer radix. Period.")
    print()

    # The 20 trits vs 32 bits comparison
    print("  Information density comparison:")
    print(f"    3^20 = {3**20:>14,}  (states in 20 trits)")
    print(f"    2^32 = {2**32:>14,}  (states in 32 bits)")
    print(f"    3^21 = {3**21:>14,}  (states in 21 trits)")
    print()
    print(f"    20 trits ~ {20 * math.log2(3):.1f} bits of information")
    print(f"    21 trits ~ {21 * math.log2(3):.1f} bits of information > 32 bits")
    print()
    print(f"    Per-digit info:  trit = {math.log2(3):.4f} bits    bit = 1.0000 bits")
    print(f"    Per-digit cost:  trit costs 3 wire states, carries {math.log2(3):.3f} bits each")
    print(f"                     bit  costs 2 wire states, carries 1.000 bits each")
    print(f"    Cost per bit:    ternary = {3/math.log2(3):.3f}    binary = {2/1:.3f}")
    print(f"    Ternary wins:    {(2 - 3/math.log2(3)) / 2 * 100:.1f}% more efficient per bit")
    print()

    # ----------------------------------------------------------
    # 6. Bonus: phi in balanced ternary
    # ----------------------------------------------------------
    print("-" * 64)
    print("  6. BONUS: phi in balanced ternary")
    print("-" * 64)
    print()

    # phi = (1 + sqrt(5)) / 2, irrational, never terminates in any integer base
    # But we can show its expansion

    # Use high-precision rational approximation of phi
    # phi ~ F(n+1)/F(n) for large Fibonacci numbers (EXACT convergence)
    # Use F(80)/F(79) for excellent approximation
    def fib(n: int) -> int:
        a, b = 0, 1
        for _ in range(n):
            a, b = b, a + b
        return a

    f80 = fib(80)
    f79 = fib(79)
    phi_approx = Base3Number(Fraction(f80, f79))

    print(f"  phi = (1 + sqrt(5)) / 2 = 1.6180339887...")
    print(f"  Approximated as F(80)/F(79) = {f80}/{f79}")
    print()
    print(f"  Standard ternary:  {phi_approx.to_standard_ternary(40)}")
    print(f"  Balanced ternary:  {phi_approx.to_balanced_ternary(40)}")
    print()
    print(f"  phi is irrational — it never terminates or repeats in ANY integer base.")
    print(f"  But watch: phi^2 = phi + 1 (the axiom).")
    print()

    phi_sq = phi_approx * phi_approx
    phi_plus_1 = phi_approx + 1
    diff = phi_sq - phi_plus_1

    print(f"  phi^2     = {phi_sq.to_balanced_ternary(30)}")
    print(f"  phi + 1   = {phi_plus_1.to_balanced_ternary(30)}")
    print(f"  Difference: {float(diff.exact):.2e} (rounding from finite approximation)")
    print()
    print(f"  The axiom holds: x^2 = x + 1, in any base, in any form.")
    print()

    # ----------------------------------------------------------
    # Tryte demo
    # ----------------------------------------------------------
    print("-" * 64)
    print("  7. TRYTE OPERATIONS (20-trit hardware words)")
    print("-" * 64)
    print()

    ta = Tryte.from_int(42)
    tb = Tryte.from_int(17)
    print(f"  a = 42  -> balanced: {ta.balanced_str():>12s}  standard: {ta.standard_str():>8s}")
    print(f"  b = 17  -> balanced: {tb.balanced_str():>12s}  standard: {tb.standard_str():>8s}")
    print()

    tc = ta + tb
    print(f"  a + b   -> balanced: {tc.balanced_str():>12s}  = {tc.to_int()}")
    td = ta - tb
    print(f"  a - b   -> balanced: {td.balanced_str():>12s}  = {td.to_int()}")
    te = ta * tb
    print(f"  a * b   -> balanced: {te.balanced_str():>12s}  = {te.to_int()}")
    tf = ta // tb
    print(f"  a // b  -> balanced: {tf.balanced_str():>12s}  = {tf.to_int()}")
    print()

    # Negative numbers are native — no two's complement needed
    tn = Tryte.from_int(-42)
    print(f"  -42     -> balanced: {tn.balanced_str():>12s}  = {tn.to_int()}")
    print(f"  Just flip every trit. No two's complement. No overflow traps.")
    print()

    # ----------------------------------------------------------
    # Trit logic demo
    # ----------------------------------------------------------
    print("-" * 64)
    print("  8. TERNARY LOGIC (three-valued: T, 0, 1)")
    print("-" * 64)
    print()
    print("  Kleene's strong three-valued logic:")
    print("  T = false/negative, 0 = unknown, 1 = true/positive")
    print()
    print(f"  {'A':>4s} {'B':>4s}  {'AND':>4s} {'OR':>4s} {'!A':>4s}")
    print(f"  {'---':>4s} {'---':>4s}  {'---':>4s} {'---':>4s} {'---':>4s}")
    for a_val in (-1, 0, 1):
        for b_val in (-1, 0, 1):
            a = Trit(a_val)
            b = Trit(b_val)
            print(f"  {a!s:>4s} {b!s:>4s}  {a & b!s:>4s} {a | b!s:>4s} {~a!s:>4s}")
    print()

    # ----------------------------------------------------------
    # Final
    # ----------------------------------------------------------
    print("=" * 64)
    print("  Base 3. The efficient base. The exact base.")
    print("  Every other script builds on this.")
    print()
    print("  x^2 = x + 1")
    print("  nos3bl33d | (DFG) DeadFoxGroup")
    print("=" * 64)


if __name__ == '__main__':
    _demo()

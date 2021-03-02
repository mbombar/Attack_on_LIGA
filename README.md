# Attack_on_LIGA
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.txt)
<img src="https://www.sagemath.org/pix/logo_sagemath+icon_oldstyle.png" width="100">
<img src="https://www.python.org/static/community_logos/python-logo-generic.svg" width="100">


A message recovery attack on LIGA, a code based cryptosystem

---

This repository provides a SageMath implementation of:

 * A static rank error channel, in the same fashion as SageMath ```StaticErrorChannel``` for the Hamming Metric.
 * A right-hand side analogue of the Welch-Berlekamp decoding algorithm for Gabidulin Codes.
 * Some useful functions to handle skew polynomial reduction.
 * A message-recovery attack on the [LIGA](https://arxiv.org/abs/1812.04892) cryptosystem.

 <!-- This is a research code, paired with the article -->

 <!-- **Decoding supercodes of Gabidulin codes and applications to cryptanalysis** -->
 <!-- *Maxime Bombar and Alain Couvreur* -->

The implementation of the LIGA cryptosystem is based on the original implementation of Julian Renner, Sven Puchinger and Antonia Wachter-Zeh, available [here](https://bitbucket.org/julianrenner/liga_pke/src).

The implementation of Gabidulin codes is based on the current implementation in SageMath introduced in SageMath version 9.1.


<!-- #### How to cite: -->

<!-- ``` bibtex -->
<!-- @misc{BC21 -->
<!--     author = {Maxime Bombar and Alain Couvreur}, -->
<!--     title = {Decoding supercodes of Gabidulin codes and applications to cryptanalysis}, -->
<!--     year = {2021}, -->
<!-- } -->
<!-- ``` -->


## Requirements

 * SageMath version >= 9.2
 * Python version >= 3.9

## How to use it

### Right-hand side Berlekamp-Welch algorithm

The provided implementation can be used directly as usual Gabidulin codes in SageMath (from version 9.1).

``` python
from gabidulin_code import GabidulinCode
from rank_channel import StaticRankErrorChannel

p = 2
s = 1
q = p^s
m = 41
n = m
k = 27

Fqm = GF(q^m)
Fq = GF(q)
twist = Fqm.frobenius_endomorphism(s)

G = GabidulinCode(Fqm, n, k, Fq, twist)
t = (G.minimum_distance() - 1)//2
Chan = StaticRankErrorChannel(G.ambient_space(), t, G.sub_field())

c = G.random_element()
y = Chan(c)

d = G.decode_to_code(y, decoder_name="RightBerlekampWelch")

assert d == c
```

### Message recovery attack on LIGA

``` sh
sage attack_on_rfl.sage
```

```sh

# ----------- Repaired Version (LIGA, toy example, zeta=2) -----------------

m=10
k=4
w=4
u=3

# ----------- Key Generation -------------------
Key Generation done

# ------------------ Encryption --------------------------
Encryption done

# ------------------- Decryption ---------------------------------
Encryted message equals the decrypted message ? --> True

# ------------------- Attack ---------------------------------

# ------------------- Step 1: Get rid of the small error ---------------------------------

# ------------------- Step 2: Remove the z dependency-----------------------------

# ------------------- Step 3: Recover the plaintext-----------------------------
I Won !!
Attack lasted 0.67 seconds !
```

The parameters can be changed in the [source code](https://github.com/mbombar/Attack_on_LIGA/blob/main/attack_on_RFL.sage#L24).

### Results

| Name     | Claimed security level | Average running time |
|----------|------------------------|----------------------|
| LIGA-128 | 128 bits               | 8 minutes            |
| LIGA-192 | 192 bits               | 27 minutes           |
| LIGA-256 | 256 bits               | 92 minutes           |

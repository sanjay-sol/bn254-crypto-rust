# bn254-crypto-rust

This repository contains a **from-first-principles, educational implementation** of the BN254 elliptic curve and its finite-field arithmetic.  
The goal is to build the entire curve **step-by-step**, the same way real zkSNARK libraries (arkworks, blst, gnark, halo2) do internally.

This project is **not meant for production use** — it is a learning journey to understand **how pairing-friendly curves work** under the hood.

---

# ✅ What’s Implemented So Far

## 1. **Base Field 𝔽ₚ (`fp.rs`)**
We implemented the prime field **Fp = ℤ/pℤ** where:

`p = 21888242871839275222246405745257275088548364400416034343698204186575808495617`


This is the BN254 scalar field and base field used in Ethereum zkSNARKs.

### Features implemented:
- Modular addition  
- Modular subtraction  
- Modular multiplication  
- Modular negation  
- Modular inverse (Extended Euclidean Algorithm)  
- Modular exponentiation  
- Property tests + random tests  
- Operator overloading (`+`, `-`, `*`, `-x`)

The field ensures:
- All results are automatically reduced mod p  
- No intermediate overflow leaks outside of the field  

---

## 2. **G1 Group (`g1.rs`)**
This is the elliptic curve group on BN254 defined over 𝔽ₚ:

`E: y² = x³ + 3`


### Coordinate system:
We use **Jacobian coordinates (X:Y:Z)**:

`x = X / Z²`
`y = Y / Z³`

### Why Jacobian?
- No field inversions during addition or doubling  
- Inversions are expensive (~100× slower than multiplications)  
- Jacobian avoids them → faster operations

### Implemented operations:
✅ Point at infinity  
✅ Affine conversion  
✅ Curve validity check  
✅ Point doubling  
✅ Point addition (complete, handles all edge cases)  
✅ Scalar multiplication (double-and-add)

### Tests:
- Infinity behavior  
- Affine conversion  
- Curve equation validity  
- Double = P + P  
- Add is commutative  
- Scalar multiplication 

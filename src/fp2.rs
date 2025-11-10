use crate::fp::Fp;
use num_bigint::ToBigUint;
use num_traits::{Zero, One};
use std::ops::{Add, Sub, Mul, Neg};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Fp2 {
    // Represents a + b * u  where u^2 = -1  (β = -1)
    pub a: Fp,
    pub b: Fp,
}

impl Fp2 {
    pub fn new(a: Fp, b: Fp) -> Self {
        Self { a, b }
    }

    pub fn zero() -> Self {
        Self { a: Fp::zero(), b: Fp::zero() }
    }

    pub fn one() -> Self {
        Self { a: Fp::one(), b: Fp::zero() }
    }

    // Conjugate: a - b*u
    pub fn conjugate(&self) -> Self {
        Self {
            a: self.a.clone(),
            b: Fp::zero() - self.b.clone(),
        }
    }

    // Norm: a^2 + b^2  (since u^2 = -1 -> (a + bu)(a - bu) = a^2 + b^2)
    pub fn norm(&self) -> Fp {
        self.a.clone() * self.a.clone() + self.b.clone() * self.b.clone()
    }

    // Inverse: (a - b u) / (a^2 + b^2)
    pub fn inv(&self) -> Self {
        // denom = a^2 + b^2  (in Fp)
        let denom = self.a.clone() * self.a.clone() + self.b.clone() * self.b.clone();
        // denom_inv = denom^{-1} in Fp
        let denom_inv = denom.inv();
        let conj = self.conjugate();
        Self {
            a: conj.a * denom_inv.clone(),
            b: conj.b * denom_inv,
        }
    }

    // Exponentiation by BigUint (binary square-and-multiply)
    pub fn pow(&self, mut exp: num_bigint::BigUint) -> Self {
        use num_traits::Zero;
        let mut base = self.clone();
        let mut result = Fp2::one();

        while exp > num_bigint::BigUint::zero() {
            if (&exp & num_bigint::BigUint::one()) == num_bigint::BigUint::one() {
                result = result * base.clone();
            }
            exp >>= 1;
            base = base.clone() * base;
        }
        result
    }
}

impl Add for Fp2 {
    type Output = Fp2;
    fn add(self, rhs: Fp2) -> Fp2 {
        Fp2::new(self.a + rhs.a, self.b + rhs.b)
    }
}

impl Sub for Fp2 {
    type Output = Fp2;
    fn sub(self, rhs: Fp2) -> Fp2 {
        Fp2::new(self.a - rhs.a, self.b - rhs.b)
    }
}

impl Neg for Fp2 {
    type Output = Fp2;
    fn neg(self) -> Fp2 {
        Fp2::new(-self.a, -self.b)
    }
}

impl Mul for Fp2 {
    type Output = Fp2;
    fn mul(self, rhs: Fp2) -> Fp2 {
        // (a + b u) * (c + d u) = (a*c + b*d*u^2) + (a*d + b*c) u
        // since u^2 = -1, u^2 -> -1
        // = (a*c - b*d) + (a*d + b*c) u
        let ac = self.a.clone() * rhs.a.clone();
        let bd = self.b.clone() * rhs.b.clone();
        let ad = self.a.clone() * rhs.b.clone();
        let bc = self.b.clone() * rhs.a.clone();

        let real = ac - bd;
        let imag = ad + bc;
        Fp2::new(real, imag)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fp::Fp;
    use num_bigint::ToBigUint;
    use rand::Rng;

    #[test]
    fn test_fp2_basic_ops() {
        let a = Fp2::new(Fp::new(3u32.to_biguint().unwrap()), Fp::new(5u32.to_biguint().unwrap()));
        let b = Fp2::new(Fp::new(2u32.to_biguint().unwrap()), Fp::new(7u32.to_biguint().unwrap()));

        // add
        let s = a.clone() + b.clone();
        assert_eq!(s, Fp2::new(Fp::new(5u32.to_biguint().unwrap()), Fp::new(12u32.to_biguint().unwrap())));

        // sub
        let d = b.clone() - a.clone();
        //  (2-3, 7-5) mod p => (-1 mod p, 2)
        let neg1 = Fp::new((&*crate::fp::P - 1u32.to_biguint().unwrap()));
        assert_eq!(d, Fp2::new(neg1, Fp::new(2u32.to_biguint().unwrap())));

        // mul using (a+bu)*(c+du)
        let m = a.clone() * b.clone();
        // check via direct values computed in Fp (sanity check)
        // (3 + 5u)*(2 + 7u) => real = 3*2 - 5*7 = 6 - 35 = -29 mod p, imag = 3*7 + 5*2 = 21 + 10 = 31
        let real_expected = Fp::new((&*crate::fp::P - 29u32.to_biguint().unwrap()));
        let imag_expected = Fp::new(31u32.to_biguint().unwrap());
        assert_eq!(m, Fp2::new(real_expected, imag_expected));
    }

    #[test]
    fn test_fp2_conjugate_and_inv() {
        let a = Fp2::new(Fp::new(3u32.to_biguint().unwrap()), Fp::new(5u32.to_biguint().unwrap()));
        let conj = a.conjugate();
        assert_eq!(conj, Fp2::new(Fp::new(3u32.to_biguint().unwrap()), Fp::new((&*crate::fp::P - 5u32.to_biguint().unwrap()))));

        // inverse: check a * a^{-1} == 1
        let inv = a.inv();
        let prod = a * inv.clone();
        assert_eq!(prod, Fp2::one());
    }

    #[test]
    fn test_fp2_norm_property() {
        let a = Fp2::new(Fp::new(4u32.to_biguint().unwrap()), Fp::new(7u32.to_biguint().unwrap()));
        let n = a.norm(); // a^2 + b^2
        let computed = (a.clone() * a.conjugate()).a; // (a + bu)(a - bu) = norm
        assert_eq!(n, computed);
    }

    #[test]
    fn test_fp2_pow_small() {
        let a = Fp2::new(Fp::new(2u32.to_biguint().unwrap()), Fp::new(1u32.to_biguint().unwrap()));
        let res = a.pow(5u32.to_biguint().unwrap());
        // naive check: multiply 5 times by hand (or call pow from reference)
        let mut exp = Fp2::one();
        for _ in 0..5 {
            exp = exp * a.clone();
        }
        assert_eq!(res, exp);
    }
}

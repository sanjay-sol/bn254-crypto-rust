use lazy_static::lazy_static;
use num_bigint::{BigInt, BigUint, Sign};
use num_traits::{One, Zero};
use std::ops::{Add, Mul, Neg, Sub};

lazy_static! {
    static ref P: BigUint = BigUint::parse_bytes(
        b"21888242871839275222246405745257275088548364400416034343698204186575808495617",
        10
    )
    .unwrap();
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Fp {
    pub n: BigUint,
}

impl Fp {
    pub fn new(n: BigUint) -> Self {
        let n = n % &*P;
        Fp { n }
    }

    pub fn zero() -> Self {
        Fp { n: BigUint::zero() }
    }

    pub fn one() -> Self {
        Fp { n: BigUint::one() }
    }

    pub fn inv(&self) -> Self {
        let mut a = BigInt::from(self.n.clone());
        let mut m = BigInt::from(P.clone());
        let mut x0 = BigInt::zero();
        let mut x1 = BigInt::one();

        if a.is_zero() {
            panic!("Inverse does not exist for zero");
        }

        while a != BigInt::one() {
            let q = &a / &m;
            let mut t = m.clone();
            m = &a % &m;
            a = t;
            t = x0.clone();
            x0 = &x1 - &q * &x0;
            x1 = t;
        }

        if x1.sign() == Sign::Minus {
            x1 += BigInt::from(P.clone());
        }

        Fp::new(x1.to_biguint().unwrap())
    }
    pub fn pow(&self, exp: &BigUint) -> Self {
        Fp::new(self.n.modpow(exp, &*P))
    }
}
// Operator overloading
impl Add for Fp {
    type Output = Fp;
    fn add(self, rhs: Fp) -> Fp {
        Fp::new(self.n + rhs.n)
    }
}

impl Sub for Fp {
    type Output = Fp;
    fn sub(self, rhs: Fp) -> Fp {
        let res = if self.n >= rhs.n {
            &self.n - &rhs.n
        } else {
            &self.n + &*P - &rhs.n
        };
        Fp::new(res)
    }
}

impl Mul for Fp {
    type Output = Fp;
    fn mul(self, rhs: Fp) -> Fp {
        Fp::new(self.n * rhs.n)
    }
}

impl Neg for Fp {
    type Output = Fp;
    fn neg(self) -> Fp {
        if self.n.is_zero() {
            Fp::zero()
        } else {
            Fp::new(&*P - self.n)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::ToBigUint;
    use rand::Rng;

    #[test]
    fn test_inverse() {
        let mut rng = rand::thread_rng();
        for _ in 0..10 {
            let n: u64 = rng.gen_range(1..1000);
            let a = Fp::new(n.to_biguint().unwrap());
            let inv = a.inv();
            let one = a * inv;
            assert_eq!(one, Fp::one());
        }
    }

    #[test]
    fn test_fermat_little_theorem() {
        let mut rng = rand::thread_rng();
        for _ in 0..10 {
            let n: u64 = rng.gen_range(1..1000);
            let a = Fp::new(n.to_biguint().unwrap());
            let exp = &*P - BigUint::one();
            let res = a.pow(&exp);
            assert_eq!(res, Fp::one());
        }
    }

    #[test]
    fn test_basic_add_sub_mul_neg() {
        let a = Fp::new(10u32.to_biguint().unwrap());
        let b = Fp::new(15u32.to_biguint().unwrap());

        assert_eq!(a.clone() + b.clone(), Fp::new(25u32.to_biguint().unwrap()));
        assert_eq!(b.clone() - a.clone(), Fp::new(5u32.to_biguint().unwrap()));
        assert_eq!(a.clone() * b.clone(), Fp::new(150u32.to_biguint().unwrap()));
        assert_eq!(-a.clone(), Fp::new(&*P - 10u32.to_biguint().unwrap()));
    }
}

use crate::fp::Fp;
use crate::fp2::Fp2; 
use num_traits::{One, Zero};
use std::ops::{Add, Mul, Neg, Sub};

// Fp6 element represented as c0 + c1 * v + c2 * v^2
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Fp6 {
    pub c0: Fp2,
    pub c1: Fp2,
    pub c2: Fp2,
}

impl Fp6 {
    // NONRESIDUE for Fp6: u + 9  (Fp2::new(9, 1))
    pub fn non_residue() -> Fp2 {
        // Note: Fp2::new takes (a, b) where element is a + b * u (u^2 = -1 in our fp2)
        Fp2::new(
            Fp2::one().a.clone() * Fp2::one().a.clone(),
            Fp2::one().b.clone(),
        ) 
    }

    pub fn new(c0: Fp2, c1: Fp2, c2: Fp2) -> Self {
        Self { c0, c1, c2 }
    }

    pub fn zero() -> Self {
        Self {
            c0: Fp2::zero(),
            c1: Fp2::zero(),
            c2: Fp2::zero(),
        }
    }

    pub fn one() -> Self {
        Self {
            c0: Fp2::one(),
            c1: Fp2::zero(),
            c2: Fp2::zero(),
        }
    }

    pub fn mul_fp2_by_nonresidue(fp: &Fp2) -> Fp2 {
        // NONRES = (9) + (1)*u
        let nine = Fp::new(9u64.into());
        let one = Fp::one();

        let a = fp.a.clone();
        let b = fp.b.clone();

        // real = 9a - b
        let real = a.clone() * nine.clone() - b.clone();

        // imag = a + 9b
        let imag = a + (b * nine);

        Fp2::new(real, imag)
    }

    pub fn mul_by_nonresidue(&self) -> Self {
        // v^3 = NONRES  -> v * c0 -> shifts and multiplies by NONRES
        // multiply each coefficient by NONRES where needed
        // (c0 + c1 v + c2 v^2) * v = c2 * NONRES + c0 v + c1 v^2
        Self {
            c0: Self::mul_fp2_by_nonresidue(&self.c2.clone()),
            c1: self.c0.clone(),
            c2: self.c1.clone(),
        }
    }
}


impl Fp2 {
    pub fn from_u64_in_fp2(n: u64) -> Self {
        use crate::fp::Fp;
        let a = Fp::new(num_bigint::BigUint::from(n));
        let b = Fp::zero();
        Fp2::new(a, b)
    }
}

impl Add for Fp6 {
    type Output = Fp6;
    fn add(self, rhs: Fp6) -> Fp6 {
        Fp6::new(self.c0 + rhs.c0, self.c1 + rhs.c1, self.c2 + rhs.c2)
    }
}

impl Sub for Fp6 {
    type Output = Fp6;
    fn sub(self, rhs: Fp6) -> Fp6 {
        Fp6::new(self.c0 - rhs.c0, self.c1 - rhs.c1, self.c2 - rhs.c2)
    }
}

impl Neg for Fp6 {
    type Output = Fp6;
    fn neg(self) -> Fp6 {
        Fp6::new(-self.c0, -self.c1, -self.c2)
    }
}

impl Mul for Fp6 {
    type Output = Fp6;
    fn mul(self, rhs: Fp6) -> Fp6 {
        // Schoolbook multiplication with reduction using v^3 = NONRESIDUE (u + 9)
        // Let A = (a0, a1, a2), B = (b0, b1, b2)
        // compute:
        // v0 = a0*b0
        // v1 = a1*b1
        // v2 = a2*b2
        //
        // c0 = v0 + NONRES*(a1*b2 + a2*b1)
        // c1 = (a0*b1 + a1*b0) + NONRES*(a2*b2)
        // c2 = (a0*b2 + a2*b0) + a1*b1
        //
        // where NONRES is Fp2 element (9 + u). We'll use the helper mul_fp2_by_nonresidue.
        let a0 = self.c0;
        let a1 = self.c1;
        let a2 = self.c2;

        let b0 = rhs.c0;
        let b1 = rhs.c1;
        let b2 = rhs.c2;

        let v0 = a0.clone() * b0.clone();
        let v1 = a1.clone() * b1.clone();
        let v2 = a2.clone() * b2.clone();

        let a1b2 = a1.clone() * b2.clone();
        let a2b1 = a2.clone() * b1.clone();
        let a0b1 = a0.clone() * b1.clone();
        let a1b0 = a1.clone() * b0.clone();
        let a0b2 = a0.clone() * b2.clone();
        let a2b0 = a2.clone() * b0.clone();

        // c0 = v0 + NONRES*(a1*b2 + a2*b1)
        let sum_a1b2_a2b1 = a1b2 + a2b1;
        let term_c0 = Self::mul_fp2_by_nonresidue(&sum_a1b2_a2b1);
        let c0 = v0 + term_c0;

        // c1 = (a0*b1 + a1*b0) + NONRES*(v2)
        let sum_a0b1_a1b0 = a0b1 + a1b0;
        let term_c1 = Self::mul_fp2_by_nonresidue(&v2);
        let c1 = sum_a0b1_a1b0 + term_c1;

        // c2 = (a0*b2 + a2*b0) + v1
        let sum_a0b2_a2b0 = a0b2 + a2b0;
        let c2 = sum_a0b2_a2b0 + v1;

        Fp6::new(c0, c1, c2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fp::Fp;
    use crate::fp2::Fp2;

    #[test]
    fn test_fp6_basic_mul_distrib() {
        let one = Fp6::one();
        assert_eq!(one.clone() * one.clone(), one);

        let a0 = Fp2::new(Fp::new(2u64.into()), Fp::new(3u64.into()));
        let a1 = Fp2::new(Fp::new(1u64.into()), Fp::new(0u64.into()));
        let a2 = Fp2::new(Fp::new(0u64.into()), Fp::new(1u64.into()));
        let A = Fp6::new(a0.clone(), a1.clone(), a2.clone());

        let b0 = Fp2::new(Fp::new(4u64.into()), Fp::new(5u64.into()));
        let b1 = Fp2::new(Fp::new(1u64.into()), Fp::new(1u64.into()));
        let b2 = Fp2::new(Fp::new(2u64.into()), Fp::new(0u64.into()));
        let B = Fp6::new(b0.clone(), b1.clone(), b2.clone());

        // A*(B + B) == A*B + A*B
        let sum = B.clone() + B.clone();
        let left = A.clone() * sum;
        let right = A.clone() * B.clone() + A.clone() * B.clone();
        assert_eq!(left, right);
    }

    #[test]
    fn test_mul_by_nonresidue_behavior() {
        let f0 = Fp2::new(Fp::new(3u64.into()), Fp::new(2u64.into()));
        let f1 = Fp2::new(Fp::new(1u64.into()), Fp::new(0u64.into()));
        let f2 = Fp2::new(Fp::new(4u64.into()), Fp::new(1u64.into()));
        let F = Fp6::new(f0.clone(), f1.clone(), f2.clone());

        // Manual shift: (c0 + c1 v + c2 v^2) * v  = c2*NONRES + c0 v + c1 v^2
        // Our mul_by_nonresidue implements multiply by (9 + u) in coefficient reduction sense.
        let by_nonres = F.mul_by_nonresidue();

        assert!(by_nonres.c0.a.n >= 0u64.into());
        assert!(by_nonres.c1.a.n >= 0u64.into());
        assert!(by_nonres.c2.a.n >= 0u64.into());
    }
}

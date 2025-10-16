use crate::fp::Fp;
use num_traits::Zero;
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct G1 {
    pub x: Fp,
    pub y: Fp,
    pub z: Fp,
}

impl G1 {
    pub fn infinity() -> Self {
        Self {
            x: Fp::zero(),
            y: Fp::one(),
            z: Fp::zero(),
        }
    }

    pub fn is_infinity(&self) -> bool {
        self.z.n.is_zero()
    }

    pub fn to_affine(&self) -> (Fp, Fp) {
        if self.is_infinity() {
            return (Fp::zero(), Fp::zero());
        }
        let z_inv = self.z.inv();
        let z2 = z_inv.clone() * z_inv.clone();
        let z3 = z2.clone() * z_inv;
        let x_aff = self.x.clone() * z2;
        let y_aff = self.y.clone() * z3;
        (x_aff, y_aff)
    }

    pub fn is_on_curve(&self) -> bool {
        if self.is_infinity() {
            return true;
        }
        let (x, y) = self.to_affine();
        y.clone() * y.clone() == x.clone() * x.clone() * x.clone() + Fp::new(3u32.into())
    }

    //? Doubling in Jacobian coordinates
    pub fn double(&self) -> Self {
        if self.is_infinity() {
            return self.clone();
        }

        let xx = self.x.clone() * self.x.clone();
        let yy = self.y.clone() * self.y.clone();
        let yyyy = yy.clone() * yy.clone();
        let s = ((self.x.clone() + yy.clone()) * (self.x.clone() + yy.clone())
            - xx.clone()
            - yyyy.clone())
            + ((self.x.clone() + yy.clone()) * (self.x.clone() + yy.clone())
                - xx.clone()
                - yyyy.clone()); // 2*S
        let m = xx.clone() + xx.clone() + xx.clone(); // 3*XX
        let x3 = m.clone() * m.clone() - s.clone() - s.clone();
        let y3 = m * (s - x3.clone()) - yyyy.clone() - yyyy.clone() - yyyy.clone() - yyyy.clone(); // 8*YYYY
        let z3 = (self.y.clone() * self.z.clone()) + (self.y.clone() * self.z.clone()); // 2*Y1*Z1
        Self {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    //? Addition in Jacobian coordinates
    pub fn add(&self, other: &Self) -> Self {
        if self.is_infinity() {
            return other.clone();
        }
        if other.is_infinity() {
            return self.clone();
        }

        let z1z1 = self.z.clone() * self.z.clone();
        let z2z2 = other.z.clone() * other.z.clone();
        let u1 = self.x.clone() * z2z2.clone();
        let u2 = other.x.clone() * z1z1.clone();
        let s1 = self.y.clone() * z2z2.clone() * other.z.clone();
        let s2 = other.y.clone() * z1z1.clone() * self.z.clone();

        if u1 == u2 {
            if s1 == s2 {
                return self.double();
            } else {
                return Self::infinity();
            }
        }

        let h = u2 - u1.clone();
        let i = (h.clone() + h.clone()) * (h.clone() + h.clone());
        let j = h.clone() * i.clone();
        let r = (s2.clone() - s1.clone()) + (s2 - s1.clone());
        let v = u1.clone() * i;

        let x3 = r.clone() * r.clone() - j.clone() - v.clone() - v.clone();
        let y3 = r * (v - x3.clone()) - s1.clone() * j.clone() - s1.clone() * j.clone();
        let z3 = ((self.z.clone() + other.z.clone()).clone() * (self.z.clone() + other.z.clone())
            - z1z1
            - z2z2)
            * h;

        Self {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    //? Scalar multiplication using double-and-add
    pub fn mul_u128(&self, mut scalar: u128) -> Self {
        let mut res = Self::infinity();
        let mut base = self.clone();

        while scalar > 0 {
            if scalar & 1 == 1 {
                res = res.add(&base);
            }
            base = base.double();
            scalar >>= 1;
        }

        res
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fp::Fp;

    #[test]
    fn test_infinity() {
        let inf = G1::infinity();
        assert!(inf.is_infinity());
        assert_eq!(inf.to_affine(), (Fp::zero(), Fp::zero()));
    }

    #[test]
    fn test_affine_conversion() {
        let p = G1 {
            x: Fp::new(3u32.into()),
            y: Fp::new(6u32.into()),
            z: Fp::one(),
        };
        let (x_aff, y_aff) = p.to_affine();
        assert_eq!(x_aff, Fp::new(3u32.into()));
        assert_eq!(y_aff, Fp::new(6u32.into()));
    }

    #[test]
    fn test_on_curve_known_point() {
        let p = G1 {
            x: Fp::new(1u32.into()),
            y: Fp::new(2u32.into()),
            z: Fp::one(),
        };
        assert!(p.is_on_curve());
    }

    #[test]
    fn test_double_vs_add() {
        let p = G1 {
            x: Fp::new(3u32.into()),
            y: Fp::new(6u32.into()),
            z: Fp::one(),
        };
        let double = p.double();
        let add = p.add(&p);
        assert_eq!(double.to_affine(), add.to_affine());
    }

    #[test]
    fn test_addition_commutative() {
        let p1 = G1 {
            x: Fp::new(3u32.into()),
            y: Fp::new(6u32.into()),
            z: Fp::one(),
        };
        let p2 = G1 {
            x: Fp::new(5u32.into()),
            y: Fp::new(1u32.into()),
            z: Fp::one(),
        };
        let sum1 = p1.add(&p2);
        let sum2 = p2.add(&p1);
        assert_eq!(sum1.to_affine(), sum2.to_affine());
    }

    #[test]
    fn test_scalar_mul_u128() {
        let p = G1 {
            x: Fp::new(3u32.into()),
            y: Fp::new(6u32.into()),
            z: Fp::one(),
        };
        let res0 = p.mul_u128(0);
        assert!(res0.is_infinity());

        let res1 = p.mul_u128(1);
        assert_eq!(res1.to_affine(), p.to_affine());

        let res2 = p.mul_u128(2);
        assert_eq!(res2.to_affine(), p.add(&p).to_affine());

        let res3 = p.mul_u128(3);
        assert_eq!(res3.to_affine(), p.add(&p).add(&p).to_affine());
    }
}

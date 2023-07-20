use std::str::FromStr;

use super::{gadget_parameters::GadgetParametersRing, gadget_ring::find_solution_gadget_ring};
use qfall_math::{
    integer::{MatPolyOverZ, PolyOverZ, Z},
    integer_mod_q::{MatPolynomialRingZq, PolyOverZq, PolynomialRingZq, Zq},
    traits::{Concatenate, GetEntry, GetNumColumns, GetNumRows, Pow, SetCoefficient, SetEntry},
};

pub fn gen_short_basis_for_trapdoor_ring(
    params: &GadgetParametersRing,
    a: &MatPolynomialRingZq,
    r: &MatPolyOverZ,
    e: &MatPolyOverZ,
) -> MatPolyOverZ {
    let sa_l = gen_sa_l(e, r);
    let sa_r = gen_sa_r(params, a);
    let k = PolyOverZq::from(&params.modulus).get_degree();
    let mut poly_degrees = MatPolyOverZ::new(1, k);
    for i in 0..k {
        let mut x_i = PolyOverZ::default();
        x_i.set_coeff(i, 1).unwrap();
        poly_degrees.set_entry(0, i, x_i).unwrap();
    }
    poly_degrees.tensor_product(sa_l * sa_r)
}

fn gen_sa_l(e: &MatPolyOverZ, r: &MatPolyOverZ) -> MatPolyOverZ {
    let out = e.concat_vertical(r).unwrap();

    let identity_lower_right = MatPolyOverZ::identity(out.get_num_columns(), out.get_num_columns());

    let out = out.concat_vertical(&identity_lower_right).unwrap();
    let identity_left = MatPolyOverZ::identity(out.get_num_rows(), 2);

    identity_left.concat_horizontal(&out).unwrap()
}

fn gen_sa_r(params: &GadgetParametersRing, a: &MatPolynomialRingZq) -> MatPolyOverZ {
    let ident = MatPolyOverZ::identity(2, 2);
    ident.concat_vertical(&compute_w(params, a)).unwrap()
}

/// gW = - a[I_2|0]
fn compute_w(params: &GadgetParametersRing, a: &MatPolynomialRingZq) -> MatPolyOverZ {
    let minus_one = PolynomialRingZq::from((&PolyOverZ::from(-1), &params.modulus));
    let rhs_1: PolynomialRingZq = a.get_entry(0, 0).unwrap();
    let rhs_2: PolynomialRingZq = a.get_entry(0, 1).unwrap();
    let w_1 =
        find_solution_gadget_ring(&(&minus_one * &rhs_1), &params.k, &params.base).transpose();
    let w_2 =
        find_solution_gadget_ring(&(&minus_one * &rhs_2), &params.k, &params.base).transpose();
    w_1.concat_horizontal(&w_2).unwrap()
}

fn compute_s(params: &GadgetParametersRing) -> MatPolyOverZ {
    let id_k = MatPolyOverZ::identity(&params.k, &params.k);
    let mut sk = &params.base * id_k;
    for i in 0..(sk.get_num_rows() - 1) {
        sk.set_entry(i + 1, i, PolyOverZ::from(-1)).unwrap();
    }
    sk = if params.base.pow(&params.k).unwrap() == Z::from(&params.q) {
        // compute s in the special case where the modulus is a power of base
        // i.e. the last column can remain as it is
        sk
    } else {
        // compute s for any arbitrary modulus
        // represent modulus in `base` and set last row accordingly
        let mut q = Z::from(&params.q);
        for i in 0..(sk.get_num_rows()) {
            let q_i = Zq::from((&q, &params.base)).get_value();
            sk.set_entry(i, sk.get_num_columns() - 1, PolyOverZ::from(q_i))
                .unwrap();
            q = q - q_i;
            q = q.div_exact(&params.base).unwrap();
        }
        sk
    };
    MatPolyOverZ::identity(&params.n, &params.n).tensor_product(&sk)
}

#[cfg(test)]
mod test_short_basis {
    use super::gen_short_basis_for_trapdoor_ring;
    use crate::sample::g_trapdoor::{
        gadget_parameters::GadgetParametersRing, gadget_ring::gen_trapdoor_ring_lwe,
    };
    use qfall_math::integer::{PolyOverZ, Z};
    use qfall_math::integer_mod_q::{MatPolynomialRingZq, Modulus, PolynomialRingZq};
    use qfall_math::rational::Q;
    use qfall_math::traits::GetEntry;

    #[test]
    fn working() {
        let params =
            GadgetParametersRing::init_default(3, &Modulus::try_from(&Z::from(16)).unwrap());
        let a_bar = PolyOverZ::sample_uniform(&params.n, 0, &params.q).unwrap();

        let (a, r, e) = gen_trapdoor_ring_lwe(&params, &a_bar, &Q::from(5)).unwrap();

        let short_base = gen_short_basis_for_trapdoor_ring(&params, &a, &r, &e);
        let short_base = MatPolynomialRingZq::from((&short_base, &params.modulus));

        let res = a * short_base;
        let res1: PolynomialRingZq = res.get_entry(0, 0).unwrap();
        let res2: PolynomialRingZq = res.get_entry(0, 0).unwrap();
        println!("{}", res1.get_poly());
        println!("{}", res2.get_poly());
    }
}

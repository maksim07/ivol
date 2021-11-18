#![cfg(test)]
use crate::black_scholes;
use crate::black_scholes::{BlackScholesParams, simulate_put};
use crate::black_scholes::simulate_call;

const BUMP: f64 = 0.00001;
const EPS: f64 = 0.0001;

#[test]
fn test_call_premium() {
    let mut cp = black_scholes::call_premium(&BlackScholesParams {
        price: 100f64,
        div_yield: 0f64,
        strike: 110f64,
        vol: 0.12,
        rate: 0.05,
        time_to_expiry: 1.5
    });
    assert!((cp - 4.95).abs() < 0.01);

    cp = black_scholes::call_premium(&BlackScholesParams {
        price: 100f64,
        div_yield: 0.01,
        strike: 110f64,
        vol: 0.12,
        rate: 0.05,
        time_to_expiry: 1.5
    });
    assert!((cp - 4.27).abs() < 0.01);
}

#[test]
fn test_put_premium() {
    let mut cp = black_scholes::put_premium(&BlackScholesParams {
        price: 150f64,
        div_yield: 0.1,
        strike: 200f64,
        vol: 0.11,
        rate: 0.01,
        time_to_expiry: 1.5
    });
    assert!((cp - 67.92).abs() < 0.01);

    cp = black_scholes::put_premium(&BlackScholesParams {
        price: 150f64,
        div_yield: 0.1,
        strike: 160f64,
        vol: 0.06,
        rate: 0.03,
        time_to_expiry: 0.5
    });
    assert!((cp - 14.95).abs() < 0.01);
}

#[test]
fn test_delta() {
    let bs_params = BlackScholesParams {
        price: 1045.0,
        div_yield: 0.06,
        strike: 1020.0,
        vol: 0.24,
        rate: 0.025,
        time_to_expiry: 0.5
    };

    let delta = black_scholes::call_delta(&bs_params);
    assert!((delta - 0.528).abs() < 0.01);

    let bump = |c: &BlackScholesParams| {BlackScholesParams{price: c.price + BUMP, ..*c}};
    let diff = simulate_call(&bs_params, bump) / BUMP;
    assert!((diff - delta).abs() < EPS);

    let pdelta = black_scholes::put_delta(&bs_params);
    let pdiff = simulate_put(&bs_params, bump) / BUMP;
    assert!((pdiff - pdelta).abs() < EPS);
}

#[test]
fn test_gamma() {
    let bs_params = BlackScholesParams {
        price: 4801.0,
        div_yield: 0.03,
        strike: 4900.0,
        vol: 0.2,
        rate: 0.02,
        time_to_expiry: 0.01
    };

    let gamma = black_scholes::gamma(&bs_params);
    assert!((gamma - 0.002).abs() < 0.001);
}

#[test]
fn test_rho() {
    let bs_params = BlackScholesParams {
        price: 873.0,
        div_yield: 0.01,
        strike: 950.0,
        vol: 0.4,
        rate: 0.08,
        time_to_expiry: 1.0
    };

    let rho = black_scholes::call_rho(&bs_params);
    assert!((rho - 3.56567).abs() < 0.001);

    let diff = simulate_call(&bs_params, |c| {BlackScholesParams{rate: c.rate + BUMP, ..*c}}) / BUMP / 100.0;
    assert!((diff - rho).abs() < EPS);
}

#[test]
fn test_vega() {
    let bs_params = BlackScholesParams {
        price: 290.0,
        div_yield: 0.1,
        strike: 250f64,
        vol: 0.15,
        rate: 0.03,
        time_to_expiry: 1.0
    };

    let vega = black_scholes::vega(&bs_params);
    assert!((vega - 0.93902) < 0.001);

    let diff = simulate_call(&bs_params, |c| {BlackScholesParams{vol: c.vol + BUMP, ..*c}}) / BUMP / 100.0;
    assert!((diff - vega).abs() < EPS);
}

#[test]
fn test_theta() {
    let bs_params = BlackScholesParams {
        price: 104.0,
        div_yield: 0.03,
        strike: 95.0,
        vol: 0.43,
        rate: 0.02,
        time_to_expiry: 1.0
    };

    let theta = black_scholes::call_theta(&bs_params);
    assert!((theta + 6.908).abs() < 0.001);

    let diff = simulate_call(&bs_params, |c| BlackScholesParams {time_to_expiry: c.time_to_expiry - BUMP, ..*c}) / BUMP;
    assert!((diff - theta).abs() < EPS);
}

#[test]
fn test_phi() {
    let bs_params = BlackScholesParams {
        price: 604.0,
        div_yield: 0.07,
        strike: 620.0,
        vol: 0.3,
        rate: 0.02,
        time_to_expiry: 1.0
    };

    let bump = |c: &BlackScholesParams| BlackScholesParams{div_yield: c.div_yield + BUMP, ..*c};

    let call_phi = black_scholes::call_phi(&bs_params);
    let call_diff = simulate_call(&bs_params, bump) / BUMP / 100.0;
    assert!((call_phi - call_diff).abs() < EPS);

    let put_phi = black_scholes::put_phi(&bs_params);
    let put_diff = simulate_put(&bs_params, bump) / BUMP / 100.0;
    assert!((put_phi - put_diff).abs() < EPS)
}

#[test]
fn test_call_impl_vol() {
    let bs_params = BlackScholesParams {
        price: 78.0,
        div_yield: 0.03,
        strike: 70.0,
        vol: 0.3,
        rate: 0.04,
        time_to_expiry: 1.0
    };

    let call_prem = black_scholes::call_premium(&bs_params);
    let iv = black_scholes::call_impl_vol(&call_prem, &bs_params);
    assert!((iv - bs_params.vol).abs() < 0.000001);
}
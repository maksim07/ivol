#![cfg(test)]
use crate::black_scholes;
use crate::black_scholes::BlackScholesParams;

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

    let call_premium1 = black_scholes::call_premium(&bs_params);
    let call_premium2 = black_scholes::call_premium(&BlackScholesParams{price: bs_params.price + 1.0, ..bs_params});
    assert!((call_premium2 - call_premium1 - delta).abs() < 0.01)
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

    let call_premium = black_scholes::call_premium(&bs_params);
    let rho = black_scholes::call_rho(&bs_params);
    assert!((rho - 3.56567).abs() < 0.001);

    let bs_params2 = BlackScholesParams{rate: bs_params.rate + 0.01, .. bs_params};
    let call_premium2 = black_scholes::call_premium(&bs_params2);
    assert!((call_premium2 - call_premium - rho).abs() < 0.1);
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

    let call_prem = black_scholes::call_premium(&bs_params);
    let vega = black_scholes::vega(&bs_params);
    assert!((vega - 0.93902) < 0.001);

    let bs_params_bumped = BlackScholesParams{
        vol: bs_params.vol + 0.01,
        ..bs_params
    };
    let call_prem_bumped = black_scholes::call_premium(&bs_params_bumped);
    assert!((call_prem + vega - call_prem_bumped).abs() < 0.01);

}
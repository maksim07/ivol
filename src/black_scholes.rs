//! module with Black/Scholes formula implementation and implied volatility calculation.
use std::f64::consts::PI;
use rv::prelude::*;
use nrfind::*;

/// Precision for Newton-Raphson method
const EPS : f64 = 0.0000001;
/// Maximum number of iterations for Newton-Raphson method
const ITER: i32 = 60000;

/// Parameters of BlackScholes
///
#[derive(Debug)]
pub struct BlackScholesParams {
    /// spot price of option underlying asset
    pub price: f64,
    /// annual dividend yield
    pub div_yield: f64,
    /// option strike
    pub strike: f64,
    /// volatility percent (in decimal)
    pub vol: f64,
    /// risk free rate
    pub rate: f64,
    /// time to expiry in years (decimal)
    pub time_to_expiry: f64
}


/// Function that calculates call option premium using Black/Scholes
///
pub fn call_premium(bs_params: &BlackScholesParams) -> f64 {
    generic_black_scholes(true, bs_params)
}

/// Function that calculates put option premium using Black/Scholes
///
pub fn put_premium(bs_params: &BlackScholesParams) -> f64 {
    generic_black_scholes(false, bs_params)
}

/// Returns call or put price based on put/call price using put/call parity formula
///
pub fn callput_price(is_call_price: bool, market_price: &f64, bs_params: &BlackScholesParams) -> f64 {
    let sign = if is_call_price {1.0} else {-1.0};
    let dprice = bs_params.price * (-bs_params.div_yield * bs_params.time_to_expiry).exp();
    let dstrike = bs_params.strike * (-bs_params.rate * bs_params.time_to_expiry).exp();
    market_price - sign * (dprice - dstrike)
}

/// Gives premium change for the call option simulation
///
/// * source - original contract
/// * bump - contract simulation function that creates new contract out of original
///
/// # Example
///
/// ```
/// use ivol::black_scholes::*;
///
/// const SMALL_BUMP: f64 = 0.00001;
/// let bs_params = BlackScholesParams{price: 345.0, strike: 330.0, div_yield: 0.06, rate: 0.025, vol: 0.34, time_to_expiry: 1.0};
///
/// let delta = simulate_call(&bs_params, |c| {BlackScholesParams{price: c.price + SMALL_BUMP, ..*c}}) / SMALL_BUMP;
/// let real_call_delta = call_delta(&bs_params);
/// assert!((delta - real_call_delta).abs() < 0.0001);
/// ```
///
/// The `delta` value in the example should be roughly the same as call's Delta sensitivity.
pub fn simulate_call<F>(source: &BlackScholesParams, bump: F) -> f64
    where F: Fn(&BlackScholesParams) -> BlackScholesParams {
    let call_prem = call_premium(source);
    let new_call = bump(source);
    let call_prem_bumped = call_premium(&new_call);
    call_prem_bumped - call_prem
}

/// Same as [`simulate_call`] but for put options
pub fn simulate_put<F>(source: &BlackScholesParams, bump: F) -> f64
    where F: Fn(&BlackScholesParams) -> BlackScholesParams {
    let put_prem = put_premium(source);
    let new_put = bump(source);
    let put_prem_bumped = put_premium(&new_put);
    put_prem_bumped - put_prem
}


/// Delta sensitivity for call options
///
pub fn call_delta(bs_params: &BlackScholesParams) -> f64 {
    generic_delta(true, bs_params)
}

/// Delta sensitivity for put options
///
pub fn put_delta(bs_params: &BlackScholesParams) -> f64 {
    generic_delta(false, bs_params)
}

/// Gamma sensitivity for call/put options
///
pub fn gamma(bs_params: &BlackScholesParams) -> f64 {
    let n: Gaussian = Gaussian::standard();
    let d1 = d1(bs_params);

    ((-bs_params.div_yield * bs_params.time_to_expiry).exp() / (bs_params.price * bs_params.vol * bs_params.time_to_expiry.sqrt())) * n.pdf(&d1)
}

/// Rho sensitivity for call options
///
pub fn call_rho(bs_params: &BlackScholesParams) -> f64 {
    generic_rho(true, bs_params)
}

/// Rho sensitivity for put options
///
pub fn put_rho(bs_params: &BlackScholesParams) -> f64 {
    generic_rho(false, bs_params)
}

/// Theta sensitivity for call options
pub fn call_theta(bs_params: &BlackScholesParams) -> f64 {
    generic_theta(true, bs_params)
}

/// Theta sensitivity for put options
pub fn put_theta(bs_params: &BlackScholesParams) -> f64 {
    generic_theta(false, bs_params)
}

/// Phi (dividend yield risk) calculation for call options
///
pub fn call_phi(bs_params: &BlackScholesParams) -> f64 {
    generic_phi(true, bs_params)
}

/// Phi (dividend yield risk) calculation for put options
///
pub fn put_phi(bs_params: &BlackScholesParams) -> f64 {
    generic_phi(false, bs_params)
}

/// Function that calculates option's Vega
///
pub fn vega(bs_params: &BlackScholesParams) -> f64 {
    0.01 * dtv_dvol(bs_params)
}

/// Calculates implied volatility from call market price and other option parameters
///
pub fn call_impl_vol(call_market_price: &f64, bs_params: &BlackScholesParams) -> Result<f64, f64> {
    let func = |v: f64| call_premium(&BlackScholesParams{vol: v, ..*bs_params}) - *call_market_price;
    let vol_deriv = |v: f64| dtv_dvol(&BlackScholesParams{vol: v, ..*bs_params});
    let vol_guess = approx_vol(call_market_price, bs_params);
    let root = find_root(&func, &vol_deriv, vol_guess, EPS, ITER);
    process_impl_vol_result(root, &vol_guess)
}

/// Calculates implied volatility from put market price and other option parameters
///
pub fn put_impl_vol(put_market_price: &f64, bs_params: &BlackScholesParams) -> Result<f64,f64> {
    let func = |v: f64| put_premium(&BlackScholesParams { vol: v, ..*bs_params }) - *put_market_price;
    let vol_deriv = |v: f64| dtv_dvol(&BlackScholesParams { vol: v, ..*bs_params });
    let call_market_price = callput_price(false, put_market_price, bs_params);
    let vol_guess = approx_vol(&call_market_price, bs_params);
    let root = find_root(&func, &vol_deriv, vol_guess, EPS, ITER);
    process_impl_vol_result(root, &vol_guess)
}

#[inline]
fn process_impl_vol_result(r: Result<f64, f64>, vol_guess: &f64) -> Result<f64, f64> {
    match r {
        Ok(v) => Ok(v),
        Err(e) => {
            if e.is_normal() {Err(e)} else {Err(*vol_guess)}
        }
    }
}

/// Calculates approximate volatility based on price, strike, rate and time to maturity (Corrado and Miller).
/// This is initial value for Newton-Raphson method.
fn approx_vol(market_price: &f64, bs_params: &BlackScholesParams) -> f64 {
    let dprice = bs_params.price * (-bs_params.div_yield * bs_params.time_to_expiry).exp();
    let dstrike = bs_params.strike * (-bs_params.rate * bs_params.time_to_expiry).exp();
    let psdiff = dprice - dstrike;
    let mdiff = market_price - psdiff / 2.0;
    let under_root = mdiff.powi(2) - psdiff.powi(2) / PI;
    let root = if under_root > 0.0 {under_root.sqrt()} else {0.0};
    let sum = (2.0 * PI).sqrt() / (dprice + dstrike) * (mdiff + root);

    sum / bs_params.time_to_expiry.sqrt()
}

#[test]
fn test_appox_vol() {
    let bs_params = BlackScholesParams {
        price: 200.0, strike: 180.0, vol: 0.34, div_yield: 0.05, rate: 0.08, time_to_expiry: 0.5
    };

    let premium = call_premium(&bs_params);
    let bumped_premium = premium * 1.05;

    let vol_approx = approx_vol(&bumped_premium, &bs_params);
    assert!((vol_approx - 0.34).abs() < 0.1);
}

/// Generic function for phi (dividend yield risk) sensitivity calculation for options
///
#[inline]
fn generic_phi(is_call: bool, bs_params: &BlackScholesParams) -> f64 {
    0.01 * dtv_ddiv(is_call, bs_params)
}

/// Black/Scholes derivative with respect to divs
///
pub fn dtv_ddiv(is_call: bool, bs_params: &BlackScholesParams) -> f64 {
    let n: Gaussian = Gaussian::standard();
    let sign = if is_call {1.0} else {-1.0};
    let dprice = bs_params.price * (-bs_params.div_yield * bs_params.time_to_expiry).exp();
    let d1 = sign * d1(bs_params);
    -sign * bs_params.time_to_expiry * dprice * n.cdf(&d1)
}

/// Generic theta calculation function
///
#[inline]
fn generic_theta(is_call: bool, bs_params: &BlackScholesParams) -> f64 {
    let n: Gaussian = Gaussian::standard();
    let sign = if is_call {1.0} else {-1.0};
    let dprice =  bs_params.price * (-bs_params.div_yield * bs_params.time_to_expiry).exp();
    let dstrike =  bs_params.strike * (-bs_params.rate * bs_params.time_to_expiry).exp();
    let d1 = d1(bs_params);
    let d2 = d2(bs_params);
    let v = bs_params.vol / (2.0 * bs_params.time_to_expiry.sqrt());

    sign * bs_params.div_yield * dprice * n.cdf(&(sign * d1)) - sign * bs_params.rate * dstrike * n.cdf(&(sign * d2)) - dprice * v * n.pdf(&d1)
}

/// Generic function for calculating Rho
///
#[inline]
fn generic_rho(is_call: bool, bs_params: &BlackScholesParams) -> f64 {
    0.01 * dtv_drate(is_call, bs_params)
}

/// Black/Scholes derivative with respect to interest rate
///
#[inline]
pub fn dtv_drate(is_call: bool, bs_params: &BlackScholesParams) -> f64 {
    let sign = if is_call {1.0} else {-1.0};
    let n: Gaussian = Gaussian::standard();
    let cdf_arg = sign * d2(bs_params);
    sign * bs_params.strike * (-bs_params.rate * bs_params.time_to_expiry).exp() * bs_params.time_to_expiry * n.cdf(&cdf_arg)
}

/// Delta calculation for both put and calls
///
#[inline]
fn generic_delta(is_call: bool, bs_params: &BlackScholesParams) -> f64 {
    let sign = if is_call {1.0} else {-1.0};
    let n: Gaussian = Gaussian::standard();
    let d1 = d1(bs_params);
    sign * (-bs_params.div_yield * bs_params.time_to_expiry).exp() * n.cdf(&(sign * d1))
}


/// BlackScholes derivative for volatility
///
#[inline]
pub fn dtv_dvol(bs_params: &BlackScholesParams) -> f64 {
    let n: Gaussian = Gaussian::standard();
    let d1 = d1(bs_params);
    if !d1.is_normal() {
        println!("Not normal d1 {:?}", bs_params);
    }
    bs_params.price * (-bs_params.div_yield * bs_params.time_to_expiry).exp() * n.pdf(&d1) * bs_params.time_to_expiry.sqrt()
}

/// Generic Black/Scholes calculation for both call and put options
///
#[inline]
fn generic_black_scholes(is_call: bool, bs_params: &BlackScholesParams) -> f64 {
    let sign = if is_call {1.0} else {-1.0};
    let n: Gaussian = Gaussian::standard();
    let d1 = sign * d1(bs_params);
    let d2 = sign * d2(bs_params);
    let d = (-bs_params.rate * bs_params.time_to_expiry).exp();
    let dd = (-bs_params.div_yield * bs_params.time_to_expiry).exp();

    sign * bs_params.price * n.cdf(&d1) * dd - sign * bs_params.strike * n.cdf(&d2) * d
}

/// D1 sub-formula of Black/Scholes
///
#[inline]
fn d1(bs_params: &BlackScholesParams) -> f64 {
    let mut d : f64 = bs_params.price / bs_params.strike;
    d = d.ln() + (bs_params.rate - bs_params.div_yield + bs_params.vol * bs_params.vol / 2.0) * bs_params.time_to_expiry;
    d / (bs_params.vol * bs_params.time_to_expiry.sqrt())
}

/// D2 sub-formula of Black/Scholes
///
#[inline]
fn d2(bs_params: &BlackScholesParams) -> f64 {
    d1(bs_params) - bs_params.vol * bs_params.time_to_expiry.sqrt()
}
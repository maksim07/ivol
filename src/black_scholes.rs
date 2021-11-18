//! module with BlackScholes implementation
use rv::prelude::*;

/// Parameters of BlackScholes
///
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
/// const EPS: f64 = 0.0001;
/// let bs_params = BlackScholesParams{price: 345.0, strike: 330.0, div_yield: 0.06, rate: 0.025, vol: 0.34, time_to_expiry: 1.0};
///
/// let delta = simulate_call(&bs_params, |c| {BlackScholesParams{price: c.price + EPS, ..*c}}) / EPS;
/// let real_call_delta = call_delta(&bs_params);
/// assert!((delta - real_call_delta).abs() < EPS);
/// ```
///
/// The `delta` value in the example should be roughly the same as call Delta sensitivity.
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

/// Generic function for calculating Rho
///
fn generic_rho(is_call: bool, bs_params: &BlackScholesParams) -> f64 {
    let sign = if is_call {1.0} else {-1.0};
    let n: Gaussian = Gaussian::standard();
    let cdf_arg = sign * d2(bs_params);
    0.01 * sign * bs_params.strike * (-bs_params.rate * bs_params.time_to_expiry).exp() * bs_params.time_to_expiry * n.cdf(&cdf_arg)
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

/// Function that calculates option's Vega
///
pub fn vega(bs_params: &BlackScholesParams) -> f64 {
    let n: Gaussian = Gaussian::standard();
    let d1 = d1(bs_params);
    0.01 * bs_params.price * (-bs_params.div_yield * bs_params.time_to_expiry).exp() * n.pdf(&d1) * bs_params.time_to_expiry.sqrt()
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



//! module with BlackScholes implementation
use rv::prelude::*;

///
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


///
/// Function that calculates call option premium using Black/Scholes
///
pub fn call_premium(bs_params: &BlackScholesParams) -> f64 {
    generic_black_scholes(true, bs_params)
}

///
/// Function that calculates put option premium using Black/Scholes
///
pub fn put_premium(bs_params: &BlackScholesParams) -> f64 {
    generic_black_scholes(false, bs_params)
}

///
/// Function that calculates option's Vega
///
pub fn vega(bs_params: &BlackScholesParams) -> f64 {
    let n: Gaussian = Gaussian::standard();
    let d1 = d1(bs_params);
    0.01 * bs_params.price * (-bs_params.div_yield * bs_params.time_to_expiry).exp() * n.pdf(&d1) * bs_params.time_to_expiry.sqrt()
}

///
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

///
/// D1 sub-formula of Black/Scholes
///
#[inline]
fn d1(bs_params: &BlackScholesParams) -> f64 {
    let mut d : f64 = bs_params.price / bs_params.strike;
    d = d.ln() + (bs_params.rate - bs_params.div_yield + bs_params.vol * bs_params.vol / 2.0) * bs_params.time_to_expiry;
    d / (bs_params.vol * bs_params.time_to_expiry.sqrt())
}

///
/// D2 sub-formula of Black/Scholes
///
#[inline]
fn d2(bs_params: &BlackScholesParams) -> f64 {
    d1(bs_params) - bs_params.vol * bs_params.time_to_expiry.sqrt()
}

///
/// Tests
///
#[cfg(test)]
mod test {
    use super::BlackScholesParams;

    #[test]
    fn test_call_premium() {
        let mut cp = super::call_premium(&BlackScholesParams {
            price: 100f64,
            div_yield: 0f64,
            strike: 110f64,
            vol: 0.12,
            rate: 0.05,
            time_to_expiry: 1.5
        });
        assert!((cp - 4.95).abs() < 0.01);

        cp = super::call_premium(&BlackScholesParams {
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
        let mut cp = super::put_premium(&BlackScholesParams {
            price: 150f64,
            div_yield: 0.1,
            strike: 200f64,
            vol: 0.11,
            rate: 0.01,
            time_to_expiry: 1.5
        });
        assert!((cp - 67.92).abs() < 0.01);

        cp = super::put_premium(&BlackScholesParams {
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
    fn test_vega() {
        let bs_params = BlackScholesParams {
            price: 290.0,
            div_yield: 0.1,
            strike: 250f64,
            vol: 0.15,
            rate: 0.03,
            time_to_expiry: 1.0
        };

        let call_prem = super::call_premium(&bs_params);
        let vega = super::vega(&bs_params);
        assert!((vega - 0.93902) < 0.001);

        let bs_params_bumped = BlackScholesParams{
            vol: bs_params.vol + 0.01,
            ..bs_params
        };
        let call_prem_bumped = super::call_premium(&bs_params_bumped);
        assert!((call_prem + vega - call_prem_bumped).abs() < 0.01);

    }
}

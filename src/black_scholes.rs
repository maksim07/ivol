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
    pub time_to_exp: f64
}


///
/// Method that calculates call option premium using Black/Scholes
///
pub fn call_premium(bs_params: &BlackScholesParams) -> f64 {
    let n: Gaussian = Gaussian::standard();
    let d1 = d1(bs_params);
    let d2 = d2(bs_params);
    let d = (-bs_params.rate * bs_params.time_to_exp).exp();
    let dd = (-bs_params.div_yield * bs_params.time_to_exp).exp();

    bs_params.price * n.cdf(&d1) * dd - bs_params.strike * n.cdf(&d2) * d
}

///
/// Method that calculates put option premium using Black/Scholes
///
pub fn put_premium(bs_params: &BlackScholesParams) -> f64 {
    let n: Gaussian = Gaussian::standard();
    let d1 = -d1(bs_params);
    let d2 = -d2(bs_params);
    let d = (-bs_params.rate * bs_params.time_to_exp).exp();
    let dd = (-bs_params.div_yield * bs_params.time_to_exp).exp();

    bs_params.strike * n.cdf(&d2) * d - bs_params.price * n.cdf(&d1) * dd
}


///
/// D1 sub-formula of Black/Scholes
///
#[inline]
fn d1(bs_params: &BlackScholesParams) -> f64 {
    let mut d : f64 = bs_params.price / bs_params.strike;
    d = d.ln() + (bs_params.rate - bs_params.div_yield + bs_params.vol * bs_params.vol / 2f64) * bs_params.time_to_exp;
    d / (bs_params.vol * bs_params.time_to_exp.sqrt())
}

///
/// D2 sub-formula of Black/Scholes
///
#[inline]
fn d2(bs_params: &BlackScholesParams) -> f64 {
    d1(bs_params) - bs_params.vol * bs_params.time_to_exp.sqrt()
}

#[cfg(test)]
mod test {
    use crate::black_scholes::BlackScholesParams;

    #[test]
    fn test_call_premium() {
        let mut cp = super::call_premium(&BlackScholesParams {
            price: 100f64,
            div_yield: 0f64,
            strike: 110f64,
            vol: 0.12,
            rate: 0.05,
            time_to_exp: 1.5
        });
        assert!((cp - 4.95).abs() < 0.01);

        cp = super::call_premium(&BlackScholesParams {
            price: 100f64,
            div_yield: 0.01,
            strike: 110f64,
            vol: 0.12,
            rate: 0.05,
            time_to_exp: 1.5
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
            time_to_exp: 1.5
        });
        assert!((cp - 67.92).abs() < 0.01);

        cp = super::put_premium(&BlackScholesParams {
            price: 150f64,
            div_yield: 0.1,
            strike: 160f64,
            vol: 0.06,
            rate: 0.03,
            time_to_exp: 0.5
        });
        assert!((cp - 14.95).abs() < 0.01);
    }
}

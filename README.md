# ivol

[![Tests status](https://github.com/maksim07/ivol/workflows/Rust/badge.svg?event=push)](https://github.com/maksim07/ivol/actions/workflows/rust.yml)
[![Build Status](https://app.travis-ci.com/maksim07/ivol.svg?branch=master)](https://app.travis-ci.com/maksim07/ivol)


Crate contains some financial API, which includes Black/Scholes formula, implied volatility calculation, european option sensitivities calculation.

## Examples

### **black_choles** module

```rust
use ivol::black_scholes::*;

let bs_params = BlackScholesParams {
    price: 4792.0,
    strike: 4400.0,
    rate: 0.025,
    time_to_expiry: 0.5,
    vol: 0.23,
    div_yield: 0.04
};

// call/put premium
let call_premium = call_premium(&bs_params);
let put_premium = put_premium(&bs_params);

// greeks
let call_delta = call_delta(&bs_params);
let put_delta = put_delta(&bs_params);
let vega = vega(&bs_params);
let put_rho = put_rho(&bs_params);
let call_rho = call_rho(&bs_params);
let call_theta = call_theta(&bs_params);
let put_theta = put_theta(&bs_params);
let gamma = gamma(&bs_params);

// dividend risk greeks
let call_phi = call_phi(&bs_params);
let put_phi = put_phi(&bs_params);


// let's calculate implied volatility from call option price
let iv = call_impl_vol(&call_premium, &bs_params).unwrap();
assert!((iv - bs_params.vol).abs() < 0.0000001);

// ... and from put price
let iv2 = put_impl_vol(&put_premium, &bs_params).unwrap();
assert!((iv2 - bs_params.vol).abs() < 0.0000001);
```

# Licence and version

* Current version: 0.0.2
* License: [Apache-2.0](LICENCE)

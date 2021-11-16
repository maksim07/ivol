use clap::{Arg, App};
use ivol::black_scholes::*;

fn main() {

    // float values validator
    let float_validator = |v: String| {
        match v.parse::<f64>() {
            Ok(_) => Ok(()),
            Err(_) => Err("value must be a float number".to_string())
        }
    };

    // parsing args
    let matches = App::new("ivol cli")
        .version("0.01")
        .about("Command Line Tool for ivol module functions")
        .arg(Arg::with_name("price")
            .short("p")
            .long("price")
            .value_name("PRICE")
            .help("option's underlying price")
            .required(true)
            .validator(float_validator)
            .takes_value(true))
        .arg(Arg::with_name("strike")
            .short("s")
            .long("strike")
            .value_name("STRIKE")
            .help("option's strike value")
            .required(true)
            .validator(float_validator)
            .takes_value(true))
        .arg(Arg::with_name("rate")
            .short("r")
            .long("rate")
            .value_name("RATE")
            .help("risk free rate")
            .required(true)
            .validator(float_validator)
            .takes_value(true))
        .arg(Arg::with_name("div_yield")
            .short("d")
            .long("div")
            .value_name("DIV")
            .help("annual dividend yield")
            .default_value("0")
            .validator(float_validator)
            .takes_value(true))
        .arg(Arg::with_name("time_to_expiry")
            .short("t")
            .long("time_to_expiry")
            .value_name("TIME")
            .help("time to expiry in years")
            .default_value("1")
            .validator(float_validator)
            .takes_value(true))
        .arg(Arg::with_name("volatility")
            .short("v")
            .long("volatility")
            .value_name("VOL")
            .help("volatility in percent (decimal)")
            .required(true)
            .validator(float_validator)
            .takes_value(true))
        .get_matches();

    // extracting values from args
    let price: f64 = matches.value_of("price").unwrap().parse().unwrap();
    let strike: f64 = matches.value_of("strike").unwrap().parse().unwrap();
    let vol: f64 = matches.value_of("volatility").unwrap().parse().unwrap();
    let rate: f64 = matches.value_of("rate").unwrap().parse().unwrap();
    let time_to_expiry: f64 = matches.value_of("time_to_expiry").unwrap().parse().unwrap();
    let div_yield: f64 = matches.value_of("div_yield").unwrap().parse().unwrap();

    let bs_params = BlackScholesParams {
        price,
        strike,
        rate,
        div_yield,
        vol,
        time_to_expiry
    };

    // calculating call and put option premiums
    let call_premium = call_premium(&bs_params);
    let put_premium = put_premium(&bs_params);
    let call_delta = call_delta(&bs_params);
    let put_delta = put_delta(&bs_params);
    let vega = vega(&bs_params);

    println!("Option call premium is {} and put premium is {} with the greeks:\n    \
              Call Delta = {}\n    \
              Put Delta = {}\n    \
              Vega = {}", &call_premium, &put_premium, &call_delta, &put_delta, &vega);
}
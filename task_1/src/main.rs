fn main() {
    let mut eps: f32 = 1.0;
    let mut mantissa_bits = 0;

    while 1.0 + (eps / 2.0) != 1.0 {
        eps /= 2.0;
        mantissa_bits += 1;
    }

    let mut max_exp = 0;
    let mut val: f32 = 1.0;

    while !val.is_infinite() {
        val *= 2.0;
        max_exp += 1;
    }

    let mut min_normal_exp = 0;
    let mut val: f32 = 1.0;

    while (val / 2.0).is_normal() {
        val /= 2.0;
        min_normal_exp -= 1;
    }

    let mut min_exp = 0;
    let mut val: f32 = 1.0;

    while (val / 2.0) != 0.0 {
        val /= 2.0;
        min_exp -= 1;
    }

    println!(
        "f32 eps: {:.23}\nf32 mantissa size: {}\nf32 max exp: {}\nf32 min normal exp: {}\nf32 min exp: {}",
        eps, mantissa_bits, max_exp, min_normal_exp, min_exp
    );
    assert!(eps == f32::EPSILON);
    // https://doc.rust-lang.org/std/primitive.f32.html#associatedconstant.MANTISSA_DIGITS
    // Note that the size of the mantissa in the bitwise representation is one smaller than this since the leading 1 is not stored explicitly.
    assert!(mantissa_bits == f32::MANTISSA_DIGITS - 1);
    assert!(max_exp == f32::MAX_EXP);
    // https://doc.rust-lang.org/std/primitive.f32.html#associatedconstant.MIN_EXP
    // 'One greater than the minimum possible normal power of 2 exponent for a significand bounded by 1 <= x < 2'
    assert!(min_normal_exp == f32::MIN_EXP - 1);

    let vals = [
        ("1.0", 1.0),
        ("1.0 + (eps/2.0)", 1.0 + (eps / 2.0)),
        ("1.0 + eps", 1.0 + eps),
        ("1.0 + (eps/2.0) + eps", 1.0 + (eps / 2.0) + eps),
    ];

    println!("{:<25} | {:<25}", "expr", "value");
    println!("{}", "-".repeat(53));

    for (desc, val) in vals {
        println!("{:<25} | {:.23}", desc, val);
    }

    println!();

    let mut eps: f64 = 1.0;
    let mut mantissa_bits = 0;

    while 1.0 + (eps / 2.0) != 1.0 {
        eps /= 2.0;
        mantissa_bits += 1;
    }

    let mut max_exp = 0;
    let mut val: f64 = 1.0;

    while !val.is_infinite() {
        val *= 2.0;
        max_exp += 1;
    }

    let mut min_normal_exp = 0;
    let mut val: f64 = 1.0;

    while (val / 2.0).is_normal() {
        val /= 2.0;
        min_normal_exp -= 1;
    }

    let mut min_exp = 0;
    let mut val: f64 = 1.0;

    while (val / 2.0) != 0.0 {
        val /= 2.0;
        min_exp -= 1;
    }

    println!(
        "f64 eps: {:.52}\nf64 mantissa size: {}\nf64 max exp: {}\nf64 min normal exp: {}\nf64 min exp: {}",
        eps, mantissa_bits, max_exp, min_normal_exp, min_exp
    );

    assert!(eps == f64::EPSILON);
    // https://doc.rust-lang.org/std/primitive.f32.html#associatedconstant.MANTISSA_DIGITS
    // Note that the size of the mantissa in the bitwise representation is one smaller than this since the leading 1 is not stored explicitly.
    assert!(mantissa_bits == f64::MANTISSA_DIGITS - 1);
    assert!(max_exp == f64::MAX_EXP);
    // https://doc.rust-lang.org/std/primitive.f64.html#associatedconstant.MIN_EXP
    // 'One greater than the minimum possible normal power of 2 exponent for a significand bounded by 1 <= x < 2'
    assert!(min_normal_exp == f64::MIN_EXP - 1);

    let vals = [
        ("1.0", 1.0),
        ("1.0 + (eps/2.0)", 1.0 + (eps / 2.0)),
        ("1.0 + eps", 1.0 + eps),
        ("1.0 + (eps/2.0) + eps", 1.0 + (eps / 2.0) + eps),
    ];

    println!("{:<25} | {:<54}", "expr", "value");
    println!("{}", "-".repeat(82));

    for (desc, val) in vals {
        println!("{:<25} | {:.52}", desc, val);
    }
}

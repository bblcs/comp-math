use std::f64::consts::PI;
use std::ops::{Add, Mul, Sub};

#[derive(Copy, Clone, Debug)]
struct Complex {
    re: f64,
    im: f64,
}

impl Add for Complex {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            re: self.re + rhs.re,
            im: self.im + rhs.im,
        }
    }
}

impl Sub for Complex {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            re: self.re - rhs.re,
            im: self.im - rhs.im,
        }
    }
}

impl Mul for Complex {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            re: self.re * rhs.re - self.im * rhs.im,
            im: self.re * rhs.im + self.im * rhs.re,
        }
    }
}

impl Complex {
    fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    fn exp(theta: f64) -> Self {
        Self {
            re: theta.cos(),
            im: theta.sin(),
        }
    }
}

fn dft(input: &[Complex], inverse: bool) -> Vec<Complex> {
    let n = input.len();
    let mut output = vec![Complex::new(0.0, 0.0); n];
    let sign = if inverse { 1.0 } else { -1.0 };

    for (k, out) in output.iter_mut().enumerate() {
        for (j, inp) in input.iter().enumerate() {
            let theta = sign * 2.0 * PI * (k as f64) * (j as f64) / (n as f64);
            let w = Complex::exp(theta);
            *out = out.add(inp.mul(w));
        }
        if inverse {
            out.re /= n as f64;
            out.im /= n as f64;
        }
    }
    output
}

// Au'(x) + Bu(x) = D
#[derive(Copy, Clone)]
struct BoundaryCondition {
    a: f64,
    b: f64,
    d: f64,
}

fn thomas<F>(
    x0: f64,
    xn: f64,
    n: usize,
    bc_left: BoundaryCondition,
    bc_right: BoundaryCondition,
    f: F,
) -> (Vec<f64>, Vec<f64>)
where
    F: Fn(f64) -> f64,
{
    let h = (xn - x0) / (n as f64);
    let mut x = vec![0.0; n + 1];
    for (i, el) in x.iter_mut().enumerate() {
        *el = x0 + (i as f64) * h;
    }

    let mut alpha = vec![0.0; n + 1];
    let mut beta = vec![0.0; n + 1];

    let b0 = bc_left.a - bc_left.b * h;
    let c0 = bc_left.a;
    let d0 = bc_left.d * h + bc_left.a * h * h / 2.0 * f(x[0]);

    alpha[1] = c0 / b0;
    beta[1] = -d0 / b0;

    for i in 1..n {
        let ai = 1.0;
        let bi = 2.0;
        let ci = 1.0;
        let di = h * h * f(x[i]);

        let den = bi - ai * alpha[i];
        alpha[i + 1] = ci / den;
        beta[i + 1] = (ai * beta[i] - di) / den;
    }

    let mut y = vec![0.0; n + 1];
    let num = bc_right.d * h - bc_right.a * h * h / 2.0 * f(x[n]) + bc_right.a * beta[n];
    let den = bc_right.a + bc_right.b * h - bc_right.a * alpha[n];
    y[n] = num / den;

    for i in (0..n).rev() {
        y[i] = alpha[i + 1] * y[i + 1] + beta[i + 1];
    }
    (x, y)
}

fn solve_fourier(n: usize) -> (Vec<f64>, Vec<f64>) {
    let a = 0.0;
    let b = 2.0 * PI;
    let h = (b - a) / (n as f64);

    let mut x = vec![0.0; n];
    let mut f_vals = vec![Complex::new(0.0, 0.0); n];

    for i in 0..n {
        x[i] = a + (i as f64) * h;
        f_vals[i] = Complex::new(x[i].sin(), 0.0);
    }

    let f_hat = dft(&f_vals, false);
    let mut y_hat = vec![Complex::new(0.0, 0.0); n];

    for i in 1..n {
        let mut k = i as f64;
        if i > n / 2 {
            k = (i as f64) - (n as f64);
        }
        let k2 = k * k;
        y_hat[i].re = -f_hat[i].re / k2;
        y_hat[i].im = -f_hat[i].im / k2;
    }
    y_hat[0] = Complex::new(0.0, 0.0);

    let y_complex = dft(&y_hat, true);
    let y = y_complex.into_iter().map(|c| c.re).collect();

    (x, y)
}

fn calc_max_error(x: &[f64], y: &[f64], exact_f: impl Fn(f64) -> f64) -> f64 {
    x.iter()
        .zip(y.iter())
        .map(|(&xi, &yi)| (yi - exact_f(xi)).abs())
        .fold(0.0, f64::max)
}

fn main() {
    let test_set = &[
        1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192,
    ];
    let exact_sol = |x: f64| -x.sin();
    let conditions = vec![
        (
            "function left right",
            BoundaryCondition {
                a: 0.,
                b: 1.,
                d: 0.,
            },
            BoundaryCondition {
                a: 0.,
                b: 1.,
                d: 0.,
            },
        ),
        (
            "function left, derivative right",
            BoundaryCondition {
                a: 0.,
                b: 1.,
                d: 0.,
            },
            BoundaryCondition {
                a: 1.,
                b: 0.,
                d: 1.,
            },
        ),
        (
            "derivative left, function right",
            BoundaryCondition {
                a: 1.,
                b: 0.,
                d: -1.,
            },
            BoundaryCondition {
                a: 0.,
                b: 1.,
                d: 0.,
            },
        ),
    ];

    println!("\n====================thomas==================");
    for (name, bc_l, bc_r) in conditions {
        println!("\ncond: {}", name);
        println!("{:>4} | {:>8} | {:>8}", "n", "max err", "times");

        let mut prev_err = 0.0;
        for &n in test_set {
            let (x, y) = thomas(0.0, PI, n, bc_l, bc_r, |x| x.sin());
            let err = calc_max_error(&x, &y, exact_sol);

            if prev_err > 0.0 {
                println!("{:>4} | {:>8.2e} | {:>8.4}", n, err, prev_err / err);
            } else {
                println!("{:>4} | {:>8.2e} | {:>8}", n, err, "-");
            }
            prev_err = err;
        }
    }

    println!("\n===================fourier==================");
    println!("{:>4} | {:>8}", "n", "max err");
    for &n in test_set {
        let (x, y) = solve_fourier(n);
        let err = calc_max_error(&x, &y, exact_sol);
        println!("{:>4} | {:>8.2e}", n, err);
    }
}

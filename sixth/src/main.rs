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
            *out = out.add(inp.mul(Complex::exp(
                sign * 2.0 * PI * (k as f64) * (j as f64) / (n as f64),
            )));
        }
        if inverse {
            out.re /= n as f64;
            out.im /= n as f64;
        }
    }
    output
}

/// Au'(x) + Bu(x) = D
#[derive(Copy, Clone)]
struct BoundaryCondition {
    a: f64,
    b: f64,
    d: f64,
}

/// d/dt [g(t) du/dt] + h(t) du/dt + s(t) u = f(t)
#[derive(Clone, Copy)]
struct BVP2<G, H, S, F> {
    g: G,
    h: H,
    s: S,
    f: F,
}

fn thomas<G, H, S, F>(
    bvp: BVP2<G, H, S, F>,
    bc_left: BoundaryCondition,
    bc_right: BoundaryCondition,
    x0: f64,
    xn: f64,
    n: usize,
) -> (Vec<f64>, Vec<f64>)
where
    G: Fn(f64) -> f64,
    H: Fn(f64) -> f64,
    S: Fn(f64) -> f64,
    F: Fn(f64) -> f64,
{
    let tau = (xn - x0) / (n as f64);
    let mut x = vec![0.0; n + 1];
    for (i, ix) in x.iter_mut().enumerate() {
        *ix = x0 + (i as f64) * tau;
    }

    //  a_i * u_{i-1} - b_i * u_i + c_i * u_{i+1} = d_i
    let mut a = vec![0.0; n + 1];
    let mut b = vec![0.0; n + 1];
    let mut c = vec![0.0; n + 1];
    let mut d = vec![0.0; n + 1];

    a[0] = 0.0;
    b[0] = bc_left.a / tau - bc_left.b;
    c[0] = bc_left.a / tau;
    d[0] = bc_left.d;

    for i in 1..n {
        let ti = x[i];

        let pre_g = (bvp.g)(ti - tau / 2.0);
        let post_g = (bvp.g)(ti + tau / 2.0);
        let hi = (bvp.h)(ti);
        let si = (bvp.s)(ti);

        a[i] = pre_g / (tau * tau) - hi / (2.0 * tau);
        c[i] = post_g / (tau * tau) + hi / (2.0 * tau);
        b[i] = a[i] + c[i] - si;
        d[i] = (bvp.f)(ti);
    }

    a[n] = bc_right.a / tau;
    b[n] = bc_right.a / tau + bc_right.b;
    c[n] = 0.0;
    d[n] = -bc_right.d;

    let mut alpha = vec![0.0; n + 2];
    let mut beta = vec![0.0; n + 2];

    for i in 0..=n {
        let den = b[i] - a[i] * alpha[i];
        alpha[i + 1] = c[i] / den;
        beta[i + 1] = (a[i] * beta[i] - d[i]) / den;
    }

    let mut u = vec![0.0; n + 1];

    u[n] = beta[n + 1];

    for i in (0..n).rev() {
        u[i] = alpha[i + 1] * u[i + 1] + beta[i + 1];
    }

    (x, u)
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
    // y'' = sin(x)
    // d/dt [1 * du/dt] + 0 * du/dt + 0 * u = sin(x)
    let bvp = BVP2 {
        g: |_| 1.0,
        h: |_| 0.0,
        s: |_| 0.0,
        f: |x: f64| x.sin(),
    };

    for (name, bc_l, bc_r) in conditions {
        println!("\ncond: {}", name);
        println!("{:>4} | {:>8} | {:>8}", "n", "max err", "times");

        let mut prev_err = 0.0;
        for &n in test_set {
            let (x, y) = thomas(bvp, bc_l, bc_r, 0.0, PI, n);
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

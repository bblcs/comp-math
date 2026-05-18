let gauss f a b =
  let roots, weights =
    ( [|
        (* 8 *)
        -0.96028986;
        -0.79666648;
        -0.52553241;
        -0.18343464;
        0.18343464;
        0.52553241;
        0.79666648;
        0.96028986;
      |],
      [|
        0.10122854;
        0.22238103;
        0.31370665;
        0.36268378;
        0.36268378;
        0.31370665;
        0.22238103;
        0.10122854;
      |] )
  in
  let mid = (a +. b) /. 2.0 in
  let half = (b -. a) /. 2.0 in
  let rec sum i acc =
    if i >= 8 then acc
    else
      let x = mid +. (half *. roots.(i)) in
      sum (i + 1) (acc +. (weights.(i) *. f x))
  in
  half *. sum 0 0.0

let simpson f a b n =
  let h = (b -. a) /. float n in
  let rec sum i acc =
    if i >= n then acc
    else
      let x = a +. (float i *. h) in
      let t = f (x -. h) +. (4.0 *. f x) +. f (x +. h) in
      sum (i + 2) (acc +. t)
  in
  h /. 3.0 *. sum 1 0.0

let fa x =
  if x < 1e-10 then Float.pi
  else if x > 1.0 -. 1e-10 then 5.0 *. Float.pi
  else sin (Float.pi *. (x ** 5.0)) /. ((x ** 5.0) *. (1.0 -. x))

let fb x = exp (-.sqrt x +. sin (x /. 10.0))
let fbg t = fb (t /. (1.0 -. t)) *. (1.0 /. ((1.0 -. t) ** 2.0))

let rec runge f a b n prev_res =
  let res = simpson f a b n in
  let err = abs_float (res -. prev_res) /. 15.0 in
  if err < 1e-8 then (res, n) else runge f a b (n * 2) res

let () =
  let inita = simpson fa 0.0 1.0 2 in
  let ans, nodes = runge fa 0.0 1.0 4 inita in
  Printf.printf "I_a = %.10f (using %d intervals)\n" ans nodes;
  let ans = gauss fa 0.0 1.0 in
  Printf.printf "I_a = %.10f (gauss)\n" ans;
  let initb = simpson fb 0.0 2000.0 1000 in
  let ans, nodes = runge fb 0.0 2000.0 2500 initb in
  Printf.printf "I_b = %.10f (using %d intervals)\n" ans nodes;
  let ans = gauss fbg 0.0 1.0 in
  Printf.printf "I_b = %.10f (gauss with x = t / (1 - t))\n" ans;
  let rec loop n prev_err =
    if n <= 40000 then begin
      let res = simpson (fun x -> exp x) 0.0 1.0 n in
      let err = abs_float (exp 1.0 -. 1.0 -. res) in
      let ratio = if prev_err > 0.0 then prev_err /. err else 0.0 in
      if n = 2 then Printf.printf "n = %5d | err = %e\n" n err
      else
        Printf.printf "n = %5d | err = %e | changed by %.2f times\n" n err ratio;
      loop (n * 2) err
    end
  in
  loop 2 0.0;
  Printf.printf "eps %e"
    (let rec find eps =
       if 1.0 +. (eps /. 2.0) = 1.0 then eps else find (eps /. 2.0)
     in
     find 1.0)

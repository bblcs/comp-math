let f x = 1.0 /. (1.0 +. (25.0 *. (x ** 2.0)))
let idist n i = (2.0 *. float_of_int i /. float_of_int n) -. 1.0
let uni_nodes n = Array.init (n + 1) (idist n)

let cheba_nodes n =
  Array.init (n + 1) (fun k ->
      cos
        (float_of_int ((2 * k) + 1) *. Float.pi /. (2.0 *. float_of_int (n + 1))))

let div_diff xs ys =
  let n = Array.length xs - 1 in
  let coefs = Array.copy ys in
  for j = 1 to n do
    for i = n downto j do
      coefs.(i) <- (coefs.(i) -. coefs.(i - 1)) /. (xs.(i) -. xs.(i - j))
    done
  done;
  coefs

let newton xs coefs x =
  let n = Array.length xs in
  let rec loop i acc =
    if i < 0 then acc else loop (i - 1) (coefs.(i) +. ((x -. xs.(i)) *. acc))
  in
  loop (n - 2) coefs.(n - 1)

let collect_errors n use_chebyshev =
  let xs = if use_chebyshev then cheba_nodes n else uni_nodes n in
  let ys = Array.map f xs in
  let coefs = div_diff xs ys in

  let test_points = 1000 in
  List.init (test_points + 1) (fun i ->
      let x = idist test_points i in
      let px = newton xs coefs x in
      let err = abs_float (f x -. px) in
      (x, err))

let () =
  let oc = open_out "data.csv" in
  Printf.fprintf oc "n,type,x,err\n";
  for n = 3 to 10 do
    let err_uni = collect_errors n false in
    List.iter
      (fun (x, err) -> Printf.fprintf oc "%d,uni,%1.15e,%1.15e\n" n x err)
      err_uni;
    let err_cheba = collect_errors n true in
    List.iter
      (fun (x, err) -> Printf.fprintf oc "%d,cheba,%1.15e,%1.15e\n" n x err)
      err_cheba
  done;
  close_out oc

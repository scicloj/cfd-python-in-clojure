(ns cfd.one-d
  (:require
   [clojure.set]
   [fastmath.core :as fm :refer [exp PI pow]]
   [utils.num-clj :as num-clj]))

(defn get-dx
  "Compute the grid spacing dx based on parameters.
  Expects keys: :x-start, :x-end, and :nx (number of spatial points)."
  [{:keys [x-start x-end nx nt c]
    :or   {x-start 0.0
           x-end   2.0}}]
  (/ (- x-end x-start) (- nx 1)))

(defn get-dt
  [{:keys [sigma dx]
    :or   {sigma 0.5}}]
  (* sigma dx))

(defn create-array-x
  "Creates a float array of x coordinates."
  [{:keys [x-start x-end nx nt c]
    :or   {x-start 0
           x-end   2
           nx      50}
    :as   _params}]
  (num-clj/linspace {:start x-start
                     :stop  x-end
                     :num   nx}))

(defn create-array-u
  "Creates an initial float array of u using a condition function."
  [{:keys [array-x nx nt c condition-fn]
    :or   {condition-fn #(float (if (and (>= % 0.5) (<= % 1.0)) 2 1))}
    :as   params}]
  (let [array-x (or array-x (create-array-x params))
        nx      (or nx (alength array-x))
        u       (float-array nx)]
    (dotimes [i nx]
      (let [x-val (aget array-x i)
            u-val (float (condition-fn x-val))]
       (aset u i u-val)))
    u))

(defn update-convection-u
  "Performs one time-step update on the 1D array Convection u.

  For indices 1 to nx-1 the update u is(with given co-eff):
  1. when :linear:
   u[i] = u[i] - c * dt/dx * (u[i] - u[i-1])

   *note*: then parameter :c is required.

  2. when :nonlinear:
   u[i] = u[i] - u[i] * dt/dx * (u[i] - u[i-1])

  The parameter map (params) must include at least:
  - :dt (the time step)
  - :x-start, :x-end, and :nx (so that dx can be computed via get-dx)

  Optional ones in the param map:
  - :co-eff - either :linear(default) or :nonlinear"
  [array-u {:keys [c dx dt nx co-eff sigma]
            :or   {co-eff :linear
                   c      1.0}
            :as   params}]
  (let [co-eff-fn (case co-eff
                    :linear (constantly c)
                    :nonlinear identity
                    identity)
        dx    (or dx (get-dx params))
        dt    (or dt (get-dt {:sigma sigma :dx dx}))
        nx    (or nx (alength array-u))
        un    (float-array array-u)]
    (dotimes [i (- nx 2)]
      (let [i      (inc i)
            un-i   (aget un i)
            un-i-1 (aget un (dec i))]
        (aset array-u i (float (- un-i
                                  (* (co-eff-fn un-i)
                                     dt
                                     (/ 1 dx)
                                     (- un-i un-i-1)))))))
    (aset array-u 0 (aget un 1))
    (aset array-u (dec nx) (aget un (- nx 2)))
    array-u))

(defn update-diffusion-u
  "u[i] = u[i] + nu * dt / dx^2 * (u[i] + u[i-1] - 2 * u[i])"
  [array-u {:keys [sigma nu nt]
            :or   {sigma 0.2}}]
  (let [nx (alength array-u)
        dx (get-dx {:nx nx})
        dt (/ (* sigma (* dx dx)) nu)
        un (float-array array-u)]
    (dotimes [i (- nx 2)]
      (let [i   (inc i)
            idx (inc i)]
        (aset array-u i
          (float (+ (aget un i)
                    (* (/ (* nu dt) (* dx dx))
                       (- (+ (aget un idx)
                             (aget un (dec i)))
                          (* 2 (aget un i)))))))))
    array-u))

(defn get-burger-array-un [un-i+1 un-i un-i-1{:keys [dx nx nt nu dt nx]}]
  (float (+ un-i
            (- (* un-i
                  dt
                  (/ 1 dx)
                  (- un-i un-i-1)))
            (* nu
               dt
               (/ 1 (* dx dx))
               (+ un-i+1
                  (- (* 2 un-i))
                  un-i-1)))))

(defn update-computational-burger-u
  [array-u {:keys [dx nx nt nu dt nx] :as params}]
  (let [un (float-array array-u)]
    (dotimes [i (- nx 2)]
      (let [i (inc i)]
       (aset array-u i (get-burger-array-un
                         (aget un (inc i))
                         (aget un i)
                         (aget un (dec i))
                         params))))
    (aset array-u 0 (get-burger-array-un
                      (aget un 1)
                      (aget un 0)
                      (aget un (- nx 2))
                      params))
    (aset array-u (dec nx) (aget array-u 0))
    array-u))

(defn simulate
  "Runs the simulation for nt time steps.

  Parameters:
  - array-u: the initial u array
  - params: a map containing simulation parameters (including :nt and :dt)

  Returns the final u array after nt updates"
  [array-u {:keys [nt nx mode]
            :or   {mode :convection}
            :as   params}]
  (let [update-fn (case mode
                    :diffusion update-diffusion-u
                    :burger update-computational-burger-u
                    update-convection-u)]
   (loop [n 0]
     (if (= n nt)
       array-u
       (do (update-fn array-u params) (recur (inc n)))))))

(defn simulate-accumulate
  "Runs the simulation for nt time steps and accumulate those results.

  Parameters:
  - array-u: the initial u array
  - params: a map containing simulation parameters (including :nt and :dt)

  Returns the vector arrays of u array after nt updates"
  [array-u {:keys [nt nx mode]
            :or   {mode :convection}
            :as   params}]
  (let [update-fn (case mode
                    :diffusion update-diffusion-u
                    :burger update-computational-burger-u
                    update-convection-u)]
    (loop [n    0
           !cum (transient [(vec array-u)])]
      (if (= n nt)
        (persistent! !cum)
        (recur (inc n) (conj! !cum (vec (update-fn array-u params))))))))

;; Burgers' Equation Analytical Approach
;; todo: Refactor to implement symbolic math properly(i.e. Emmy)

(defn get-phi-first [{:keys [x nu t]}]
  (exp
    (/
     (- (pow (- x (* 4.0 t) (* 2.0 PI)) 2))
     (* (* 4.0 nu) (+ t 1.0)))))

(defn get-phi-second [{:keys [x nu t]}]
  (exp
    (/
     (- (pow (- x (* 4.0 t)) 2))
     (* (* 4.0 nu) (+ t 1.0)))))

(defn get-phi [{:keys [x nu t] :as params}]
  (+ (get-phi-first params)
     (get-phi-second params)))

(defn get-phi-prime [{:keys [x nu t] :as params}]
  (- (+
      (* (get-phi-second params)
         (/ 1.0 (* (* 4.0 nu) (+ t 1.0)))
         (- (* 2.0 x) (* 8.0 t)))
      (* (/ 1.0 (* (* 4.0 nu) (+ t 1.0)))
         (- (* 2.0 x) (* 8.0 t) (* 4.0 PI))
         (get-phi-first params)))))

(defn burgers-u [{:keys [x nu t] :as params}]
  (float (+ (- (* 2.0
                  nu
                  (/ (get-phi-prime params) (get-phi params))))
            4.0)))

(defn update-analytical-burger-u
  [{:keys [dx nx nt nu dt nx x-start x-end]
    :or   {x-start 0
           x-end   (* 2.0 PI)}}]
  (let [x       (num-clj/linspace {:start x-start :stop x-end :num nx})
        array-u (make-array Float nx)]
    (dotimes [i nx]
      (let [x-i (aget x i)]
        (aset array-u i (burgers-u {:t dt :x x-i :nu nu}))))
    array-u))

(comment)

(ns cfd.one-d
  (:refer-clojure :exclude [* + - /])
  (:require
   [fastmath.core :as fm :refer
    [* + - / exp long-div PI pow]]))

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
  "Creates an float array of x coordinates."
  [{:keys [x-start x-end nx nt c]
    :or   {x-start 0
           x-end   2}
    :as   params}]
  (let [dx  (get-dx params)
        arr (float-array nx)]
    (dotimes [i nx]
      (aset arr i (float (+ x-start (* i dx)))))
    arr))

(defn create-array-u
  "Creates an initial float array of u using a condition function."
  [{:keys [array-x nx nt c condition-fn]
    :or   {condition-fn #(float (if (and (>= % 0.5) (<= % 1.0)) 2 1))}
    :as   params}]
  (let [array-x (or array-x (create-array-x params))
        nx      (or nx (alength array-x))
        u       (float-array nx)]
    (dotimes [i nx]
      (aset u i (condition-fn (aget array-x i))))
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
  [array-u {:keys [c dt co-eff sigma]
            :or   {co-eff :linear}
            :as   params}]
  (let [co-eff-fn (case co-eff
                    :linear (constantly c)
                    (fn [{:keys [u i]}] (aget u i)))
        dx    (get-dx params)
        dt    (or dt (get-dt {:sigma sigma :dx dx}))
        nx    (alength array-u)
        un    (float-array array-u)]
    (dotimes [i (dec nx)]
      (let [idx (inc i)]
        (aset array-u idx (float (- (aget un idx)
                                    (* (co-eff-fn {:u un :i idx})
                                       (/ dt dx)
                                       (- (aget un idx) (aget un i))))))))
    array-u))

(defn update-diffusion-u
  [array-u {:keys [sigma nu nt]}]
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

(defn simulate
  "Runs the simulation for nt time steps.

  Parameters:
  - array-u: the initial u array
  - params: a map containing simulation parameters (including :nt and :dt)

  Returns the final u array after nt updates"
  [array-u {:keys [nt mode]
            :or   {mode :convection}
            :as   params}]
  (let [update-fn (case mode
                    :diffusion update-diffusion-u
                    update-convection-u)]
   (loop [n         0
          current-u array-u]
     (if (= n nt)
       current-u
       (recur (inc n) (update-fn current-u params))))))

;; Burgers' Equation
;; todo: Refactor to implement symbolic math properly(i.e. Emmy)

(defn linspace [{:keys [start stop num]}]
  (let [arr  (float-array num)
        step (/ (- stop start) (dec num))]
    (dotimes [i num]
      (aset arr i (float (* i step))))
    arr))

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

(defn update-burger-u
  [{:keys [dx wnx nt nu t]}]
  (let [dt      (* dx nu)
        x       (linspace {:start 0 :stop (* 2 PI) :num nx})
        array-u (make-array Float nx)]
    (dotimes [i nx]
      (let [x-i (aget x i)]
        (aset array-u i (burgers-u {:t t :x x-i :nu nu}))))
    array-u))

(comment)

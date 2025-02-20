(ns cfd.one-d)

(defn get-dx
  "Compute the grid spacing dx based on parameters.
  Expects keys: :x-start, :x-end, and :nx (number of spatial points)."
  [{:keys [x-start x-end nx nt c]
    :or   {x-start 0.0
           x-end   2.0}}]
  (/ (- x-end x-start) (- nx 1)))

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

(defn update-u
  "Performs one time-step update on the 1D array u.

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
  [array-u {:keys [c dt co-eff]
            :or   {co-eff :linear}
            :as   params}]
  (let [co-eff-fn (case co-eff
                    :linear (constantly c)
                    (fn [{:keys [u i]}] (aget u i)))
        dx    (get-dx params)
        nx    (alength array-u)
        un    (float-array array-u)]
    (dotimes [i (dec nx)]
      (let [idx (inc i)]
        (aset array-u idx (float (- (aget un idx)
                                    (* (co-eff-fn {:u un :i idx})
                                       (/ dt dx)
                                       (- (aget un idx) (aget un i))))))))
    array-u))

(defn simulate
  "Runs the simulation for nt time steps.

  Parameters:
  - array-u: the initial u array
  - params: a map containing simulation parameters (including :nt and :dt)

  Returns the final u array after nt updates"
  [array-u {:keys [nt] :as params}]
  (loop [n         0
         current-u array-u]
    (if (= n nt)
      current-u
      (recur (inc n) (update-u current-u params)))))

(comment)

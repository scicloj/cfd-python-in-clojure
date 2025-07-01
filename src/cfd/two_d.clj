(ns cfd.two-d
  (:require
   [tech.v3.datatype :as dt]
   [tech.v3.tensor :as dtt]
   [utils.num-clj :as num-clj]))

(defn ^:WIP create-2d-spacial-array->dt [{:keys [dt-type
                                                 nx x-start x-end
                                                 ny y-start y-end]
                                          :or   {dt-type :float32}}]
  {:array-x (dt/make-container dt-type
              (num-clj/linspace {:start x-start :stop x-end :num nx}))
   :array-y (dt/make-container dt-type
              (num-clj/linspace {:start y-start :stop y-end :num ny}))})

(defn create-2d-spacial-array [{:keys [nx x-start x-end
                                       ny y-start y-end]}]
  {:array-x (num-clj/linspace {:start x-start :stop x-end :num nx})
   :array-y (num-clj/linspace {:start y-start :stop y-end :num ny})})

(defn- default-condition-fn
  "Given x-val and y-val, creates initial flow velocity

  (0.5<=x<=1 && 0.5<=y<=1)=2 else 1"
  [x-val y-val]
  (if (and (< 0.5 x-val) (<= x-val 1)
           (< 0.5 y-val) (<= y-val 1))
    2 1))

(defn create-init-u
  "Create u(flow velocity) for x, y grid"
  [{:keys [nx ny array-x array-y condition-fn]
    :or {condition-fn default-condition-fn}}]
  (let [array-x-len     (or nx (alength array-x))
        array-y-len     (or ny (alength array-y))
        spatial-array-u (make-array Float/TYPE array-x-len array-y-len)]
    (dotimes [y-idx (- (alength spatial-array-u) 1)]
      (let [row-array (aget spatial-array-u y-idx)]
        (dotimes [x-idx (- (alength row-array) 1)]
          (let [x-val (aget array-x x-idx)
                y-val (aget array-y y-idx)
                u-val (condition-fn x-val y-val)]
            (aset spatial-array-u y-idx x-idx u-val)))))
    spatial-array-u))

(defn- update-convection-u
  [array-u {:keys [c dx dy dt nx ny]
            :or   {c 1.0}
            :as   params}]
  (let [un (object-array array-u)]
    (dotimes [y-idx (- ny 1)]
      (when (pos? y-idx)
        (dotimes [x-idx (- nx 1)]
          (when (pos? x-idx)
           (let [u-j-i   (aget un y-idx x-idx)
                 u-j-i-1 (aget un y-idx (dec x-idx))
                 u-j-1-i (aget un (dec y-idx) x-idx)]
             (aset array-u y-idx x-idx
               (float (- u-j-i
                         (* c
                            (/ dt dx)
                            (- u-j-i u-j-i-1))
                         (* c
                            (/ dt dy)
                            (- u-j-i u-j-1-i))))))))))

   ;; boundary condition
   (aset array-u 0 (float-array nx 1))
   (aset array-u (dec ny) (float-array nx 1))
   (dotimes [y-idx (- ny 1)]
     (aset array-u y-idx 0 (float 1))
     (aset array-u y-idx (dec nx) (float 1))))
  array-u)

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
                    :convection update-convection-u)]
    (loop [n 0]
      (if (= n nt)
        array-u
        (do (update-fn array-u params) (recur (inc n)))))))

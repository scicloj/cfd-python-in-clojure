^:kindly/hide-code
(ns conferences.scinoj-light-1.siyoung-talk
  (:require
   [cfd.one-d :as one-d]
   [fastmath.core :as fm :refer [exp sin PI pow]]
   [scicloj.kindly.v4.api :as kindly]
   [scicloj.kindly.v4.kind :as kind]
   [utils.notebook :refer [md tex cumulated-line-plot ts-vega-lite-plot-map vega-lite-plot]]))

;; ---
;; format:
;;   revealjs:
;;     incremental: false
;; ---
;;
;;
;;
;; ## 1D Viscous Fluid Flow Data Analysis Using Burgers’ Equation {style="font-size: 2.4rem; margin-top: 100px;"}
;; <p style="text-align: end; margin-top: 60px;">Siyoung B</p>
;;
;;
;;
;; ## About me
;; #### I am Siyoung _[she-young]_
;; <div style="display: flex; gap: 10px; width: 100%;">
;; <div style="transform: scale(0.7); margin: -70px 0 0 -140px;">
;; - A developer, loves Clojure and loves to fidget all day
;;    - cat, hiking, cycling, walking, knitting & sewing
;; - Likes to look up to see the night sky
;; - Studied astrophysics in undergrad
;;    - participated research of simulating binary stars' collision
;; </div>
;; <img src="https://avatars.githubusercontent.com/u/1319016?v=4" style="width:300px;" />
;; </div>
;;
;;
;;
;; ## Fluid Dynamics & CFD{style="font-size: 2rem; margin-bottom: 20px;"}
;;
;; #### Fluid Dynamics{style="margin-top: 40px;"}
;;
;; <div style="font-size: 1.4rem;">
;; - The physics of how liquids and gases move and flow
;; - **Many Variables**: Involves factors like friction, pressure, heat, and momentum etc.
;; - **Complex Phenomena**: Includes complicated movements like turbulence
;; - **Crucial Understanding**: Important for fields like airplane design, healthcare, and plumbing and even more!
;; </div>
;;
;; #### Computational Fluid Dynamics (CFD){style="margin-top: 40px;"}
;; <div style="font-size: 1.4rem;">
;; - Uses computers and math to simulate and predict fluid movement
;; - _**Clojure Gap**_: No CFD tools using Clojure yet(?)
;; - _**My Initiative**_: This project aims to create CFD tools using Clojure
;; </div>
;;
;;
;;
;; ## Past Research{style="font-size: 2rem; margin-bottom: 50px;"}
;;
;; <div style="display: flex; justify-content: space-between; width: 50%; height: 70%; gap: 20px;">
^:kindly/hide-code
(kind/fn {:sim-vid-src "notebooks/conferences/scinoj_light_1/resources/575n65_4plots.mp4"}
  {:kindly/f (fn [{:keys [sim-vid-src]}] (kind/video {:src sim-vid-src}))})
;;
^:kindly/hide-code
(kind/fn {:sim-vid-src "notebooks/conferences/scinoj_light_1/resources/105_100k.mp4"}
  {:kindly/f (fn [{:keys [sim-vid-src]}] (kind/video {:src sim-vid-src}))})
;;
;; </div>
;;
;;
;;
;; ## Where to start{style="font-size: 2rem;"}
;;
;; <div style="font-size: 1.4rem;">
;; - **Past CFD Experience**: As an end user, preparing, running then analyzing out of existing tools in astrophysics research
;; - **Knowledge Refresh**: Not much of formal CFD knowledge, and dated since university
;; </div>
;;
;; #### Learning Resource: <a href="https://github.com/barbagroup/CFDPython" target="_blank">CFD Python</a>{style="margin-top: 40px;"}
;; <div style="font-size: 1.4rem;">
;; - Utilizing Prof. Lorena Barba's "CFD Python" materials for relearning foundational CFD
;; - Developed at Boston University with Python code examples for teaching
;; - Broken down into 12 steps of learning materials to start with a simplified concept
;; - <a href="https://github.com/scicloj/cfd-python-in-clojure" target="_blank">CFD Python in Clojure</a>
;;    - Using "CFD Python" as a guide to implement CFD in Clojure
;;    - Currently in progress
;; </div>
;;
;;
;;
;; ## 1D Shock Interaction & Evolution{style="font-size: 2rem;"}
;;
^:kindly/hide-code (def nx 101)
^:kindly/hide-code (def nt 500)
^:kindly/hide-code (def nu 0.07)
^:kindly/hide-code (def dx (* 2.0 PI (/ 1 (- nx 1))))
^:kindly/hide-code (def dt (* dx nu))
^:kindly/hide-code (def x-start 0)
^:kindly/hide-code (def x-end (* 2.0 PI))
^:kindly/hide-code (def array-x (one-d/create-array-x {:x-start x-start :x-end x-end :nx nx}))
^:kindly/hide-code (def init-params {:nx      nx
                                     :nt      nt
                                     :dx      dx
                                     :dt      dt
                                     :nu      nu
                                     :array-x array-x
                                     :mode    :burger})
;;
;; <span style="font-size: 1.2rem;"> initial condition - Step function $u(x, 0) = \begin{cases} 1, \text{if } x < 2.0\\ 0 \text{, otherwise} \end{cases}$</span>
;; <div style="display: flex;">
^:kindly/hide-code (def init-step-u (one-d/create-array-u {:array-x      array-x
                                                           :condition-fn #(if (< % 2.0) 1 0)}))
^:kindly/hide-code (def cumulated-step-u (one-d/simulate-accumulate init-step-u init-params))
^:kindly/hide-code (def cum-plot-params (assoc init-params :cum-array-y cumulated-step-u))
^:kindly/hide-code (cumulated-line-plot cum-plot-params)
^:kindly/hide-code (vega-lite-plot {:plot-map (ts-vega-lite-plot-map [-0.1 1.1] cum-plot-params)})
;; </div>
;; <div style="font-size: 1.4rem; margin-top: -35px;">
;; - **Non-linear convection**: Faster fluid tends to "bunch up" with slower fluid in a complex way, steepening changes
;;    - **Non-linear**: Effects not always simple nor directly proportional
;;    - **Convection**: Movement affected by some changes i.e. heat, flow etc.
;; - **Viscosity**: The "stickiness" of the fluid that resists flow and smooths out speed differences
;; </div>
;;
;;
;;
;; ## Effect of Viscosity on Step Structure{style="font-size: 2rem;"}
;;
;; <div style="font-size: 1.4rem; margin-bottom: 20px;">
;; - **Viscosity** = _"Stickiness"_
;; - $\nu (=viscosity) = [0.03, 0.07, 0.1, 0.4]$
;; - Lower viscosity &rarr; sharper shock features
;; - Higher viscosity &rarr; increased smoothing of shock features
;; </div>
^:kindly/hide-code (def viscosity-vec [0.03 nu 0.1 0.4])
^:kindly/hide-code (def viscosity-result-arr-u-seq
                     (->> viscosity-vec
                          (map #(let [init-u       (float-array init-step-u)
                                      nu           %
                                      result-arr-u (one-d/simulate init-u (assoc init-params :nu nu :mode :burger))]
                                  (map (fn [x u] (hash-map :nu nu :x x :y u)) array-x result-arr-u)))
                          (apply concat)))
^:kindly/hide-code
(vega-lite-plot {:data          {:values (into [] viscosity-result-arr-u-seq)}
                 :x-label       "X"
                 :y-label       "Velocity"
                 :encoding-opts {:color {:field "nu" :type "nominal" :title "Kinematic Viscosity"}}
                 :height        260
                 :width         400})
;;
;;
;;
;; ## Implementation in Clojure 1{style="font-size: 2rem;"}
;;
;; #### Initial setup{style="margin-top: 20px;"}
;;
;; <div style="display: flex; justify-content: space-between; font-size: 1.4rem; margin-top: 15px;">
;; ```{.clojure style="width: 48%;"}
;; (defn create-initial-x [{:keys [x-start x-stop nx]}]
;;  (let [arr  (float-array nx)
;;        step (/ (- x-stop x-start) (dec nx))]
;;    (dotimes [i nx]
;;      (aset arr i (float (* i step))))
;;    arr))
;; ```
;; ```{.clojure style="width: 48%;"}
;; (defn create-initial-fluid-velocity
;;   [{:keys [array-x condition-fn] :as _params}]
;;   (let [nx      (alength array-x)
;;         array-u (float-array nx)]
;;     (dotimes [i nx]
;;       (let [x-val (aget array-x i)
;;             u-val (float (condition-fn x-val))]
;;         (aset array-u i u-val)))
;;     array-u))
;; ```
;; </div>
;;
;;
;; <div style="font-size: 1.4rem; margin-top: 40px;">
;; - Chose Java primitive arrays(float-array) for performance
;; - Mutable, non-persistent approach for speed and memory efficiency for future large-scale simulation
;; - **Pros**: Low memory overhead, Fast access and updates, Better suited for large numerical grid
;; - **Cons**: Breaks Clojure's idiomatic immutability, manual memory handling and index tracking, harder to debug and reason functionality
;; </div>
;;
;;
;;
;; ## Implementation in Clojure 2{style="font-size: 2rem;"}
;;
;; #### Python &rarr; Clojure{style="margin-top: 20px;"}
;; <div style="display: flex; justify-content: space-between;">
;;
;; ```{.python style="width: 50%; font-size: 1.4rem;"}
;; u = numpy.ones(nx)
;; u[int(.5 / dx):int(1 / dx + 1)] = 2
;; un = numpy.ones(nx)
;;
;; for n in range(nt):
;;     un = u.copy()
;;     for i in range(1, nx):
;;         u[i] = un[i] - c * dt / dx * (un[i] - un[i-1])
;; ```
;;
;; &rarr;
;;
;; ```{.clojure style="width: 50%; font-size: 1.4rem;"}
;; (def arr-u (create-initial-fluid-velocity init-params))
;;
;; (defn linear-convection-mode [idx arr-u un {:keys [c dt dx] :as _init-params}]
;;   (aset arr-u (inc idx)
;;         (float (- (aget un (inc idx)) (* c (/ dt dx) (- (aget un (inc idx)) (aget un idx)))))))
;;
;; (defn update-u [arr-u {:keys [c dt dx mode] :as params}]
;;   (let [un        (float-array arr-u)
;;         update-fn (case mode
;;                     :linear-convection linear-convection-mode
;;                     :burgers burgers-mode
;;                     other-modes)]
;;    (dotimes [idx (dec nx)]
;;      (update-fn idx arr-u un params))
;;    arr-u))
;;
;; (defn simulate [arr-u {:keys [nt] :as init-params}]
;;   (loop [n 0]
;;     (if (= n nt)
;;       arr-u
;;       (do (update-u arr-u init-params) (recur (inc n))))))
;; ```
;; </div>
;;
;;
;;
;; ## Implementation in Clojure 3{style="font-size: 2rem;"}
;; #### Mathematical Equation &rarr; Computational Equation{style="margin-top: 20px;"}
;;
;; <div style="display: flex; flex-direction: column; align-items: center;">
;; <div style="transform: scale(0.7); margin: -70px 0; width: 100vw;">
(tex "\\frac{\\partial u}{\\partial t} + u \\frac{\\partial u}{\\partial x} = \\nu\\frac{\\partial^2 u}{\\partial x^2}")
;; &rarr; _discretization, forward difference for time, backward difference for space etc...._ &rarr;
(tex
  (str "  u_i^{n+1} = u_i^n - u_i^n \\frac{\\Delta t}{\\Delta x}(u_i^n - u_{i-1}^n) + "
       "\\nu\\frac{\\Delta t}{\\Delta x^2}(u_{i+1}^n + u_{i-1}^n - 2u_i^n)"))
;; </div>
;; &darr;
;; ```{.clojure style="font-size: 1.4rem;"}
;; (defn get-burgers-arr-un [un-i+1 un-i un-i-1 {:keys [nu dx dt] :as params}]
;;   (float (+ un-i
;;             (- (* un-i dt (/ 1 dx) (- un-i un-i-1)))
;;             (* nu dt (/ 1 (* dx dx)) (+ un-i+1 (- (* 2 un-i)) un-i-1)))))
;; ```
;; </div>
;;
;;
;;
;; ## What's next{style="font-size: 2rem;"}
;;
;; <div style="font-size: 1.6rem;">
;; - Extend from 1D to 2D/3D
;; - Introduce pressure terms to meet Navier-Stokes equations
;; - Add boundary conditions and validation checks
;; - Explore further to implement w/ Clojure's functional and idiomatic way
;; - Scale up and test
;; </div>
;;
;;
;;
;; ## Thank you :) {style="text-align: center; margin-top: 200px;"}
;; ### Questions?{style="text-align: center; margin: 40px 0;"}
;;
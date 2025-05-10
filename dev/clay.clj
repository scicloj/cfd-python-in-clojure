(ns clay
  (:require
   [scicloj.clay.v2.api :as clay]))

(def make-config
  {:base-source-path    nil #_"notebooks" ;;files watch doesn't work
   :live-reload         true
   :format              [:quarto :html]
   :base-target-path    "docs"
   :clean-up-target-dir true
   :source-path         ["notebooks/index.clj"
                         "notebooks/steps/step_01.clj"
                         "notebooks/steps/step_02.clj"
                         "notebooks/steps/cfl_condition.clj"
                         "notebooks/steps/step_03.clj"
                         "notebooks/steps/step_04.clj"]
   :book                {:title "CFD Python in Clojure"}})

(def revealjs-live-reload-make-config
  {:base-source-path nil
   :live-reload      true
   :format           [:quarto :revealjs]
   :source-path      "clj-file-path"
   :quarto           {:format {:revealjs {:theme :serif}}}})

(comment
  (clay/config)

  (clay/stop!)

  (clay/browse!)

  ;; notebooks livereload
  (clay/make! make-config)

  ;; revealjs live reload
  (do (clay/stop!)
      (clay/make! revealjs-live-reload-make-config)))
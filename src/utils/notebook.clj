(ns utils.notebook
  (:require
   [scicloj.kindly.v4.api :as kindly]
   [scicloj.kindly.v4.kind :as kind]))

(def md (comp kindly/hide-code kind/md))

(def tex (comp kindly/hide-code kind/tex))

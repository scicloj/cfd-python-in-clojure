^:kindly/hide-code
(ns index
  (:require
   [scicloj.kindly.v4.api :as kindly]
   [scicloj.kindly.v4.kind :as kind]))

^:kindly/hide-code
(def md (comp kindly/hide-code
          kind/md))

(md
  "
# CFD Python in Clojure
")
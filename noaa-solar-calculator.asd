;;;; noaa-solar-calculator.asd

(asdf:defsystem #:noaa-solar-calculator
  :description "CL Implementation of http://www.esrl.noaa.gov/gmd/grad/solcalc/main.js"
  :author "Coseltech / NOAA.gov"
  :license "Original in Public Domain, this implementation GPLv3"
  :serial t
  :components ((:file "package")
               (:file "noaa-solar-calculator")))


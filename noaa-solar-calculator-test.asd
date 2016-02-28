(defsystem noaa-solar-calculator-test
  :author "CoselTech <coseltech@pareidolia.nl>"
  :license "GPLv3"
  :description "Tests for NOAA Solar Calculator."
  :depends-on (:noaa-solar-calculator
               :fiveam)
  :components ((:module "t"
                :serial t
                :components
                ((:file "noaa-solar-calculator")))))

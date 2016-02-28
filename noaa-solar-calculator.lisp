;;;; noaa-solar-calculator.lisp

(in-package #:noaa-solar-calculator)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defparameter *previous-read-default-float-format* *read-default-float-format*)
  (setf *read-default-float-format* 'double-float))

(defun jd (universal-time &optional timezone)
  (multiple-value-bind (seconds minute hour day month year) (decode-universal-time universal-time timezone)
    (declare (ignore seconds minute hour)) 
    (multiple-value-bind (year month) (if (<= month 2) (values (1- year) (+ month 12)) (values year month))
      (let* ((a (floor (/ year 100)))
             (b (+ (- 2 a)
                   (floor (/ a 4)))))
        (- (+ (floor (* 365.25 (+ year 4716)))
              (floor (* 30.6001 (+ month 1)))
              day
              b)
           1524.5)))))

(defun time-local (universal-time &optional timezone)
  "Result in minutes"
  (multiple-value-bind (seconds minute hour) (decode-universal-time universal-time timezone) 
    (+
     (* hour 60.0)
     minute
     (/ seconds 60.0))))

(defun tz (universal-time)
  (- (nth-value 8 (decode-universal-time universal-time))))

(defun time-julian-cent (jd)
  (/ (- jd 2451545.0) 36525.0))

(defun obliquity-correction (tee)
  "Result in degrees"
  (let ((e0 (mean-obliquity-of-ecliptic tee))
        (omega (- 125.04 (* 1934.136 tee))))
    (+ e0 (* 0.00256 (cos (deg-to-rad omega))))))

(defun mean-obliquity-of-ecliptic (tee)
  "Result in degrees"
  (let ((seconds (- 21.448 (* tee (+ 46.815 (* tee (- 5.9e-4 (* tee 0.001813))))))))
    (+ 23.0
       (/ (+ 26.0 (/ seconds 60.0))
          60.0))))

(defun geom-mean-long-sun (tee)
  "Result in degrees"
  (let* ((l0 (+ 280.46646 (* tee (+ 36000.76983 (* tee 3.032e-4))))))
    (if (plusp l0)
        (multiple-value-bind (quot rem) (floor l0 360.0)
          (declare (ignore quot))
          (if (zerop rem) 360.0 rem))
        (mod l0 360.0))))


(defun eccentricity-earth-orbit (tee)
  "Result is unitless"
  (- 0.016708634
     (* tee
        (+ 4.2037e-5
           (* 1.267e-7 tee)))))

(defun geom-mean-anomaly-sun (tee)
  "Result in degrees"
  (+ 357.52911
     (* tee
        (- 35999.05029
           (* 1.537e-4 tee)))))

(defun sun-declination (tee)
  "Result in degrees"
  (let* ((e (obliquity-correction tee))
        (lambd (sun-apparent-long tee))
        (sint (* (sin (deg-to-rad e))
                 (sin (deg-to-rad lambd)))))
    (rad-to-deg (asin sint))))

(defun sun-apparent-long (tee)
  "Result in degrees"
  (let ((o (sun-true-long tee))
        (omega (- 125.04 (* 1934.136 tee))))
    (- o 0.00569 (* 0.00478 (sin (deg-to-rad omega)))) ))

(defun sun-true-long (tee)
  "Result in degrees"
  (+ (geom-mean-long-sun tee)
     (sun-eq-of-center tee)))

(defun sun-rad-vector (tee)
  "Result in AUs"
  (let ((v (sun-true-anomaly tee))
        (e (eccentricity-earth-orbit tee))) 
    (/ (* 1.000001018
          (- 1 (* e e)))
       (+ 1 (* e
               (cos (deg-to-rad v)))))))

(defun sun-true-anomaly (tee)
  "Result in degrees"
  (+ (geom-mean-anomaly-sun tee)
     (sun-eq-of-center tee)))

(defun sun-eq-of-center (tee)
  "Result in degrees"
  (let* ((m (geom-mean-anomaly-sun tee))
         (mrad (deg-to-rad m))
         (sinm (sin mrad))
         (sin2m (sin (* 2 mrad)))
         (sin3m (sin (* 3 mrad))))
    (+ (* sinm (- 1.914602 (* tee (+ 0.004817 (* 1.4e-5 tee)))))
       (* sin2m (- 0.019993 (* 1.01e-4 tee))) (* sin3m 2.89e-4))))

(defun equation-of-time (tee)
  "Result in minutes of time"
  (let* ((epsilon (obliquity-correction tee))
         (l0 (geom-mean-long-sun tee))
         (e (eccentricity-earth-orbit tee))
         (m (geom-mean-anomaly-sun tee))
         (y (expt (tan (/ (deg-to-rad epsilon) 2.0)) 2))
         (sin2l0 (sin (* 2.0 (deg-to-rad l0))))
         (sinm (sin (deg-to-rad m)))
         (cos2l0 (cos (* 2.0 (deg-to-rad l0))))
         (sin4l0 (sin (* 4.0 (deg-to-rad l0))))
         (sin2m (sin (* 2.0 (deg-to-rad m)))))
    (* 4.0
       (rad-to-deg 
        (+ (- (* y sin2l0)
              (* 2.0 e sinm))
           (* 4.0 e y sinm cos2l0)
           (- (* 0.5 y y sin4l0))
           (- (* 1.25 e e sin2m)))))))

(defun deg-to-rad (degs)
  (/ (* pi degs)
     180.0))

(defun clamp (lo hi num)
  (min hi (max lo num)))

(defun rad-to-deg (rads)
  (/ (* 180.0 rads)
     pi))

(defun refraction-correction (exoatm-elevation)
  "Atmospheric refraction correction"
  (if (> exoatm-elevation 85.0)
      0.0
      (/ (let ((te (tan (deg-to-rad exoatm-elevation))))
           (cond
             ((> exoatm-elevation 5.0)
              (+ (- (/ 58.1 te)
                    (/ 0.07 (expt te 3)))
                 (/ 0.000086
                    (expt te 5))))
             ((> exoatm-elevation -0.575)
              (+ 1735.0
                 (* exoatm-elevation
                    (+ -518.2
                       (* exoatm-elevation
                          (+ 103.4
                             (* exoatm-elevation
                                (+ -12.79
                                   (* exoatm-elevation
                                      0.711)))))))))
             (t (/ -20.774 te))))
         3600.0)))

(defun az-el (universal-time lat lon &optional timezone)
  "Calculates multiple values, in order:  azimuth and elevation (in degrees) dark-p, exo-atmospheric angle (degrees), equation of time (minutes) and solar declination (degrees). Timezone argument such that time + timezone = UTC"
  (let* ((jday (jd universal-time timezone))
         (tl (time-local universal-time timezone))
         (tz (if timezone (- timezone) (tz universal-time)))
         (total (- (+ jday
                      (/ tl 1440.0))
                   (/ tz 24.0)))
         (tee (time-julian-cent total)))
    (%az-el tee tl lat lon tz)))

(defun %az-el (tee tl lat lon tz)
  (let* ((eq-time (equation-of-time tee))
         (theta (sun-declination tee))
         (solar-time-fix (- (+ eq-time
                               (* 4.0 lon))
                            (* 60.0 tz)))
         (true-solar-time (multiple-value-bind (quot rem) (floor (+ tl solar-time-fix) 1440)
                            (declare (ignore quot))
                            (cond 
                              ((zerop rem) 1440)
                              (t rem))))
         (hour-angle (- (/ true-solar-time 4.0) 180.0))
         (csz (clamp -1.0 1.0 (+ (* (sin (deg-to-rad lat))
                                    (sin (deg-to-rad theta)))
                                 (* (cos (deg-to-rad lat))
                                    (cos (deg-to-rad theta))
                                    (cos (deg-to-rad  hour-angle))))))
         (zenith (rad-to-deg (acos csz)))
         (az-denom (* (cos (deg-to-rad lat))
                      (sin (deg-to-rad zenith))))
         (azimuth (mod (if (> (abs az-denom) 0.001)
                           (let ((az-rad (clamp -1.0 1.0 (/ (- (* (sin (deg-to-rad lat))
                                                                  (cos (deg-to-rad zenith)))
                                                               (sin (deg-to-rad theta)))
                                                            az-denom))))
                             (* (if (plusp hour-angle)
                                    -1
                                    1)
                                (- 180.0
                                   (rad-to-deg (acos az-rad)))))
                           (if (plusp lat)
                               180.0
                               0.0))
                       360))
         (exoatm-elevation (- 90.0 zenith))
         (solar-zen (- zenith (refraction-correction exoatm-elevation)))
         (dark-p (> solar-zen 108))
         (elevation (- 90.0 solar-zen)))
    (values azimuth elevation dark-p exoatm-elevation eq-time theta)))

(defun sol-noon (universal-time lon &optional timezone)
  (let* ((jday (jd universal-time timezone))
         (tz (if timezone (- timezone) (tz universal-time)))
         (tnoon (time-julian-cent (- jday
                                     (/ lon 360.0))))
         (eq-time (equation-of-time tnoon))
         (sol-noon-offset (- 720.0
                             (* lon 4)
                             eq-time))
         (newt (time-julian-cent (+ jday
                                    (/ sol-noon-offset 1440.0))))
         (eqtime-new (equation-of-time newt))
         (sol-noon-local (+ (- 720
                               (* lon 4)
                               eqtime-new)
                            (* tz 60.0))))
    (mod sol-noon-local 1440.0)))

(define-condition no-sun-rise-set-today ()())

(defun hour-angle-sunrise (lat solar-dec)
  "Result in radians"
  (let* ((lat-rad (deg-to-rad lat))
         (sd-rad (deg-to-rad solar-dec))
         (ha
          (acos
           (- (/ (cos (deg-to-rad 90.833))
                 (* (cos lat-rad) (cos sd-rad)))
              (* (tan lat-rad) (tan sd-rad))))))
    (when (complexp ha)
      (error 'no-sun-rise-set-today))
    ha))

(defun leap-year-p (year)
  (or (and (= (mod year 4) 0) (not (= (mod year 100) 0))) (= (mod year 400) 0)))
 
(defun doy-from-jd (jd)
  (let* ((z (floor (+ jd 0.5)))
         (f (- (+ jd 0.5) z))
         (a (if (< z 2299161)
                z
                (let ((alpha (floor (/ (- z 1867216.25) 36524.25))))
                  (+ z 1 alpha (- (floor (/ alpha 4)))))))
         (b (+ a 1524))
         (c (floor (/ (- b 122.1) 365.25)))
         (d (floor (* 365.25 c)))
         (e (floor (/ (- b d) 30.6001)))
         (day (+ (- b d (floor (* 30.6001 e))) f))
         (month (if (< e 14) (- e 1) (- e 13)))
         (year (if (> month 2) (- c 4716) (- c 4715)))
         (k (if (leap-year-p year) 1 2)))
    (+ (- (floor (/ (* 275 month) 9)) (* k (floor (/ (+ month 9) 12))))
       day (- 30))))

(defun rise-set-utc (mode jday lat lon)
               (let* ((tee (time-julian-cent jday))
                      (eq-time (equation-of-time tee))
                      (solar-dec (sun-declination tee))
                      (hour-angle (hour-angle-sunrise lat solar-dec))
                      (delta (+ lon
                                (rad-to-deg (ecase mode
                                              (:rise hour-angle)
                                              (:set (- hour-angle)))))))
                 (- 720 (* 4.0 delta) eq-time)))

(defun jd-of-next-prev-rise-set (direction mode jd lat lon tz)
  (let* ((time
          (loop for time = (handler-case
                               (rise-set-utc mode jd lat lon)
                             (no-sun-rise-set-today () nil))
             while (null time)
             do (ecase direction
                  (:prev (decf jd))
                  (:next (incf jd)))
             finally (return time)))
         (time-local (+ time (* tz 60.0))))
    (multiple-value-bind (quot rem) (floor time-local 1440.0)
      (values (+ jd quot) rem))))

(defun todays-sun-rise-set (mode jday lat lon tz)
  (let* ((time-utc (rise-set-utc mode jday lat lon))
         (new-time-utc (rise-set-utc mode (+ jday
                                             (/ time-utc
                                                1440.0)) lat lon))
         (time-local (+ new-time-utc
                        (* tz 60.0)))
         (tee (time-julian-cent (+ jday
                                   (/ new-time-utc
                                      1440.0))))
         (az (nth-value 0 (%az-el tee time-local lat lon tz))))
    
    (values az time-local)))

(defun sun-rise-set (universal-time lat lon &optional timezone)
  (let* ((jday (jd universal-time timezone))
         (tz (if timezone (- timezone) (tz universal-time))))
    (flet ((rise-set (mode)
             (multiple-value-bind (s m h day month year) (decode-universal-time universal-time (- tz))
               (declare (ignore s m h))
               (handler-case
                   (multiple-value-bind (az time-local) (todays-sun-rise-set mode jday lat lon tz)
                     (list :az az
                           :time (round (+ (encode-universal-time 0 0 0 day month year (- tz))
                                           (* 60 time-local)))
                           :today-p t))
                 (no-sun-rise-set-today ()
                   (let ((doy (doy-from-jd jday)))
                     (flet ((find-rise-set (direction)
                              (let* ((other-jday (jd-of-next-prev-rise-set direction mode jday lat lon tz))
                                     (day-offset-secs (* 24 60 60 (- other-jday jday))))
                                (multiple-value-bind (az time-local) (todays-sun-rise-set mode other-jday lat lon tz)
                                  (list :az az
                                        :time (round (+ (encode-universal-time 0 0 0 day month year (- tz))
                                                        (* 60 time-local)
                                                        day-offset-secs))
                                        :today-p nil)))))
                       (if (or (and (> lat 66.4) (> doy 79) (< doy 267))
                               (and (< lat (- 66.4)) (or (< doy 83) (> doy 263))))
                           (find-rise-set (ecase mode (:rise :prev) (:set :next)))
                           (find-rise-set (ecase mode (:rise :next) (:set :prev)))))))))))
      (list :rise (rise-set :rise) :set (rise-set :set)))))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (setf *read-default-float-format* *previous-read-default-float-format*))

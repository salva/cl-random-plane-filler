(ql:quickload "lisp-magick")
(ql:quickload "matlisp")

(defpackage #:circle-pixels
  (:use #:cl :matlisp :alexandria))

(in-package circle-pixels)

(proclaim '(optimize speed))

(defparameter *dilution* 10)

;(defun draw-pixel (dw x0 y0 pixel-size r g b pw-border pw-r pw-g pw-b)
;  (let* ((radius (* 0.5 pixel-size))
;         (ox (+ x0 radius))
;         (oy (+ y0 radius))
;         (area-total (* pi radius radius))
;         (area-r (/ (* *dilution* area-total r) 256))
;         (area-g (/ (* *dilution* area-total g) 256))
;         (area-b (/ (* *dilution* area-total b) 256))
;         (r-b (sqrt area-b))
;         (r-g (sqrt (+ area-b area-g)))
;         (r-r (sqrt (+ area-b area-g area-r))))
;    (magick:draw-set-stroke-color dw pw-border)
;    (magick:draw-set-fill-color dw pw-r)
;    (magick:draw-circle dw ox oy (+ ox r-r) oy)
;    (magick:draw-set-fill-color dw pw-g)
;    (magick:draw-circle dw ox oy (+ ox r-g) oy)
;    (magick:draw-set-fill-color dw pw-b)
;    (magick:draw-circle dw ox oy (+ ox r-b) oy)))

(defparameter *rgbk-wheel-list* '((1 0 0) (0 1 0) (0 0 1) (0 0 0) (1 1 1)))
(defparameter *cmyk-wheel-list* '((0 1 1) (1 0 1) (1 1 0) (0 0 0) (1 1 1)))

(defun color-to-matrix-b (color)
  (let ((b (matlisp::zeros '(4 1))))
    (setf (matlisp::ref b 3 0) 1.0)
    (do ((comp color (cdr comp))
         (j 0 (+ j 1)))
        ((null comp) b)
      (setf (matlisp::ref b j 0) (coerce (car comp) 'double-float)))))

(defun wheel-to-matrix-A (wheel)
  (let* ((n-colors (length wheel))
         (A (matlisp::zeros (list 4 n-colors))))
    (do ((color wheel (cdr color))
         (i 0 (+ i 1)))
        ((null color) A)
      (do ((comp (car color) (cdr comp))
           (j 0 (+ j 1)))
          ((null comp))
        ;;(print (list :i i :j j :comp (car comp)))
        (setf (matlisp::ref A j i) (coerce (car comp) 'double-float)))
      (setf (matlisp::ref A 3 i) 1.0))))

(defun decompose-color-with-matrices (A b)
  (let* ((sol (matlisp::gelsy A b))
         (dims (matlisp::dims A))
         (n-colors (cadr dims))
         (negs))
    ;;(print (list :A A :b b '--> :sol sol))
    (dotimes (j n-colors)
      (when (< (matlisp::ref sol j 0) -0.01)
        (setf negs (cons j negs))))
    (if (null negs) sol
        (let* ((rows (car dims))
               (extra (length negs))
               (A+ (matlisp::zeros (list (+ rows extra) n-colors)))
               (b+ (matlisp::zeros (list (+ rows extra) 1))))
          (dotimes (j rows)
            (setf (matlisp::ref b+ j 0) (matlisp::ref b j 0))
            (dotimes (i n-colors)
              (setf (matlisp::ref A+ j i) (matlisp::ref A j i))))
          (do ((zero negs (cdr zero))
               (j rows (+ j 1)))
              ((null zero))
            (setf (matlisp::ref A+ j (car zero)) 1.0)
            (setf (matlisp::ref b+ j 0) 0.0)
            (dotimes (j 3)
              (setf (matlisp::ref A+ j (car zero)) 0.0)))
          (decompose-color-with-matrices A+ b+)))))

(defun decompose-color (color wheel)
  (let* ((A (wheel-to-matrix-A wheel))
         (b (color-to-matrix-b color))
         (sol (decompose-color-with-matrices A b))
         (n-colors (car (matlisp::dims sol)))
         (result))
    ;; (print (list :n-colors n-colors :sol sol))
    (do ((j (- n-colors 1) (- j 1)))
        ((< j 0) result)
      (setf result (cons (matlisp::ref sol j 0) result)))))
    
(defun draw-pixel-with-wheel (dw x0 y0 pixel-size r g b wheel &rest pws)
  (let* ((radius (* 0.5 pixel-size))
         (area-total (* pi radius radius))
         (inv-area-total (/ 1.0 area-total))
         (comps (decompose-color (list (* r inv-area-total)
                                       (* g inv-area-total)
                                       (* b inv-area-total))
                                 wheel)))
    ;(print (list :comps comps))
    (let* ((ox (+ x0 radius))
           (oy (+ y0 radius))
           (comps-pairs (sort (mapcar #'cons comps pws)
                              #'(lambda (a b) (<= (car a) (car b)))))
           (areas-pairs)
           (acu-area 0))
      (dolist (pair comps-pairs)
        (let ((area (* area-total (car pair))))
          (setf acu-area (+ acu-area (if (< area 0) 0 area)))
          (setf areas-pairs (cons (cons acu-area (cdr pair)) areas-pairs))))
      (dolist (pair areas-pairs)
        (let* ((area (car pair))
               (radius (sqrt area))
               (pw (cdr pair)))
          (magick:draw-set-fill-color dw pw)
          (magick:draw-circle dw ox oy (+ ox radius) oy))))))

(defun draw-pixel (dw x0 y0 pixel-size r g b pw-border pw-r pw-g pw-b)
  (print '- )
  (magick:draw-set-stroke-color dw pw-border)
  (draw-pixel-with-wheel dw
                         x0 y0
                         pixel-size
                         (* *dilution* r)
                         (* *dilution* g)
                         (* *dilution* b)
                         *rgbk-wheel-list*
                         pw-r pw-g pw-b))

(defun circletize (file-in file-out &key (pixel-size 32))
  (print "loading file")(finish-output)
  (magick:with-magick-wand (wand-in :load file-in)
    (print "file loaded")(finish-output)
    (let* ((width  (magick:get-image-width wand-in))
           (height (magick:get-image-height wand-in))
           (width-pixels  (ceiling width pixel-size))
           (height-pixels (ceiling height pixel-size)))
      ;;(print (list :width width :height height :width-pixels width-pixels :height-pixels height-pixels))
      (magick:with-drawing-wand (dw-out)
        (magick:with-pixel-wand (pw-black :comp (0 0 0))
          (magick:with-pixel-wand (pw-r :comp (255 0 0))
            (magick:with-pixel-wand (pw-g :comp (0 255 0))
              (magick:with-pixel-wand (pw-b :comp (0 0 255))
                (magick:with-pixel-data (pd-in wand-in)
                  (dotimes (x-p width-pixels)
                    ;;(print (list :x-p x-p))
                    (dotimes (y-p height-pixels)
                      (let* ((x0 (floor (* x-p width) width-pixels))
                             (y0 (floor (* y-p height) height-pixels))
                             (x1 (floor (* (+ x-p 1) width) width-pixels))
                             (y1 (floor (* (+ y-p 1) height) height-pixels))
                             (dx (- x1 x0))
                             (dy (- y1 y0))
                             (source-pixel-area (* dx dy))
                             (acu-r 0)
                             (acu-g 0)
                             (acu-b 0))
                        (dotimes (x dx)
                          (dotimes (y dy)
                            (multiple-value-bind (r b g) (magick:get-pixel pd-in (+ x0 x) (+ y0 y))
                              (setf acu-r (+ acu-r r))
                              (setf acu-g (+ acu-g g))
                              (setf acu-b (+ acu-b b)))))
                        (let ((mean-r (floor acu-r source-pixel-area))
                              (mean-g (floor acu-g source-pixel-area))
                              (mean-b (floor acu-b source-pixel-area))
                              (x0 (* x-p pixel-size))
                              (y0 (* y-p pixel-size)))
                          (draw-pixel dw-out x0 y0 pixel-size mean-r mean-g mean-b
                                      pw-black pw-r pw-g pw-b))))))))))
        (magick:with-magick-wand (wand-out :create (* pixel-size width-pixels) (* pixel-size height-pixels)
                                           :comp (0 0 0))
          (print "drawing image")
          (magick:draw-image wand-out dw-out)
          (print "saving image")
          (magick:write-image wand-out file-out))))))

(ql:quickload "alexandria")
(ql:quickload "vecto")
(defpackage #:random-filler
  (:use #:cl :alexandria))

(in-package random-filler)

(proclaim '(optimize (debug 3)))

(defconstant max-flat-filler-length 10)

(defun vectors-dim (&rest vs)
  (if vs
      (let ((s (length(car vs))))
        (dolist (v (cdr vs) s)
          (unless (eql s (length v))
            (error "vectors dimensions do not match"))))
      NIL))

(defun v-norm2 (v)
  (let ((dim (vectors-dim v))
        (acu 0))
    (dotimes (i dim acu)
      (let ((x (aref v i)))
        (setf acu (+ acu (* x x)))))))

(defun v-norm (v)
  (sqrt (v-norm2 v)))

(defun v. (v0 v1)
  (declare ((array number *) v0 v1))
  (let ((dim (vectors-dim v0 v1))
        (acu 0))
    (dotimes (i dim acu)
      (setf acu (+ acu (* (aref v0 i) (aref v1 i)))))))

(defun v+ (v &rest ws)
  (let ((dim (vectors-dim v))
        (out (copy-array v)))
    (dolist (w ws out)
      (dotimes (i dim)
        (setf (aref out i) (+ (aref out i) (aref w i)))))))

(defun v- (v &rest ws)
  (let ((dim (vectors-dim v))
        (out (copy-array v)))
    (dolist (w ws out)
      (dotimes (i dim)
        (setf (aref out i) (- (aref out i) (aref w i)))))))

(defun -v (v)
  (let* ((dim (vectors-dim v))
         (out (make-array dim)))
    (dotimes (i dim out)
      (setf (aref out i) (- (aref v i))))))

(defun v*s (v s)
  (let* ((dim (vectors-dim v))
         (out (make-array dim)))
    (dotimes (i dim out)
      (setf (aref out i) (* s (aref v i))))))

(defun v-largest-component-ix (v)
  (let ((dim (vectors-dim v))
        (largest NIL)
        (max-x -1))
    (dotimes (i dim largest)
      (let ((x (abs (aref v i))))
        (if (< max-x x)
            (setf largest i
                  max-x x))))))

(defun v-copy-and-set (v ix x)
  (let ((out (copy-array v)))
    (setf (aref out ix) x)
    out))

(defun v-hypervol (v)
  (let ((dim (vectors-dim v))
        (acu 1.0))
    (dotimes (i dim acu)
      (setf acu (* acu (aref v i))))))

(defun v-box-list (v ws)
  (let ((dim (vectors-dim v))
        (min (copy-array v))
        (max (copy-array v)))
    (dolist (w ws (cons min max))
      (dotimes (i dim)
        (let ((x (aref w i)))
          (if (> (aref min i) x)
              (setf (aref min i) x)
              (if (< (aref max i) x)
                  (setf (aref max i) x))))))))

(defun v-box (v &rest ws)
  (v-box-list v ws))

(defun v-nearest-in-box (v w0 &rest ws)
  (let* ((dim (vectors-dim v))
         (box (v-box-list w0 ws))
         (b0 (car box))
         (b1 (cdr box))
         (out (copy-array v)))
    (dotimes (i dim out)
      (let ((x (aref out i))
            (min-x (aref b0 i))
            (max-x (aref b1 i)))
        (if (< x min-x)
            (setf (aref out i) min-x)
            (if (> x max-x)
                (setf (aref out i) max-x)))))))

(defun v-farthest-in-box (v w0 &rest ws)
  (let* ((dim (vectors-dim v))
         (box (v-box-list w0 ws))
         (min (car box))
         (max (cdr box))
         (out (copy-array v)))
    (dotimes (i dim out)
      (let ((x (aref out i))
            (min-x (aref min i))
            (max-x (aref max i)))
        (setf (aref out i)
              (if (> (- x min-x) (- max-x x)) min-x max-x))))))

(defun v-dist2 (v0 v1)
  (let ((dim (vectors-dim v0 v1))
        (acu 0))
    (dotimes (i dim acu)
      (let ((x (- (aref v0 i) (aref v1 i))))
        (setf acu (* x x))))))

(defun v-dist (v0 v1)
  (sqrt (v-dist2 v0 v1)))

(defclass shape () ())
(defgeneric shape-area (shape))
(defgeneric shape-split (shape))
(defgeneric shapes-intersection-area (shape0 shape1))
(defgeneric shapes-touch-p (shape0 shape1))
(defgeneric shape-touches-point-p (shape point))
(defgeneric shapes-dist (shape0 shape1))
(defgeneric shape-max-dist2-to-point (shape point))
(defgeneric shape-dist2-to-point (shape point))
(defgeneric shape-random-point (shape))

(defclass circle (shape)
  ((o :type (array real 2) :initarg :o :initform #2(0.0 0.0) :reader circle-o)
   (r :type (real) :initarg :r :initform 1.0 :reader circle-r)))

(defmethod shape-area ((c circle))
  (with-slots (r) c
    (* pi r r)))

(defclass rectangle (shape)
  ((o0 :type (array real 2) :initarg :o0 :initform #2(0.0 0.0) :reader rectangle-o0)
   (o1 :type (array real 2) :initarg :o1 :initform #2(1.0 1.0) :reader rectangle-o1)))

(defmethod shape-area ((r rectangle))
  (with-slots (o0 o1) r
    (abs (v-hypervol (v- o0 o1)))))

(defmethod shape-max-dist2-to-point ((r rectangle) p)
  (with-slots (o0 o1) r
    (v-dist2 p (v-farthest-in-box o0 01))))

(defmethod shape-dist2-to-point ((r rectangle) p)
  (with-slots (o0 o1) r
    (v-dist2 p (v-nearest-in-box o0 o1))))

(defmethod shape-random-point ((r rectangle))
  (with-slots (o0 o1) r
    (let* ((dim (length o0))
           (out (make-array dim)))
      (dotimes (i dim out)
        (let ((r (random 1.0)))
          (setf (aref out i) (+ (* r (aref o0 i)) (* (- 1 r) (aref o1 i)))))))))

(defmethod shape-touches-point-p ((c circle) p)
  (with-slots (o r) c
    (<= (v-dist2 c p) (* r r))))

(defmethod shape-split ((r rectangle))
  (with-slots (o0 o1) r
    (let* ((axis (v-largest-component-ix (v- o0 o1)))
           (x (* 0.5 (+ (aref o0 axis) (aref o1 axis))))
           (o01 (v-copy-and-set o1 axis x))
           (o10 (v-copy-and-set o0 axis x)))
      (list (make-instance 'rectangle :o0 o0 :o1 o01)
            (make-instance 'rectangle :o0 o10 :o1 o1)))))

(defmethod shapes-dist ((cir circle) (rect rectangle))
  (with-slots (o r) cir
    (with-slots (o0 o1) rect
      (let* ((nearest (v-nearest-in-box o o0 o1))
             (d (- (v-dist o nearest) r)))
        (if (< d 0) 0 d)))))

(defmethod shapes-touch-p (shape0 shape1)
  (eql (shapes-dist shape0 shape1) 0))

(defun v-angle (v0 v1)
  (let ((dim (vectors-dim v0 v1)))
    (if (= dim 2)
        (let ((dot   (+ (* (aref v0 0) (aref v1 0))
                        (* (aref v0 1) (aref v1 1))))
              (cross (- (* (aref v0 0) (aref v1 1))
                        (* (aref v0 1) (aref v1 0)))))
          (atan cross dot))
        (let ((a0 (v-norm v0)))
          (if (= a0 0) 0
              (let* ((u0 (v*s v0 (/ 1.0 a0)))
                     (p (v. v1 u0)))
                (atan (v-norm (v- v1 (v*s u0 p))) p)))))))

(defun ia2-circle-segment (r2 a b)
  (let* ((ab (v- b a))
         (c2 (v-norm2 ab)))
    (if (<= c2 0) 0
        (let* ((c1 (v. a ab))
               (c0 (- (v-norm2 a) r2))
               (discriminant (- (* c1 c1) (* c0 c2))))
          (print (list :c1 c1 :c0 c0 :c2 c2 :discriminant discriminant :a a :ab ab))
          (if (> discriminant 0)
              (let* ((inv-c2 (/ 1.0 c2))
                     (sqrt-discriminant (sqrt discriminant))
                     (alfa0 (* inv-c2 (- (- c1) sqrt-discriminant)))
                     (alfa1 (* inv-c2 (+ (- c1) sqrt-discriminant))))
                (print (list :inv-c2 inv-c2 :sqrt-discriminant sqrt-discriminant :alfa0 alfa0 :alfa1 alfa1))
                (if (and (< alfa0 1.0) (> alfa1 0))
                    (let ((beta 0)
                          (p0)
                          (p1))
                      (if (> alfa0 0)
                          (setf p0 (v+ a (v*s ab alfa0))
                                beta (+ beta (v-angle a p0)))
                          (setf p0 a))
                      (if (<= alfa1 1)
                          (setf p1 (v+ a (v*s ab alfa1))
                                beta (+ beta (v-angle p1 b)))
                          (setf p1 b))
                      (return-from ia2-circle-segment
                        (+ (- (* (aref p0 0) (aref p1 1))
                              (* (aref p0 1) (aref p1 0)))
                           (* beta r2)))))))
          (* r2 (v-angle a b))))))

(defmethod shapes-intersection-area ((cir circle) (rect rectangle))
  (with-slots (o r) cir
    (with-slots (o0 o1) rect
      (let ((a (v- o0 o))
            (c (v- o1 o)))
        (if (or (<= (aref c 0) (- r))
                (>= (aref a 0) r)
                (<= (aref c 1) (- r))
                (>= (aref a 1) r))
            0
            (let ((b (vector (aref a 0) (aref c 1)))
                  (d (vector (aref c 0) (aref a 1)))
                  (r2 (* r r)))
              (abs (* 0.5 (+ (ia2-circle-segment r2 a b)
                             (ia2-circle-segment r2 b c)
                             (ia2-circle-segment r2 c d)
                             (ia2-circle-segment r2 d a))))))))))

(defclass filler ()
  ((box :initarg :box
        :initform (make-instance 'rectangle)
        :accessor filler-box)
   (free-area :accessor filler-free-area)))

(defgeneric filler-insert (filler shape))
(defgeneric filler-insert-unchecked (filler shape))
(defgeneric filler-shapes (filler))
(defgeneric filler-subfillers (filler)) 

(defgeneric filler-total-area (filler))

(defmethod filler-total-area (f)
  (shape-area (filler-box f)))

(defmethod initialize-instance :after ((f filler) &rest initargs)
  (declare (ignore initargs))
  (setf (filler-free-area f) (shape-area (filler-box f))))

(defclass filler-flat (filler)
  ((shapes :initform NIL :reader filler-shapes)))

(defmethod filler-subfillers ((f filler-flat)) NIL)

(defclass filler-tree (filler)
  ((subfillers :initarg :subfillers :reader filler-subfillers)))

(defgeneric filler-add-shapes-to-hashtable (filler hashtable))

(defmethod filler-shapes ((f filler-tree))
  (let ((ht (make-hash-table)))
    (filler-add-shapes-to-hashtable f ht)
    (hash-table-keys ht)))

(defmethod filler-add-shapes-to-hashtable ((f filler-tree) ht)
  (dolist (subfiller (filler-subfillers f))
    (filler-add-shapes-to-hashtable subfiller ht)))

(defmethod filler-add-shapes-to-hashtable ((f filler-flat) ht)
  (dolist (shape (filler-shapes f))
    (setf (gethash shape ht) t)))

(defmethod filler-insert ((filler filler) shape)
  (if (shapes-touch-p shape (filler-box filler))
      (filler-insert-unchecked filler shape)))

(defmethod filler-insert-unchecked ((filler filler-flat) shape)
  (with-slots (shapes free-area box) filler
    (setf shapes (cons shape shapes))
    (setf free-area (- free-area
                       (shapes-intersection-area shape box)))
    (if (< ( * 2 free-area) (shape-area box))
        (let ((old-shapes shapes)
              (subfillers (mapcar #'(lambda (box) (make-instance 'filler-flat :box box))
                                  (shape-split box))))
          (change-class filler 'ramdom-plane-filler-tree :subfillers subfillers)
          (dolist (shape old-shapes)
            (filler-insert-unchecked filler shape))))))
  
(defmethod filler-insert-unchecked ((filler filler-tree) shape)
  (let ((subfillers (filler-subfillers filler)))
    (dolist (subfiller subfillers)
      (filler-insert subfiller shape))
    (setf (filler-free-area filler)
          (reduce #'+ (mapcar #'filler-free-area subfillers)))))

(defgeneric filler-add-touching-shapes-to-hashtable (filler shape ht))
(defgeneric filler-add-touching-shapes-to-hashtable-unchecked (filler shape ht))

(defgeneric filler-touching-shapes (f s))
(defgeneric filler-nearest-shape-to-point (f p))

(defmethod filler-touching-shapes ((f filler) (s shape))
  (let ((ht (make-hash-table)))
    (filler-add-touching-shapes-to-hashtable f s ht)
    (hash-table-keys ht)))

(defmethod filler-add-touching-shapes-to-hashtable ((f filler) (s shape) ht)
  (if (shapes-touch-p s (filler-box f))
      (filler-add-touching-shapes-to-hashtable-unchecked f s ht)))

(defmethod filler-add-touching-shapes-to-hashtable-unchecked ((f filler-flat) (s shape) ht)
  (dolist (s1 (filler-shapes f))
    (if (shapes-touch-p s1 s)
        (setf (gethash s1 ht) t))))

(defmethod filler-add-touching-shapes-to-hashtable-unchecked ((f filler-tree) (s shape) ht)
  (dolist (subfiller (filler-subfillers f))
    (filler-add-touching-shapes-to-hashtable subfiller s ht)))

(defun insert-into-queue (queue v d2)
  (if queue
      (if (<= d2 (cdar queue))
          (cons (cons v d2) queue)
          (cons (car queue) (insert-into-queue (cdr queue) v d2)))
      (list (cons v d2))))

(defmethod filler-nearest-shape-to-point ((f filler) p)
  (let* ((checked-shapes (make-hash-table))
         (best-d2 (shape-max-dist2-to-point (filler-box f) p))
         (best-shape)
         (queue))
    (do* ((head NIL (pop queue))
          (f f (car head))
          (head-d2 0 (cdr head)))
         ((and head-d2 (< head-d2 best-d2)))
      (let ((subfillers (filler-subfillers f)))
        (if subfillers
            (dolist (sf subfillers)
              (let ((d2 (shape-dist2-to-point (filler-box sf) p)))
                (if (< d2 best-d2)
                    (setf queue (insert-into-queue queue sf d2)))))
            (dolist (s (filler-shapes f))
              (unless (gethash s checked-shapes)
                (let ((d2 (shape-dist2-to-point s p)))
                  (if (< d2 best-d2)
                      (setf best-d2 d2
                            best-shape s))
                  (setf (gethash s checked-shapes) t)))))))
    (values best-shape best-d2)))

(defgeneric filler-random-free-point (f))
(defgeneric filler-random-free-point-unchecked (f))

(defmethod filler-random-free-point ((f filler))
  (if (< 0 (filler-free-area f))
      (filler-random-free-point-unchecked f)))

(defgeneric filler-random-free-point-unchecked (f))

(defmethod filler-random-free-point-unchecked ((f filler-flat))
  (with-slots (box shapes) f
    (let ((p (shape-random-point box)))
      (dolist (s shapes p)
        (if (shape-touches-point-p s p)
            (return (filler-random-free-point-unchecked f)))))))
  
(defmethod filler-random-free-point-unchecked ((f filler-tree))
  (with-slots (subfillers free-area) f
    (let ((ix (random free-area)))
      (dolist (sf subfillers NIL)
        (with-slots (free-area) sf
          (if (> ix free-area)
              (setf ix (- ix free-area))
              (return (filler-random-free-point sf))))))))

(defun draw-circles (file-name circles)
  (let* ((scl 1024)
         (px0 1024)
         (py0 1024))
    (defun scl-x (x) (+ (* x scl) px0))
    (defun scl-y (y) (+ (* y scl) py0))
    (vecto:with-canvas (:width (* 2 scl) :height (* 2 scl))
      (vecto:set-rgb-stroke 1 0 0)
      (vecto:set-line-width 3)
      (vecto:set-rgb-fill 1 1 0)
      (dolist (c circles)
        (with-slots (o r) c
          (print (list :ox (aref o 0) :pox (scl-x (aref o 0))
                       :oy (aref o 1) :poy (scl-y (aref o 1))
                       :r r))
          (vecto:centered-circle-path (scl-x (aref o 0)) (scl-y (aref o 1)) (* scl r))
          (vecto:stroke)))
      (vecto:save-png file-name))))

(defun random-circle ()
  (let ((r (random 1.0))
        (o (vector (- (random 2.0) 1.0)
                   (- (random 2.0) 1.0))))
    (make-instance 'circle :r r :o o)))

(defun random-circles (n)
  (let ((circles))
    (dotimes (i n circles)
      (setf circles (cons (random-circle) circles)))))

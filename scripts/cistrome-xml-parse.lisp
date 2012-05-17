#!/usr/bin/sbcl --script
(require 'xmls)

(with-open-file (file "cistrome.xml")
  (defparameter *parsed*  (xmls:parse file)))
(defvar *a* nil)
(defvar *t* nil)
(defvar *c* nil)
(defvar *g* nil)

(defun prob->count (prob)
  (round (* 100 prob)))

(defmacro push-bas (bas func)
  `(push
    (funcall ,func (read-from-string
			    (third ,bas)))
    (symbol-value (intern (format nil "*~a*" (first ,bas))))))

(defmacro clear-bas (bas-symbol)
  `(setf (symbol-value (intern (format nil "*~a*" ,bas-symbol))) nil))

(defmacro format-bas (bas-symbol)
  `(progn
     (dolist (pos (symbol-value (intern (format nil "*~a*" ,bas-symbol))))
       (format t "~a~3t" pos))
     (format t "~%")))

(defun clear-all ()
  (progn
    (clear-bas 'a)
    (clear-bas 't)
    (clear-bas 'c)
    (clear-bas 'g)))

(defun push-all()
  (loop for i in (cddr *parsed*) do
       (format t "~%~a~%" (second (car (second i))))
       (loop for mid in (cdr i) do
	    (clear-all)
	    (if (equal "pssm" (car mid ))
		(progn
		  (dolist (pos (cddr mid ))
		    (dolist (bas (cddr pos))
		      ;; (push-bas bas #'prob->count)
		      (push-bas bas #'(lambda (x) x))
		      ))
		  (format-bas 'a)
		  (format-bas 't)
		  (format-bas 'g)
		  (format-bas 'c))))))
(push-all)

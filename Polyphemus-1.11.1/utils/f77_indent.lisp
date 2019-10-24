(defun f77_indent ()
  (fortran-mode)
  (indent-region (point-min) (point-max) ())
  (save-buffer))

(defun f77_indent_all ()
  (save-excursion
    (interactive)
    (dolist (buffer (buffer-list))
      (when (not  (string= "*scratch*" (buffer-name)))
	(set-buffer buffer)
	(fortran-mode)
	(indent-region (point-min) (point-max) nil)
	(save-buffer)))))

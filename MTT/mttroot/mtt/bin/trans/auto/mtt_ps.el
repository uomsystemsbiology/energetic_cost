(TeX-add-style-hook "mtt_ps"
 (lambda ()
    (LaTeX-add-labels
     "fig:#2")
    (TeX-add-symbols
     '("tfig" 4)
     '("fig" 4))
    (TeX-run-style-hooks
     "hyperref"
     "color"
     "graphicx")))


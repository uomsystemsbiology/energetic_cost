(TeX-add-style-hook
 "Figures"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "a4paper" "12pt")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art12"
    "a4wide"
    "graphicx"
    "subfigure"
    "fancybox"
    "color")
   (TeX-add-symbols
    '("SubFig" 3)
    '("Fig" 2))
   (LaTeX-add-labels
    "subfig:#1"
    "fig:comparison"
    "fig:pGate_HH"
    "fig:membrane_resp"
    "fig:membrane_resp_energy")))


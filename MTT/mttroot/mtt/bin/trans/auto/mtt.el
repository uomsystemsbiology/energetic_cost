(TeX-add-style-hook "mtt"
 (lambda ()
    (TeX-add-symbols
     '("BB" 2)
     '("B" 1)
     '("Ref" 1)
     '("Input" 1)
     '("txt" 2))
    (TeX-run-style-hooks
     "fancyheadings"
     "amsmath"
     "supertabular"
     "moreverb"
     "verbatim"
     "mathptm"
     "times")))


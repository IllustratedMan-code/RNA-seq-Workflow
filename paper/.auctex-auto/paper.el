(TeX-add-style-hook
 "paper"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "11pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8") ("fontenc" "T1") ("ulem" "normalem")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art11"
    "inputenc"
    "fontenc"
    "graphicx"
    "grffile"
    "longtable"
    "wrapfig"
    "rotating"
    "ulem"
    "amsmath"
    "textcomp"
    "amssymb"
    "capt-of"
    "hyperref")
   (LaTeX-add-labels
    "sec:org673c52e"
    "sec:orga234600"
    "sec:org461de2b"
    "sec:orgd2cd0a0"
    "sec:orgf0e41a5"
    "sec:org71dd47e"
    "sec:org8b88d21"
    "sec:org1492adf"
    "sec:org67be582"
    "sec:org7926636"
    "sec:org70b0d3f"
    "sec:orgd728dbd"
    "sec:orgf8494aa"
    "sec:org1bf6993"
    "sec:org60543b7"
    "sec:orgce86978"
    "sec:org18f899d"
    "sec:org001c7de"
    "sec:org5268e66"
    "sec:orgce98a64"
    "sec:orgce0e1d8"
    "sec:org37ba1d8"
    "sec:orgac833e6"
    "sec:orga86a084"
    "sec:orgd785b30"
    "sec:org75379a0"
    "sec:org012262b"
    "sec:org4aa5cb9"
    "sec:orgac3c73d"
    "sec:orgfd6a939"
    "sec:org3d6e145"
    "sec:orgf75c4fd"
    "sec:org740262e"
    "sec:org42e7e0f"
    "sec:orge600387"
    "sec:orgcde0dbf"
    "sec:org87a18f0"
    "sec:orga0a7bff"
    "sec:org6a2089e"
    "sec:orgd8ac333"))
 :latex)


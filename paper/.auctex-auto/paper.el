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
    "sec:org98627f3"
    "sec:orgf8a70d2"
    "sec:org12b1abf"
    "sec:orge9bf373"
    "Overview"
    "sec:org9e72812"
    "sec:org3397a4f"
    "sec:org0898722"
    "sec:org4512d2f"
    "sec:orgd318d42"
    "sec:org152e034"
    "sec:org87b4262"
    "sec:org5fc1448"
    "sec:org1ce2aea"
    "sec:orgeb4034f"
    "sec:orgfb66112"
    "sec:orge2c038c"
    "sec:org12523fb"
    "sec:org33c75e2"
    "sec:org746a210"
    "sec:orga1cd89e"
    "sec:org44fcae3"
    "sec:orgc90c3df"
    "sec:org64cbd17"
    "sec:orga4c0677"
    "sec:org337ecb7"
    "sec:orgfbf9d05"
    "sec:orge4eea0a"
    "sec:orgff8dcba"
    "sec:orgf313ebd"
    "sec:org015ac2f"
    "sec:orgaa0a633"
    "sec:orga9bd67a"
    "sec:org311c7e2"
    "sec:org47dc028"
    "sec:orgc85c09f"
    "sec:orge764dd9"
    "sec:org8266250"
    "sec:org1f076cc"
    "sec:orgd918233"))
 :latex)


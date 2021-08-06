(TeX-add-style-hook
 "paper"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "11pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8") ("fontenc" "T1") ("ulem" "normalem")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
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
    "sec:org52343fb"
    "sec:orgb459187"
    "sec:org2e03ff4"
    "sec:orgc609234"
    "sec:org753876c"
    "sec:orge9a89da"
    "sec:orgbb507c5"
    "sec:org50ae625"
    "sec:org4a972d3"
    "sec:org650badf"
    "sec:org4ff4065"
    "sec:org6f660e8"
    "sec:orgb19913f"
    "sec:orgb667e92"
    "sec:org526f677"
    "sec:orgce9bf6f"
    "sec:orge417252"
    "sec:org61576e7"
    "sec:orga7443f6"
    "sec:org80b8500"
    "sec:org2f435e9"
    "sec:orgb78b77b"
    "sec:org305e391"
    "sec:orgf38097f"
    "sec:orgf9dbc02"
    "sec:orgdd49d6b"
    "sec:org655fa3a"
    "sec:org9a76e36"
    "sec:org87ea6df"
    "sec:org81930ab"
    "sec:orgb02379f"
    "sec:orgcce9fb5"
    "sec:org7508fe2"
    "sec:orgfaa377a"
    "sec:org1667a82"
    "sec:org9f24c09"))
 :latex)


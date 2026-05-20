;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "REFERENCES"
 (lambda ()
   (setq TeX-command-extra-options
         "-synctex=1")
   (LaTeX-add-bibitems
    "Buerkner2019"
    "Engle1995"
    "Engle2001a"
    "Engle2002"
    "Rast2020"
    "Yao2018"))
 '(or :bibtex :latex))


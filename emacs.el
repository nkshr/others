(require 'linum)
(global-linum-mode)
(add-hook 'c++-mode-hook '(lambda () (c-set-style "bsd")))

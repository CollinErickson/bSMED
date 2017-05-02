## ------------------------------------------------------------------------
set.seed(0) # 0 for text in here, 3 to get missed peak

## ------------------------------------------------------------------------
quad_peaks_slant <- TestFunctions::add_linear_terms(function(XX) {.2+.015*TestFunctions::add_zoom(TestFunctions::rastrigin, scale_low = c(.4,.4), scale_high = c(.6,.6))(XX)^.9}, coeffs = c(.02,.01))
cf::cf(quad_peaks_slant)


## ------------------------------------------------------------------------
a <- bSMED::bSMED$new(D=2,func=quad_peaks_slant,
                      obj="func", b=3, nb=5,
                      X0=lhs::maximinLHS(20,2),
                      Xopts=lhs::maximinLHS(500,2)
                      )


## ----run_a, fig.width=13, fig.height=9-----------------------------------
a$run(1, noplot=F)

## ----run a second iter, fig.width=13, fig.height=9-----------------------
a$run(1)

## ----run a third iter, fig.width=13, fig.height=9------------------------
a$run(1)

## ----run a fourth iter, fig.width=13, fig.height=9-----------------------
a$run(1)

## ----run a fifth iter, fig.width=13, fig.height=9------------------------
a$run(1)


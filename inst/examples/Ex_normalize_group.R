## Block partition helpers (factor, l2_blocks counts, or row-index list)
## Documented internal: use ::: outside package code / examples.

l2 <- 10L

## One level per row (length l2)
lmebayesCore:::normalize_group(factor(rep(c("A", "B"), each = 5L)), l2)

## Contiguous counts summing to l2
lmebayesCore:::normalize_group(c(4L, 6L), l2)

## Explicit row-index list
lmebayesCore:::normalize_group(list(A = 1:4, B = 5:10), l2)

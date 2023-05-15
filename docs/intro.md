# Introduction

If no reference is included assembly will be 
completed using [flye](https://github.com/fenderglass/Flye) and polished with 
[medaka](https://www.github.com/nanoporetech/medaka). If a reference is provided
alignment will be done with [mini_align](https://github.com/nanoporetech/pomoxis/blob/master/scripts/mini_align)
and variant called using medaka. The workflow has a few optional extras. It can run
[prokka](https://github.com/tseemann/prokka) to annotate the resulting
consensus sequence or [ResFinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) to check it against a database of antimicrobial resistance genes.

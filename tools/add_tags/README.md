# add_tags.py script

Author(s): Miguel A. Caro

## Documentation

**add_tags.py** adds regularization tags to the XYZ file used with **gap_fit** for fine-tuned
control of regularization parameters. It can also sparsify over force observables my masking
individual (vector, i.e., all 3 components) forces. This is useful in case not all forces can
be stored into the available memory or for faster fitting during prototyping (since in typical
datasets forces make up most of the available observables).

IMPORTANT NOTICE: This script is not fool-proof and it requires the user to edit a bunch of
options at the top and knowing what they're doing. Do not use this script to preprocess your
database unless you're familiar with the concept of regularization in GAPs.

== Dependencies ==

* ASE
* Numpy

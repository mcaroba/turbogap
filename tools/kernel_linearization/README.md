# Kernel linearization

Author(s): Patricia Hernandez Leon


## Documentation

Related files included in TurboGAP:

- src/kernel_linearization.f90
    **kernel_linearization** module containing all subroutines and functions related to kernel
    linearization, both for calculations and testing purposes.

- tools/kernel_linearization/get_linearized_sparse_model.f90
    Given a trained potential, **get_linearized_sparse_model.f90** script checks whether kernel
    linearization is favourable and returns the necessary files to perform it with TurboGAP.
    Current implementation supports only energies (forces are on the way).
    It uses the **kernel_linearization** module included in TurboGAP.

- tools/kernel_linearization/test_kernel_linearization_module.f90
- tools/kernel_linearization/test_kernel_linearization_workflow.sh

## Using kernel linearization with TurboGAP

(1) Get the linearized version of the sparse contribution by running **get_linearized_sparse_model.f90**
    script (see **Running the tool** section). 
    The output file(s) are placed on the same folder as the given trained GAP: 
    [alphas file]\_linear.dat (energies), [descriptor file]\_linear (forces)

(2) Activate kernel linearization support when running TurboGAP by adding the following line to the 
    input file:

    kernel_linearization = .true.


## Running the tool 
### (a.k.a. getting the linearized sparse contribution)

You need to have a **gap_files** directory with your GAP potential, check first tools/quip_xml_to_gap
otherwise.

Compile the following code from the directory containing your gap_files folder:

    gfortran -c /turbogap/source/directory/src/kernel_linearization.f90

    gfortran -o get_linearized_sparse_model kernel_linearization.o /turbogap/source/directory/tools/kernel_linearization/get_linearized_sparse_model.f90 -lblas

Note that `/turbogap/source/directory` is the directory where you have your **TurboGAP** installation.
To get the linearized sparse contribution, just run:

    ./get_linearized_sparse_model [path to your GAP file]


## Example/tests



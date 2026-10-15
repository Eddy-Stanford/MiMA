"""mimadoc: MiMA's diagnostics and namelist reference, generated from the
Fortran sources by a static parser (Python >= 3.7, standard library only).

    python3 tools/mimadoc generate        # rewrite the generated parts of docs/
    python3 tools/mimadoc check           # CI: exit 1 if they are out of date
    python3 tools/mimadoc validate diag_table [--nml input.nml]
    python3 tools/mimadoc list [--csv F] [--json F]

Modules:
    fortran      free-form Fortran: statements, scopes, declarations, doc comments
    evaluate     string/integer expressions of the register_diag_field arguments
    conditions   the conditions of a registration: parsing, printing, evaluation
    diagnostics  the registrations, their conditions and send status
    namelists    namelist groups, variables, defaults; the input.nml reader
    diag_table   diag_table validation
    render       Markdown tables
    markers      the generated regions of the Markdown pages
"""

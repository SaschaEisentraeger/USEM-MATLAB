Stand-ins for MATLAB functionality that GNU Octave 8 lacks. cyl_setup_paths.m
adds this folder to the path only when the study runs in Octave.

  contains.m       MATLAB's contains (used by the ABAQUS readers)
  decomposition.m  only d\B is used in ShapeFunctions2D_Metric
  unique.m         Octave's unique does not return the third output of
                   unique(...,'rows','stable') (used by BC_Edge_sub); the
                   shim calls Octave's own unique and adds that output

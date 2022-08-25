1. check on project files

   (a) find . -name ".*.lock" -print

       Remove any locks found.
       Verify that a valid job is not running before removing a lock file.

   (b) proj_search has output.show.license

       Set any found with true to false.

2. To set up a project:

   If derived comparisons are desired, then run with test.create.cases set to true.
   Then:
      mpirun -np OpenParEM2D *.proj
      mv *_prototype_test_cases.csv *_derived_test_cases.csv
      Edit the file *_derived_test_cases.csv as desired.
      Set test.create.cases to false

   To check a setup project for errors:
   process.sh project.proj N

   Adjust *.proj settings and *_derived_test_cases.csv error limits as desired to complete the setup.
   Setup is complete when all tests pass.

3. To check a project for errors:

   process.sh project.proj N

4. To check errors in all projects:

   regression.sh >& regression.log

   from the regression directory.  

   Check "regression.log" for "ERROR", "NOT CONVERGED", and "terminated".
   Check "regression_results.csv" for "FAIL".  Waive small errors.

5. To add projects to the regression suite, add a line to

   regression_case_list.txt




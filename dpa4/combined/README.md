# Exercise 2
I wrote 3 functions for computing MIS:

1. `bibp_mis`: bigger index bigger priority;
2. `bdbp_mis`: bigger degree bigger priority;
3. `luby_mis_random_priority`: random priority.

Note that the size of MIS computed by each method will not be the same.

But the *maximality* and *independence* (No adjacent nodes will be added into MIS at the same time) of the sets can be surely guaranteed as they're base upon generating a Directed Acyclic Graph.

To run the program:
- Build the project;
- Run `<executable-path>/dpa4-nr2 <matrix-market-file-path> > output.txt` in your terminal.
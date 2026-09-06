# Post-Implementation

1. Remove references in the code to engineer-a-plan*.md and tise-task-breakdown.md. Instead, write comments in code to explain how it works.
2. Remove magic numbers from code. Replace with descriptively named constants that are either globally scoped or internal to the local scope depending on what is appropriate in each case. Someone who is not very familiar with C++ should be able to understand the intended behavior of the code after reading the comments. Similarly, someone who is not an expert in Quantum Mechanics (but does possess an undergraduate-level understanding of math and quantum mechanics) should be able to follow along with the math and reasoning for each step based on the comments.

As part of this task, make backup scale (in `const Real scale = std::max(std::abs(reference), 1.0);`) configurable (but still default to 1.0) in classifyAsymptote and related functions which also use it (if any). Any other configuration item in our functions like knumSamples and kRatio should similarly be an optional argument to the function in question.

For the new configurability, we should first make these values configurable via a function input (perhaps with a default value based on constants defined at the top of the file or in the header) and then consider whether these inputs should be made available to the user via the existing interface or remain only configurable in our API as function inputs. Part of this task is to create a document under docs/planning that discusses whether this new configurability should be exposed to the user or not as well as ideas on how this can be done and what parts of the code would need to change for that to be exposed to the UI.

3. Review PR bot comments to see if changes need to be made to address potential bugs.
4. Update TISE README accordingly based on changes we have made to the code as part of the Engineer A tasks.
5. ~~Review engineer B changes, merge into main, rebase the TISE-Generalization branch and resolve any conflicts.~~
6. ~~Create plan for TDSE for engineers A and B.~~
7. Wire up interface with recent changes to TISE and generate end-to-end/integration tests.
8. Verify node placement implementation is at least as capable as what we discussed in original design meetings.
9. Determine what was deemed out of scope (if anything) for each task in the engineer A plan and make both an ADR for each item and a planning doc to address it. 
10. Make updates to the SDD and other docs that need to be made after these changes, if any.
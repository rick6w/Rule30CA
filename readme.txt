Source code for testing random number generation using Wolfram's rule 30 cellular automaton.
rule30.c contains code for experimentation.
test.c is used to test rule 30 with the gcd, gorilla, and birthday test. Compile with "gcc test.c -lm".
note: The tests are written for 32-bit numbers. On the machine used to test the code, an unsigned int was used. Depending on the machine, an unsigned int may have a different number of bits, causing the tests to malfunction.
A data type such as unit_32 would be a better datatype to use were the code to be modified.

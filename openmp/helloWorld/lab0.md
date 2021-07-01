Lab 0

- read the Unit 1 of "A Hands-on Introduction to OpenMP" by Tim Mattson
- implement, compile and run your parallel hello-world program with OpenMP
- answer the following questions
- submit your answers with your program to TA via emails
- you have 24 hours to complete this lab

- how many threads are created? explain with evidences.
- does your program generate deterministic outputs? explain why this happens.
- what are the minimum and the maximum possible outputs for the program below (assuming there are 2 threads)? explain when the minimum and the maximum outputs are generated.

    int count = 0;
    #pragma omp parallel 
    { for (int i = 0; i <10000; i++) {
        count++;
      }
    }
